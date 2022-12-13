/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace MPFR routines
 *
 * Copyright 2022 Andreas Weichselbaum
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * --------------------------------------------------------------------- */

#ifndef __WB_MPFR_HH__
#define __WB_MPFR_HH__

/* -------------------------------------------------------------------- *
   wrapper class for MPFR (multiple precision floating point library)
   Wb,Feb01,12

------------------------------------------------------------------------
 Docu: see $SRC/MPFR/doc/mpfr.pdf
------------------------------------------------------------------------

 MPFR basics
  * before assigning, MPFR variables need to be initialized.
  * when done with a variable, it needs to be cleared out.
  * variable will have the same allocated space during all its life.

  * MPFR functions may create intermediate caches
    => use mpfr_free_cache() to free
      It is strongly advised to do that before terminating a thread
      to avoid memory leaks.

  * The semantics of a calculation in MPFR is specified as follows:
    Compute the requested operation exactly (with "infinite accuracy"),
    and round the result to the precision of the destination variable,
    with the given rounding mode.

  * if two variables are used to store only a few significant bits,
    and their product is stored in a variable with large precision,
    then MPFR will still compute the result with full precision.

------------------------------------------------------------------------
 General rules:
  * all MPFR functions expect output arguments before input arguments;
    typicall: output var, input vars ..., rounding
  * ok to use the same variable for both input and output in the same expression
  * functions returning int: ternary value
    0: value stored in output is the exact result of given mathematical function.
    pos. (neg): value stored in output is greater (lower) than the exact result.

-------------------------------------------------------------------------
 MPFR rounding modes (similar to GMP rounding modes
 GMP_RNDD, GMP_RNDU, GMP_RNDN, ...):

   MPFR_RNDN: nearest: round to nearest (roundTiesToEven in IEEE 754-2008)
   MPFR_RNDZ: zero   : round toward zero (roundTowardZero in IEEE 754-2008)
   MPFR_RNDU: up     : round toward plus infty (roundTowardPositive in IEEE 754-2008)
   MPFR_RNDD: down   : round toward minus infty (roundTowardNegative in IEEE 754-2008)
   MPFR_RNDA: away 0 : round away from zero.

--------------------------------------------------------------------------
 Initialization routines can clearance

   void mpfr_init(mpfr_t x)
   => initialize x to default precision and set its value to NaN

   void mpfr_inits2(mpfr_t x1, mpfr_t x2, ..., (mpfr_ptr)0)
   => initialize set of variables to default precision
   => must be terminated by null-pointer (whose type must also be mpf_ptr)

   void mpfr_init2(mpfr_t x, mpfr_prec_t prec)
   void mpfr_inits2(mpfr_prec_t prec, mpfr_t x1, mpfr_t x2, ..., (mpfr_ptr)0)
   => same as previous init() functions, but with specified exact precision

   void mpfr_clear(mpfr_t x)
   => free the space occupied by the significand of x.
      Make sure to call this function for all mpfr_t variables
      when you are done with them.

   void mpfr_clears(mpfr_t x1, mpfr_t x2, ..., (mpfr_ptr)0)
   => clear sequence of variables
   => must be terminated by null-pointer (whose type must also be mpf_ptr)

   void mpfr_set_default_prec(mpfr_prec_t prec)
   => if MPFR was compiled with --enable-thread-safe,
      the default precision is local to each thread.

   mpfr_prec_t mpfr_get_default_prec(void)
   => return the current default MPFR precision in bits.

   void mpfr_set_prec (mpfr_t x, mpfr_prec_t prec)
   mpfr_prec_t mpfr_get_prec (mpfr_t x)
   => set/get variable specific precision

 * -------------------------------------------------------------------- */

#define DEF_RND GMP_RNDN

  namespace Wb {

   gmp_randstate_t wb_rstate;

template <unsigned P=128>
class mpfr__ { 

  public:

    mpfr__(double val=0) { 
       mpfr_init2(f,P); init(FL,val);
       #pragma omp atomic
       ++num_vars;
    };

    mpfr__(const char *F, int L, double val) {
       mpfr_init2(f,P); init(F_L,val);
       #pragma omp atomic
       ++num_vars;
    };

    mpfr__(const mpfr__& b) {
       mpfr_init2(f,P); mpfr_set(f,b.f,DEF_RND);
       #pragma omp atomic
       ++num_vars;
    };

   ~mpfr__() {
      mpfr_clear(f); 

     #pragma omp atomic
      --num_vars;

      if (!num_vars) {
         Cleanup(FL);
      }
    };

    void Cleanup(const char *F=NULL, int L=0) const;

    mpfr__& init  (const char *F, int L, double val, char lflag=0);
    mpfr__& init_s(const char *F, int L, const char *s, unsigned base=10);
    mpfr__& init_d(const char *F, int L, double val, char dcheck=1);

    mpfr__& set_pi  () { mpfr_const_pi  (f,DEF_RND); return *this; };
    mpfr__& set_log2() { mpfr_const_log2(f,DEF_RND); return *this; };
    mpfr__& set_e   () { init(FL,1); mpfr_exp(f,f,DEF_RND); return *this; };

    mpfr__& set_nan ()         { mpfr_set_nan(f);    return *this; };
    mpfr__& set_inf (int s=+1) { mpfr_set_inf(f,s);  return *this; };
    mpfr__& set_zero(int s=+1) { mpfr_set_zero(f,s); return *this; };

    int isnan    () const { return mpfr_nan_p(f);     };
    int isinf    () const { return mpfr_inf_p(f);     };
    int iszero   () const { return mpfr_zero_p(f);    };
    int isnumber () const { return mpfr_number_p(f);  };   
    int isnonzero() const { return mpfr_regular_p(f); };   

    int sign() { return mpfr_sgn(f); }; 

    void swap(mpfr__<P> &b) { mpfr_swap(f,b.f); };

    int cmp(const mpfr__        &b  ) const { return mpfr_cmp   (f,b.f); };
    int cmp(const double        &val) const { return mpfr_cmp_d (f,val); };
    int cmp(const          int  &val) const { return mpfr_cmp_si(f,val); };
    int cmp(const          long &val) const { return mpfr_cmp_si(f,val); };
    int cmp(const unsigned      &val) const { return mpfr_cmp_ui(f,val); };
    int cmp(const unsigned long &val) const { return mpfr_cmp_ui(f,val); };

    bool operator> (const mpfr__        &x) const { return (cmp(x)> 0); };
    bool operator> (const double        &x) const { return (cmp(x)> 0); };
    bool operator> (const          int  &x) const { return (cmp(x)> 0); };
    bool operator> (const          long &x) const { return (cmp(x)> 0); };
    bool operator> (const unsigned      &x) const { return (cmp(x)> 0); };
    bool operator> (const unsigned long &x) const { return (cmp(x)> 0); };

    bool operator< (const mpfr__        &x) const { return (cmp(x)< 0); };
    bool operator< (const double        &x) const { return (cmp(x)< 0); };
    bool operator< (const          int  &x) const { return (cmp(x)< 0); };
    bool operator< (const          long &x) const { return (cmp(x)< 0); };
    bool operator< (const unsigned      &x) const { return (cmp(x)< 0); };
    bool operator< (const unsigned long &x) const { return (cmp(x)< 0); };

    bool operator>=(const mpfr__        &x) const { return (cmp(x)>=0); };
    bool operator>=(const double        &x) const { return (cmp(x)>=0); };
    bool operator>=(const          int  &x) const { return (cmp(x)>=0); };
    bool operator>=(const          long &x) const { return (cmp(x)>=0); };
    bool operator>=(const unsigned      &x) const { return (cmp(x)>=0); };
    bool operator>=(const unsigned long &x) const { return (cmp(x)>=0); };

    bool operator<=(const mpfr__        &x) const { return (cmp(x)<=0); };
    bool operator<=(const double        &x) const { return (cmp(x)<=0); };
    bool operator<=(const          int  &x) const { return (cmp(x)<=0); };
    bool operator<=(const          long &x) const { return (cmp(x)<=0); };
    bool operator<=(const unsigned      &x) const { return (cmp(x)<=0); };
    bool operator<=(const unsigned long &x) const { return (cmp(x)<=0); };

    bool operator!=(const mpfr__        &x) const { return (cmp(x)!=0); };
    bool operator!=(const double        &x) const { return (cmp(x)!=0); };
    bool operator!=(const          int  &x) const { return (cmp(x)!=0); };
    bool operator!=(const          long &x) const { return (cmp(x)!=0); };
    bool operator!=(const unsigned      &x) const { return (cmp(x)!=0); };
    bool operator!=(const unsigned long &x) const { return (cmp(x)!=0); };

    bool operator==(const mpfr__        &x) const { return (cmp(x)==0); };
    bool operator==(const double        &x) const { return (cmp(x)==0); };
    bool operator==(const          int  &x) const { return (cmp(x)==0); };
    bool operator==(const          long &x) const { return (cmp(x)==0); };
    bool operator==(const unsigned      &x) const { return (cmp(x)==0); };
    bool operator==(const unsigned long &x) const { return (cmp(x)==0); };

    mpfr__& operator= (double val) { return init(FL,val); };
    mpfr__& operator= (const mpfr__& b) {
       mpfr_set(f,b.f,DEF_RND);
       return *this;
    };

    mpfr__& times(const mpfr__& b, mpfr__ &c) const { 
       mpfr_mul(c.f,f,b.f,DEF_RND); return c; };
    mpfr__& divby(const mpfr__& b, mpfr__ &c) const { 
       mpfr_div(c.f,f,b.f,DEF_RND); return c; };

    mpfr__& operator*= (const mpfr__& b) {
       mpfr_mul(f,f,b.f,DEF_RND); return *this; };
    mpfr__& operator/= (const mpfr__& b) {
       mpfr_div(f,f,b.f,DEF_RND); return *this; };
    mpfr__& operator+= (const mpfr__& b) {
       mpfr_add(f,f,b.f,DEF_RND); return *this; };
    mpfr__& operator-= (const mpfr__& b) {
       mpfr_sub(f,f,b.f,DEF_RND); return *this; };

    mpfr__& operator*= (double val) { mpfr__ b(FL,val); return (*this)*=b; };
    mpfr__& operator/= (double val) { mpfr__ b(FL,val); return (*this)/=b; };
    mpfr__& operator+= (double val) { mpfr__ b(FL,val); return (*this)+=b; };
    mpfr__& operator-= (double val) { mpfr__ b(FL,val); return (*this)-=b; };

    mpfr__ operator* (const mpfr__& b) const {
       mpfr__ x; mpfr_mul(x.f,f,b.f,DEF_RND); return x; };
    mpfr__ operator/ (const mpfr__& b) const {
       mpfr__ x; mpfr_div(x.f,f,b.f,DEF_RND); return x; };
    mpfr__ operator+ (const mpfr__& b) const {
       mpfr__ x; mpfr_add(x.f,f,b.f,DEF_RND); return x; };
    mpfr__ operator- (const mpfr__& b) const {
       mpfr__ x; mpfr_sub(x.f,f,b.f,DEF_RND); return x; };

    mpfr__ operator* (int n) const {
       mpfr__ x; mpfr_mul_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator/ (int n) const {
       mpfr__ x; mpfr_div_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator+ (int n) const {
       mpfr__ x; mpfr_add_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator- (int n) const {
       mpfr__ x; mpfr_sub_si(x.f,f,n,DEF_RND); return x; };

    mpfr__ operator* (unsigned n) const {
       mpfr__ x; mpfr_mul_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator/ (unsigned n) const {
       mpfr__ x; mpfr_div_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator+ (unsigned n) const {
       mpfr__ x; mpfr_add_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator- (unsigned n) const {
       mpfr__ x; mpfr_sub_ui(x.f,f,n,DEF_RND); return x; };

    mpfr__ operator* (long n) const {
       mpfr__ x; mpfr_mul_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator/ (long n) const {
       mpfr__ x; mpfr_div_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator+ (long n) const {
       mpfr__ x; mpfr_add_si(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator- (long n) const {
       mpfr__ x; mpfr_sub_si(x.f,f,n,DEF_RND); return x; };

    mpfr__ operator* (unsigned long n) const {
        mpfr__ x; mpfr_mul_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator/ (unsigned long n) const {
        mpfr__ x; mpfr_div_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator+ (unsigned long n) const {
        mpfr__ x; mpfr_add_ui(x.f,f,n,DEF_RND); return x; };
    mpfr__ operator- (unsigned long n) const {
        mpfr__ x; mpfr_sub_ui(x.f,f,n,DEF_RND); return x; };

    mpfr__& Sqrt() { mpfr_sqrt(f,f,DEF_RND); return *this; };
    mpfr__& Exp () { mpfr_exp (f,f,DEF_RND); return *this; };
    mpfr__& Log () { mpfr_log (f,f,DEF_RND); return *this; };
    mpfr__& Log2() { mpfr_log2(f,f,DEF_RND); return *this; };
    mpfr__& Atan() { mpfr_atan(f,f,DEF_RND); return *this; };

    mpfr__& Tanh() { mpfr_tanh(f,f,DEF_RND); return *this; };
    mpfr__& Atanh() { mpfr_atanh(f,f,DEF_RND); return *this; };

    mpfr__& Rand() {
       mpfr_urandomb(f,wb_rstate); 
       return *this;
    };

    mpfr__& Log10() { mpfr_log10(f,f,DEF_RND); return *this; };
    mpfr__& Round() { mpfr_round(f,f); return *this; };

    mpfr__ sqrt() const {
       mpfr__<P> x(*this); mpfr_sqrt(x.f,x.f,DEF_RND); return x; };
    mpfr__ exp() const {
       mpfr__<P> x(*this); mpfr_exp (x.f,x.f,DEF_RND); return x; };
    mpfr__ log() const {
       mpfr__<P> x(*this); mpfr_log (x.f,x.f,DEF_RND); return x; };
    mpfr__ log2() const {
       mpfr__<P> x(*this); mpfr_log2(x.f,x.f,DEF_RND); return x; };
    mpfr__ atan() const {
       mpfr__<P> x(*this); mpfr_atan(x.f,x.f,DEF_RND); return x; };
    mpfr__ log10() const {
       mpfr__<P> x(*this);mpfr_log10(x.f,x.f,DEF_RND); return x; };
    mpfr__ round() const {
       mpfr__<P> x(*this);mpfr_round(x.f,x.f); return x; };

    explicit operator double() const { return mpfr_get_d(f,DEF_RND); };

    explicit operator wbcomplex() const;

    explicit operator bool() const {
       return (mpfr_zero_p(f) ? 0 : 1); 
    };

    double getval_d() const { return mpfr_get_d(f,DEF_RND); };

    bool operator< (double val) const { return (getval_d()< val); };
    bool operator<=(double val) const { return (getval_d()<=val); };
    bool operator> (double val) const { return (getval_d()> val); };
    bool operator>=(double val) const { return (getval_d()>=val); };

    int prec(unsigned base=2) const {
       int p=mpfr_get_prec(f);
         if (base!=2) { p = p*(::log(2)/::log(base)) + 1; }
       return p;
    };

    mpfr__ eps() {  
       mpfr__ x; 
       mpfr_set_uj_2exp(x.f, 1, -prec(2), DEF_RND);
       return x;
    };

    mpfr__& FlipSign() { 
       mpfr_setsign(f,f, mpfr_signbit(f) ? 0 : 1, DEF_RND);
       return *this;
    };

    mpfr__& Abs() {
       mpfr_setsign(f,f, 0, DEF_RND);
       return *this;
    };

    mpfr__ abs()  const { return mpfr__(*this).Abs(); };
    mpfr__ fabs() const { return mpfr__(*this).Abs(); };
    mpfr__ abs2() const {
       mpfr__ x(*this); x*=(*this);
       return x;
    };

    void print(const char *vname="", const char *istr=0, char sep=0) const {
       if (sep) printf("---------"
         "0----+----1----+----2----+----3----+----4----+----"
         "5----+----6----+----7----\n"
       );
       if (vname && vname[0]) printf("%6s = ",vname);
       else printf("         ");

       mpfr_out_str (stdout,
         10, 
          0, 
          f,
          GMP_RNDU 
       );

       if (!istr || !istr[0])  {
          long p,q; unsigned niter; double r, x=getval_d();
          int e=Rational(x,p,q,&r,&niter,NULL,6,1000000,1E-10,1E-12);
          if (!e && p) { char s[32];
             if (q!=1) snprintf(s,32," = %ld/%ld",p,q);
               else snprintf(s,32," = %ld",p);
             istr=s;
          }
       }
       printf(" %s\n", istr && istr[0] ? istr : "");
    };

    size_t to_str(char *s, size_t n_, int base=10) const;
    wbstring to_str(int base=10) const;

    void toStr(
       char *s, unsigned n, int base=10, int nd=-1, char align=0) const;

    wbstring& toStr(
       wbstring &sout,int base=10, int nd=-1, char align=0) const;

    mpfr_t f;

    static unsigned long num_vars; 

  private:

};

   typedef mpfr__<128> quad;   
   typedef mpfr__<128> quad2;  

   typedef mpfr__<192> quad3;  
   typedef mpfr__<256> quad4;  
   typedef mpfr__<256> qquad;  

   template<>
   unsigned long quad::num_vars=0;

   quad mpfr_dummy_quad;

template <unsigned P>
mpfr__<P> operator+(const mpfr__<P> &a) { return mpfr__<P>(a); };

template <unsigned P>
mpfr__<P> operator-(const mpfr__<P> &a) { return mpfr__<P>(a).FlipSign(); };

template <class T, unsigned P>
mpfr__<P> operator*(const T& x_, const mpfr__<P> &a) {
   mpfr__<P> x(x_); return (x*=a);
};

template <class T, unsigned P>
mpfr__<P> operator+(const T& x_, const mpfr__<P> &a) {
   mpfr__<P> x(x_); return (x+=a);
};

template <class T, unsigned P>
mpfr__<P> operator-(const T& x_, const mpfr__<P> &a) {
   mpfr__<P> x(x_); return (x-=a);
};

template <class T, unsigned P>
mpfr__<P> operator/(const T& x_, const mpfr__<P> &a) {
   mpfr__<P> x(x_); return (x/=a);
};

template <unsigned P> inline
int isnan(const mpfr__<P> &a) { return a.isnan(); };

template <unsigned P> inline
int isinf(const mpfr__<P> &a) { return a.isinf(); };

template <unsigned P> inline
int iszero(const mpfr__<P> &a) { return a.iszero(); };

template <unsigned P> inline
int isnumber(const mpfr__<P> &a) { return a.isnumber(); };

template <unsigned P> inline
int isfinite(const mpfr__<P> &a) { return a.isfinite(); };

template <unsigned P>
mpfr__<P> fabs (const mpfr__<P> &x) { return x.fabs(); };

template <unsigned P>
mpfr__<P> sqrt (const mpfr__<P> &x) { return x.sqrt(); };

template <unsigned P>
mpfr__<P> exp  (const mpfr__<P> &x) { return x.exp();  }

template <unsigned P>
mpfr__<P> log  (const mpfr__<P> &x) { return x.log();  }

template <unsigned P>
mpfr__<P> log2 (const mpfr__<P> &x) { return x.log2(); }

template <unsigned P>
mpfr__<P> atan (const mpfr__<P> &x) { return x.atan(); }

template <unsigned P>
mpfr__<P> log10(const mpfr__<P> &x) { return x.log10(); }

template <unsigned P>
mpfr__<P> round(const mpfr__<P> &x) { return x.round(); }

class cc_quad {

 public:

   cc_quad(long double l=0) : f(l) { if (l!=0) {
      long double r1=std::sqrt(l), r2=std::sqrt(f);

      printf("1> %22.20Lf => %22.20Lf (%.3Lg)\n",l,r1,r1*r1-l);
      printf("2> %22.20Lf => %22.20Lf (%.3Lg)\n", f,r2,r2*r2-f);

      r2=std::sqrt(f);
      printf("3> %22.20Lf => %22.20Lf (%.3Lg)\n\n",f,r2,r2*r2-f);

   }};

   cc_quad sqrt() const { cc_quad x(*this); x.f=std::sqrt(f); return x; };
   cc_quad fabs() const { cc_quad x(*this); x.f=std::fabs(f); return x; };
   cc_quad  abs() const { cc_quad x(*this); x.f=std::fabs(f); return x; };

   cc_quad operator+(const cc_quad& b) const { return (f+b.f); };
   cc_quad operator-(const cc_quad& b) const { return (f-b.f); };
   cc_quad operator*(const cc_quad& b) const { return (f*b.f); };
   cc_quad operator/(const cc_quad& b) const { return (f/b.f); };

   cc_quad operator+=(const cc_quad& b) { return (f+=b.f); };
   cc_quad operator-=(const cc_quad& b) { return (f-=b.f); };
   cc_quad operator*=(const cc_quad& b) { return (f*=b.f); };
   cc_quad operator/=(const cc_quad& b) { return (f/=b.f); };

   long double f;

 private:

};

   cc_quad sqrt(const cc_quad &a) { return a.sqrt(); }

}; 

#endif

