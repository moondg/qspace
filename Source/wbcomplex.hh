#ifndef __WB_COMPLEX_HCC__
#define __WB_COMPLEX_HCC__

class wbcomplex;

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
// complex number class
// Wb,Dec20,01

class wbcomplex {

  public:

    wbcomplex() : r(0.), i(0.) {};

    wbcomplex(const wbcomplex &z) = default; 

    wbcomplex(double r0, double i0=0.) : r(r0), i(i0) {};

   ~wbcomplex() = default; 

    wbcomplex& operator= (const wbcomplex &z) = default;

    wbcomplex& operator= (const double &dbl) {
       r=dbl; i=0.;  return *this;
    };

    explicit operator bool() const { return (r || i); };
    bool operator! () const { return (!r && !i); }

    operator double() const { 
       if (i!=0) wblog(FL,
          "ERR got double(%.4g%+.4gi) type conversion !?",r,i);
       return r;
    };

    wbcomplex& init(double r0, double i0) { r=r0; i=i0; return *this; };
    wbcomplex&  set(double r0, double i0) { r=r0; i=i0; return *this; };

    wbcomplex& setRand(double fac=2., double shift=-1.);

    wbcomplex& init_expi(double ph) {  
       r=std::cos(ph); i=std::sin(ph); 
       return *this;
    };

    mxComplexDouble& explicit_copy2(mxComplexDouble &zm) const { 
       zm.real=r; zm.imag=i;
       return zm;
    };

    bool isinf() const { return (std::isinf(r) || std::isinf(i)); };
    bool isnan() const { return (std::isnan(r) || std::isnan(i)); };

    bool isfinite() const { return (std::isfinite(r) && std::isfinite(i)); }
    bool isnormal() const { return (std::isnormal(r) && std::isnormal(i)); }

    bool operator== (const wbcomplex &z) const {
       return (float(r)==float(z.r) && float(i)==float(z.i));
    };
    bool operator!= (const wbcomplex &z) const {
       return (float(r)!=float(z.r) || float(i)!=float(z.i));
    };

    bool operator== (const double &dbl) const {
       return (i==0. && float(r)==float(dbl)); };
    bool operator!= (const double &dbl) const {
       return (i!=0. || float(r)!=float(dbl)); };

    bool operator> (const wbcomplex &z) const {
       return (r>z.r || (r==z.r && i> z.i) ? 1 : 0); };

    bool operator>=(const wbcomplex &z) const {
       return (r>z.r || (r==z.r && i>=z.i) ? 1 : 0); };

    bool operator< (const wbcomplex &z) const {
       return (r<z.r || (r==z.r && i< z.i) ? 1 : 0); };

    bool operator<=(const wbcomplex &z) const {
       return (r<z.r || (r==z.r && i<=z.i) ? 1 : 0); };

    wbcomplex& operator+= (const wbcomplex &z) {
        r+=z.r; i+=z.i;   return *this; };
    wbcomplex& operator-= (const wbcomplex &z) {
        r-=z.r; i-=z.i;   return *this; };
    wbcomplex& operator*= (const wbcomplex &z) { 
        double r_; 
        r_= r*z.r - i*z.i;
        i = r*z.i + i*z.r; r=r_; return *this;
    };
    wbcomplex& operator/= (const wbcomplex &z); 

    wbcomplex& operator*= (double dbl) { r*=dbl; i*=dbl; return *this; }
    wbcomplex& operator+= (double dbl) { r+=dbl; return *this; }
    wbcomplex& operator-= (double dbl) { r-=dbl; return *this; }

    wbcomplex operator+  (const wbcomplex &z) const {
       return wbcomplex(r+z.r, i+z.i); };
    wbcomplex operator-  (const wbcomplex &z) const {
       return wbcomplex(r-z.r, i-z.i); };
    wbcomplex operator*  (const wbcomplex &z) const {
       return wbcomplex( r*z.r-i*z.i, i*z.r+r*z.i); };
    wbcomplex operator/  (const wbcomplex &z) const; 

    wbcomplex operator+  (const double &dbl) const {
       return wbcomplex(r+dbl,i); };
    wbcomplex operator-  (const double &dbl) const {
       return wbcomplex(r-dbl,i); };
    wbcomplex operator*  (const double &dbl) const {
       return wbcomplex(r*dbl, i*dbl); };
    wbcomplex operator/  (const double &dbl) const {
       if (dbl==0.) wblog(FL,"ERR DIV/0! (%g+%gi)/%g",r,i,dbl);
       return wbcomplex(r/dbl, i/dbl);
    };

    wbcomplex operator/= (double dbl) { 
       if (dbl==0.) wblog(FL,"ERR DIV/0! (%g+%gi)/%g",r,i,dbl);
       r/=dbl; i/=dbl; return *this;
    };

    wbcomplex times(const wbcomplex &) const;

    wbcomplex ConjProd(const wbcomplex &z) const {
       return wbcomplex(r*z.r+i*z.i, -i*z.r+r*z.i);  
    };

    void Conj() { i=-i; };
    wbcomplex conj() const { return wbcomplex(r,-i); };

    double abs2() const { return double(r*r+i*i); };
    double abs() const;

    wbcomplex sqrt() const;

    wbcomplex sign() const {
       wbcomplex z; 
       double a=abs(); if (a) { z.init(r/a, i/a); }
       return z;
    };

    double arg() const; 
    double arg_degree() { return (180./M_PI)*arg(); }

    wbcomplex round() const {
       return wbcomplex( ::round(r), ::round(i) );
    };

    wbstring toStr(const char *fmt="%.4g") const;
    void put (const char *vname="ans", const char *ws="base") const;

    double r,i;

  protected:

  private:
};

inline wbcomplex wbcomplex::operator/ (const wbcomplex &b) const {
   double x, den;

   if (fabs(b.r)>=fabs(b.i)) {
      if (b.r==0) { wblog(FL,"WRN DIV/0"); return wbcomplex(Inf,Inf); }
      x=b.i/b.r;
      den = b.r + x*b.i;
      return wbcomplex((r+x*i)/den, (i-x*r)/den);
   }
   else {
      x=b.r/b.i;
      den = b.i + x*b.r;

      return wbcomplex((x*r+i)/den, (x*i-r)/den);
   }
};

inline wbcomplex& wbcomplex::operator/= (const wbcomplex &z) {

   double r_, dbl=1./z.abs2();
   if (z==0.) wblog(FL,"ERR DIV/0! (%g+i%g)/(%g+i%g)",r,i,z.r,z.i);

   r_= (r * z.r + i * z.i) * dbl;
   i = (i * z.r - r * z.i) * dbl; r = r_; return *this;
};

inline wbcomplex wbcomplex::times(const wbcomplex &z) const {
   double a,b; a=r*z.r; b=i*z.i;
   wbcomplex zz(a-b, (r+i)*(z.r+z.i)-a-b); 
   return zz;
};

inline wbcomplex& wbcomplex::setRand(double fac, double shift) {
   static char first_call=1;

   if (first_call) { wb_srand(); first_call=0; } 
   fac/=(double)RAND_MAX;

   if (shift==0.)
        { r=fac*rand();       i=fac*rand();       }
   else { r=fac*rand()+shift; i=fac*rand()+shift; }

   return *this;
};

inline double wbcomplex::abs() const {
   double a=0.; 
   double x=fabs(r), y=fabs(i);

   if (!x) { a=y; } else
   if (!y) { a=x; } else if (x>y)
        { a=y/x; a = x*std::sqrt(1.+a*a); } 
   else { a=x/y; a = y*std::sqrt(1.+a*a); } 

   return a;
};

inline wbcomplex wbcomplex::sqrt() const {
    wbcomplex z; 
    if (!r && !i) { return z; }

    double ar=fabs(r), ai=fabs(i), w, x;
    if (ar>=ai)
         { x=ai/ar; w=std::sqrt(ar)*std::sqrt(0.5*(1.+std::sqrt(1+x*x))); }
    else { x=ar/ai; w=std::sqrt(ai)*std::sqrt(0.5*(x +std::sqrt(1+x*x))); }

    if (r>=0.) { z.init( w, i/(2.*w) ); } else
    if (i>=0.) { z.init( i/(2.*w), w ); } else { z.init(-i/(2.*w),-w ); }

    return z;
};

inline double wbcomplex::arg() const { 
    double phi=0.; 

    if (!i) { if (r<0) { phi=+M_PI; }; } else  
    if (!r) { phi=( i>0? +0.5 : -0.5 )*M_PI; } 
    else {
       phi=atan(i/r);
       if (r<0) { phi+=(i>0 ? M_PI : -M_PI); }
    }
    return phi;
};

void wbcomplex::put(const char *vname, const char *ws) const {

    mxArray *a=NULL;
    if (i)
         { a=Mx::Array<wbcomplex>(1,1).copyFromTR(this); }
    else { a=Mx::Array<double   >(1,1).copyFromTR(&r  ); }

    int e=mexPutVariable(ws,vname,a);
    if (e==1) wblog(FL,
       "ERR failed to put variable %s to workspace `%s'",vname,ws);

    mxDestroyArray(a);
};

inline wbcomplex operator+(const wbcomplex &z) { return z; }
inline wbcomplex operator-(const wbcomplex &z) {
   return wbcomplex(-z.r, -z.i);
};
inline wbcomplex operator*(double d, const wbcomplex &z) {
   return wbcomplex(z.r*d, z.i*d);
};
inline wbcomplex operator+(double d, const wbcomplex &z) {
   return wbcomplex(d+z.r, +z.i);
};
inline wbcomplex operator-(double d, const wbcomplex &z) {
   return wbcomplex(d-z.r, -z.i);
};
inline wbcomplex operator/(double d, const wbcomplex &z) {
   return wbcomplex(d,0)/z;
};

inline bool isnan(const wbcomplex &z) { return z.isnan(); };

inline bool isfinite(const wbcomplex &z) { return z.isfinite(); };
inline bool isnormal(const wbcomplex &z) { return z.isnormal(); };

template <class T> inline T real(const T &d) { return d; };
template <class T> inline T imag(const T &d) { return 0; };
template <class T> inline T conj(const T &d) { return d; };

template <class T> inline T CONJ(const T &d) { return d; };
template <class T> inline T ABS (const T &a) { return (a<0 ? -a : a); };
template <class T> inline T ABS2(const T &d) { return d*d; };
template <class T> inline T NORM (const T &a){ return (a<0 ? -a : a); };
template <class T> inline T NORM2(const T &d){ return d*d; };

inline double real(const wbcomplex &z) { return z.r; };
inline double imag(const wbcomplex &z) { return z.i; };
inline wbcomplex conj(const wbcomplex &z) {
   return wbcomplex(z.r,-z.i);
};

inline double ROUND(const double &d) { return ::round(d); }
inline wbcomplex ROUND(const wbcomplex &d) { return d.round(); }

#ifdef __WB_MPFR_HH__
template <unsigned P> inline
Wb::mpfr__<P> ROUND(const Wb::mpfr__<P> &q) { return q.round(); };
#endif

template <class T> inline
T SQRT(const T &d) { return std::sqrt(d); }

template <> inline
wbcomplex SQRT(const wbcomplex &d) { return d.sqrt(); }

#ifdef __WB_MPFR_HH__
template <unsigned P> inline
Wb::mpfr__<P> SQRT(const Wb::mpfr__<P> &q) { return q.sqrt(); };
#endif

template <> inline
double ABS(const double &a) { return std::fabs(a); };
template <> inline
double NORM(const double &a) { return std::fabs(a); };

template <> inline
wbcomplex CONJ(const wbcomplex &a) { return a.conj(); };

inline double ABS  (const wbcomplex &a) { return a.abs();  };
inline double ABS2 (const wbcomplex &a) { return a.abs2(); };
inline double NORM (const wbcomplex &a) { return a.abs();  };
inline double NORM2(const wbcomplex &a) { return a.abs2(); };

inline bool ISINF(const double &a) { return std::isinf(a); }
inline bool ISINF(const wbcomplex &a) {
   return (std::isinf(a.r) || std::isinf(a.i)); }

inline bool ISNAN(const double &a) { return std::isnan(a); }
inline bool ISNAN(const wbcomplex &a) {
   return (std::isnan(a.r) || std::isnan(a.i)); }

inline double REAL(const double &a) { return a; }
inline double REAL(const wbcomplex &a) { return a.r; }

inline double IMAG(const double &a __attribute__ ((unused))) { return 0; }
inline double IMAG(const wbcomplex &a) { return a.i; }

#ifdef __WB_MPFR_HH__

template <unsigned P> inline
Wb::mpfr__<P> ABS(const Wb::mpfr__<P> &a) { return a.abs(); };
template <unsigned P> inline
Wb::mpfr__<P> NORM(const Wb::mpfr__<P> &a) { return a.abs(); };

template <unsigned P> inline
Wb::mpfr__<P> ABS2(const Wb::mpfr__<P> &a) { return a.abs2(); };
template <unsigned P> inline
Wb::mpfr__<P> NORM2(const Wb::mpfr__<P> &a) { return a.abs2(); };

#endif

template <class T> inline double ABSDIFF_F(
   const T &d1, const T &d2
){ return double(fabs(float(d1)-float(d2))); };

template <> inline double ABSDIFF_F(
   const wbcomplex &d1, const wbcomplex &d2
){
   float dr=float(d1.r)-float(d2.r), di=float(d1.i)-float(d2.i);
   return std::sqrt(double(dr*dr+di*di));
};

inline wbcomplex exp(const wbcomplex &z) {
   double r=exp(z.r);
   return wbcomplex(r*cos(z.i), r*sin(z.i));
};

inline wbcomplex expi(const double &x) {
   return wbcomplex(cos(x), sin(x));
};

inline wbstring toStr(const wbcomplex &z, const char* fmt="%.4g");
inline wbstring toStr(const double &d, const char* fmt="%.4g");

#endif

