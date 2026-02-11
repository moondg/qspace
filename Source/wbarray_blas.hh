/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbarray (array class, col-major)
 *
 * Copyright 2024 Andreas Weichselbaum
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

#ifndef __WB_ARRAY_BLAS_HH__
#define __WB_ARRAY_BLAS_HH__

// ----------------------------------------------------------------- //
// also link wbblas.hh for wrapper routines to Fortran library
// extern "C" { ... }
// Wb,Sep16,09
// ----------------------------------------------------------------- //

namespace Wb {

template<class TA, class TB, class TC>  
wbarray<TC>& MatProd(
    const wbarray<TA> &A, const wbarray<TB> &B,   
    wbarray<TC> &C,                               
    char aflag='N', char bflag='N',
    TA afac=1., TC cfac=0., char cforce=0
);

template<class T>  
void VMatProd(
    const wbarray<T> &A0, const wbarray<T> &B0, wbarray<T> &C,
    char aflag='N', char bflag='N', const T afac=1., const T cfac=0.,
    const char cforce=1
);

template<class T>  
void MatVProd(
    const wbarray<T> &A, const wbvector<T> &v, wbarray<T> &x,
    char aflag='N', char bflag='N', const T afac=1., const T cfac=0.,
    const char cforce=1
);

template<class T> T VMatVprod( 
const wbarray<T> &v1, const wbarray<T> &M, const wbarray<T> &v2);

template<class T> T VMatVprod( 
const wbvector<T> &v1, const wbarray<T> &M, const wbvector<T> &v2);

template<class TA, class TB, class TC>
void DMatProd(
    const wbarray<TA> &A0, const wbarray<TB> &B0, wbarray<TC> &C,
    char aflag='N', char bflag='N', const TA afac=1., const TC cfac=0.,
    const char cforce=1
);

template<class TA, class TB, class TC>
void DMatProd(
    const wbvector<TA> &A0, const wbarray<TB> &B0, wbarray<TC> &C,
    char aflag='N', char bflag='N', const TA afac=1., const TC cfac=0.,
    const char cforce=1
){
    wbarray<TA> A; A.init2Vec(A0.data,A0.len,'r');
    DMatProd(A,B0,C,aflag,bflag,afac,cfac,cforce);
};

template<class TA, class TB, class TC>
void DMatProd(
    const wbarray<TA> &A0, const wbvector<TB> &B0, wbarray<TC> &C,
    char aflag='N', char bflag='N', const TA afac=1., const TC cfac=0.,
    const char cforce=1
){
    wbarray<TB> B; B.init2Vec(B0.data,B0.len,'r');
    DMatProd(A0,B,C,aflag,bflag,afac,cfac,cforce);
};

}; 

   template<class T>
   wbarray<T>& wbInverse(
      const char *F, int L,
      const wbarray<T> &M QS_UNUSED_VAR, wbarray<T> &Minv
   ){
      wblog(F_L,"ERR %s() not yet defined for T=%s",FCT,TSTR(T));
      return Minv;
   }

   template<class T>
   wbarray<T>& wbInverse(const char *F, int L, wbarray<T>&M);

   void wbEigenS( 
       const wbarray<double>&M,                      
       wbarray<double>&U, wbvector<double> &E        
   );

   void wbEigenS( 
       const wbarray<wbcomplex>&M,                   
       wbarray<wbcomplex>&U, wbvector<double> &E     
   );

   void wbEigen_CS( 
       wbarray<wbcomplex> &M,                        
       wbarray<wbcomplex>&U, wbvector<wbcomplex> &E, 
       char wjob='R', char issym='H', char tnorm=1, char qflag=0
   );

   template<class T>
   void wbEigen_CS_regen(
      const wbarray<T> &M0, wbarray<T> &V, wbvector<T> &E,
      char wjob
   );

   template<class T> 
   void wbSVD(
       const wbarray<T> &A, 
       wbarray<T> &U, wbvector<double> &S, wbarray<T> &Vc, 
       const WBINDEX &I=WBINDEX(),
       char bare=0
   );

   template <class TQ>
   void wbEigenS (
       const wbsparray<double>&M,
       wbsparray<double>&U, wbvector<double> &E
   ){
       wbarray<double> Mf,Uf; M.toFull(Mf);
       wbEigenS(Mf,Uf,E); U.init(FL,Uf,1e-14);
   };

template<class T>
class opFlags {
  public:

    opFlags() : p(0), type(0) {};
    opFlags(const opFlags &f) : p(f.p), type(0) {};

    opFlags(const char *F, int L, char f)
     : p(0), type(0) { init(F,L,f); };

    opFlags(const opFlags &f1, const opFlags &f2)
     : p(0), type(0) { init(f1,f2); };

    opFlags<T>& operator= (char f) { init(f); return *this; }
    opFlags<T>& operator= (const opFlags<T> &x) { p=x.p; return *this; }

    opFlags(char f) : p(0), type(0) { init(FL,f); }; 
    operator char() { return p; } 

    opFlags<T>& init(const char *F, int L, char f);
    opFlags<T>& init(char f) { return init(0,0,f); };

    opFlags<T>& init(const opFlags<T> &f1, const opFlags<T> &f2){
    p=(f1.p ^ f2.p);  return *this; };

    char trans() const { return (p & TRANS_PAT); }
    char conj () const { return (p & CONJ_PAT); }
    char idle () const { return (!(p & TRANS_PAT) && !(p & CONJ_PAT)); }

    void applyOrRef(
      const char *F, int L,
      const wbarray<T> &A0, wbarray<T> &A, const T fac=1
    ) const;

    template <class T2> 
    T2& applyOrRef(const char *F, int L, const T2 &A0, const T fac=1) const;

    template <class T2> 
    const T2& applyConjOrRef(const char *F, int L, const T2 &A0, const T fac=1) const;

    template <class T2>
    char apply(const char *F, int L, T2 &A, const T fac=1) const;

    char toChar() const {
       char c, tflag=(p & TRANS_PAT), cflag=(p & CONJ_PAT);
       if (tflag)
            c=(cflag ? 'C' : 'T');
       else c=(cflag ? 'c' : 'N');
       return c;
    }

    mxArray *toMx() const { return numtoMx(this->toChar()); };

    void info(const char *F, int L, const char *istr="") {
       wblog(F,L,"<i> opFlags<%s> %s: %d -> [%d %d]",
       TSTR(T), istr && istr[0] ? istr : "",
       p, trans(), conj());
    };

    char p;   

  private:

    T type; 

    static char TRANS_PAT;
    static char CONJ_PAT;
};

template <class T> char opFlags<T>::TRANS_PAT = 1; 
template <class T> char opFlags<T>::CONJ_PAT  = 2; 

template <class T> inline
opFlags<T>& opFlags<T>::init(
   const char *file, int line, char f
){
   p=0;

   if (f>16) {
      if (!strchr("ENCTc",f)) wblog(__FL__,"ERR invalid flag %c<%d>",f,f);
      if (strchr("CT",f)) p |= TRANS_PAT;
      if (f=='E') type=(T)'E'; 

      if (typeid(T)==typeid(wbcomplex))
      if (strchr("Cc",f)) p |= CONJ_PAT;
   }
   else {
      p=f;
      if (p & CONJ_PAT && typeid(T)!=typeid(wbcomplex))
      p &= ~CONJ_PAT;
   }

   return *this;
};

template <class T> inline
void opFlags<T>::applyOrRef(
   const char *F, int L,
   const wbarray<T> &A0, wbarray<T> &A, const T fac
) const {

   unsigned r=A0.SIZE.len;
   char p2=p;

   if (p  && ((p & ~TRANS_PAT) & ~CONJ_PAT)) 
   wblog(F,L,"ERR invalid flags (%d)",p);

   if (!r) { A.init(); return; }
   if (r>2 && r%2) wblog(F,L,"ERR invalid rank-%d object (odd)",r);

   if (r==1) p2 &= ~TRANS_PAT;

   if (p2 || fac!=T(1)) {
      if (p2 & TRANS_PAT) {
          wbperm P; P.initTranspose(r);
          A0.permute(A,P); 
      } else A=A0;

      if (p2 & CONJ_PAT) A.Conj();

      if (fac!=T(1)) A*=fac;
   }
   else A.init2ref(A0); 
};

template <class T>
template <class T2> inline
T2& opFlags<T>::applyOrRef(
  const char *F, int L, const T2 &A, const T fac
) const {

   if (p  && ((p & ~TRANS_PAT) & ~CONJ_PAT)) 
   wblog(F_L,"ERR invalid flags (%d)",p);

   if (!p && fac==1) return A;

   T2 *a = new T2;
   if (p & TRANS_PAT) { A.transpose(F_L,*a); } else { (*a)=A; }
   if (p & CONJ_PAT) a->Conj();
   if (fac!=1) (*a)*=fac;

   return (*a);
};

template <class T>
template <class T2> inline
const T2& opFlags<T>::applyConjOrRef(
  const char *F, int L, const T2 &A, const T fac
) const {

   char c=(p & CONJ_PAT);
   if (p && (p & ~TRANS_PAT & ~CONJ_PAT)) 
      wblog(F_L,"ERR invalid flags (%d)",p);

   if (!c && fac==1) return A;
   else {
      T2 *a = new T2(A);
         if (c) a->Conj();
         if (fac!=1) (*a)*=fac;
      return (*a);
   }
};

template <class T>
template <class T2> inline
char opFlags<T>::apply(
  const char *F, int L, T2 &A, const T fac
) const {

   if (p  && ((p & ~TRANS_PAT) & ~CONJ_PAT)) 
   wblog(F_L,"ERR invalid flags (%d)",p);

   if (p || fac!=1) {
      if (p & TRANS_PAT) A.Transpose(F_L);
      if (p & CONJ_PAT) A.Conj();
      if (fac!=1) A*=fac;
      return 1;
   }
   return 0;
};

#endif

