/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : clebsch (for abelian and non-abelian symmetries)
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

#ifndef __WB_CLEBSCH_GORDAN_HH__
#define __WB_CLEBSCH_GORDAN_HH__

//====================================================================//
// CStore => global CG container: CStore<gTQ> gCS; // CG_GLOBAL
//
//    global container for Clebsch Gordan coefficients
//     - based on CData for single (Q1,Q2->Q) mapping.
//     - based on QVec, qset, cdata below.
//
// cdata : wbarray<double>    // full data arrays!
// cdata : wbsparray<double>  // switching to sparse // Wb,Dec22,11
//
//    actual storage of ClebschGordan coefficients
//     - of rank 3(+1) for standard CG (3+1=4 in case of multiplicity)
//     - cdata also includes (contracted) derivatives of ANY rank
//
// QVec : wbvector<QType>
//
//    single set of types of quantum numbers with QType
//     - the most basic symmetry unit
//
// qset : wbvector<TQ>
//
//    actual set of quantum numbers interpreted through QVec
//
// Wb,Sep20,09 ; Wb,May01,10
//
// ------------------------------------------------------------------ //
// NB! sign convention on CGCs #CD_NORM_QS  #CD_NORM_EXTERN
//
// using CRef, ensure that CData is simply normalized to 1
// internally within gCS. This ensures that for each OM set
//  * overlaps can simply be computed based on cgw data!
//  * OM within gCS can be also determined through cdata.norm2().
// return value: number of trailing dimensions to be interpreted
// as OM space.
//
// Changing cgw from vector to // CGW_ARRAY
//    -> wbMatrix // Wb,Jun08,20
//       -> wbarray  // Wb,Apr23,21
//
//  * 1st (fast) index contracts onto OM index of CData // wdim1()
//  * 2nd (slow) index onto OM index of A.data // wdim2()
//    #CGW_INDEX_ORDER
//
// => normalization convention is cgw^T*cgw=1 // #CGW_NORM_QS
//    such that cgw represents an isometry
//    where cgw may have have fewer (but not more) columns than rows
//    so it may hold cgw*cgw^T != 1
//
// Benefits of this normalization
//  - OM normalization convention of CData   // #CD_NORM_QS
//    required for simple orthonormal OM basis decomposition
//    remains intact also for CData_actual = CData * cgw
//    irrespective of direction of legs
//  - therefore A.data already carries correct 3j weight factors for simple
//    contraction, such that e.g., A.norm2() of some QSpace A
//    is simply carried by |A.data|^2
//
// NB! for `external' purposes only (readability) #CD_NORM_EXTERN
// and only for rank-2 tensors (i.e, scalar operators,
// which never have OM!), use |cgw|=sqrt(qdim)
// to reflect actual matrix elements in A.data, e.g., see
// QSpace::NormCGW() -> CGR->NormSignW -> normExt()!
// @QSpace/display.m also assumes this convention throughout.
//
// ------------------------------------------------------------------ //
// Other (not so useful) attempts to normalize cgw
// (with corresponding inverse factor multiplied onto data{:}):
//
//     |cgw| := sqrt(max(size(cdata))) // former version of CD_NORM_EXTERN
//
// (+) rank-2 tensor are exact identities (or untaries for 1J symbols)
// (+) rank-3 tensors are mostly normalized to max(|CGC|)
//     e.g. (SS|0) is normalized the same way as (S0|S)
//     whereas their standard normalization [see NormStd()]
//     would be sqrt(d_out), here 1 or (2S+1), respectively!
//
// (-) standard rank-3 identity A tensors do not have cgd = Id
// (-) does not work specifically well for rank>3 tensors.
//     e.g. consider
//         contract( (1,1|2), 3, (1,1|2), 3 )   (i)
//     vs. contract( (1,1|0), 3, (1,1|0), 3 )   (ii)
//     while the output CData space is the same (1,1,1,1)
//     in either case, already the normalization of the
//     CData requires *different* normalization factors 1/sqrt(d3)
//     depending on the *contracted* 3rd index containing a multiplet
//     of dimension d3. This information is no longer available
//     once the CGCs have been contracted!
//
// Wb,Oct16,14 ; Wb,Nov04,14
//
// ================================================================== //
// RStore => global constructor class for non-abelian Lie groups
//
// RStore <TQ,TD> gRS; // global
//  * map[QType] :=> genRG_struct buf[]
//
//    genRG_struct<TQ,TD>
//     *  map[qset] :=> genRG_base<TQ,RTD>  RSet[]
//
//       genRG_base<TQ,RTD>
//        * basic symmetry settings such as
//          generators, commutator relations (CR), ...
//
// Wb,Jun07,10
//
// ================================================================== //
/* CHANGE LOG

// Wb,Aug28,20 // PERMIT_RANK1_Q0 ----------------------------------- //

      Permit rank-1 QSpace if scalar, i.e. for q=0 only;
      e.g. for SU(2), can contract (S,0|S) at '13'
      for any spin-S multiplet, resulting in a non-zero weight
      since by CD_NORM_QS with d=2S+1 the multiplet dimension,
      trace('S,0;S,'13') = 1/sqrt(d) * d = sqrt(d)
      [e.g. see @QSpace/getIdentity3.m with --rho flag]
      however, contract('S,S;S,'13') = 0, i.e., rank-1 QSpace is
      not permitted for non-zero (non-scalar) symmetry label q!=0.

      NB! a rank-1 QSpace may be generated, e.g., by contraction
      to obtain expectation values [HAM/calc_SdotS.m] after using
      @QSpace/getIdentity3.m

      NB! contract('0,S,S*,'13','S,S*','12')
        = sqrt(d) * trace('S,0;S,'13')
        = 1.
      since 2nd line acquires additional sqrt(d) factor from 'S,S*'
       => X-symbol has x3=1 - ok!

   => permit rank-1 QSpace as output, but not necessarily as input
      yet with exceptions, such as getDimQS so that at least QSpace
      can be displayed using @QSpace/display.m

// Wb,June 2014 ----------------------------------------------------- //

    * removed CG_SU2 (merged with treatment of QT_SUN)
    * removed zflag from qset [e.g. nq = n(q) + n(qz)]
    * removed zflags from QType and qset
    * removed qset::RemoveZLabels (still got QSpace::RemoveZLabels!)
    * removed following variables from CData as these are redundant
      with data in RStore => genRG_struct => genRG_base:
      - qset<TQ> Q1,Q2,Q;     // (Q1,Q2) -> Q
      - wbMatrix<TQ> Z1,Z2,Z; // (Z1,Z2) -> Z
      use gRS container to retrieve these where required
    * renamed CData => cdata; CGSpace => CData

   => see Archive/clebsch_140627.[ch]* for version including zflag

//====================================================================*/

#ifndef RTD
  #ifdef __WB_MPFR_HH__

    #define MTI unsigned long 
    #define RTD Wb::quad

  #else

    #define MTI unsigned long
    #define RTD double

  #endif
#endif

#ifndef gTQ

   #define gTQ int
#endif

#ifndef gTD
#define gTD double
#endif

#define CDATA_TQ    CData<TQ,RTD>
#define cdata__     cdata<RTD>
#define SPARR_RTD   wbsparray<RTD>

   template <class TQ, class TD> class CData;
   template <class TQ> class CRef;

   class QVec;

   int get_CG_VERBOSE(const char *F=0, int L=0);
   int CG_VERBOSE=get_CG_VERBOSE();

   unsigned get_CG_FIXIT(const char *F=0, int L=0);
   unsigned CG_FIX=get_CG_FIXIT();

   #define cgfix_CID3  1  
   #define cgfix_CTR   2  

   unsigned get_CG_PREVIEW(const char *F=0, int L=0);
   unsigned CG_PREVIEW=get_CG_PREVIEW();

   #define RC_LOAD_ERR_AB_MSG "RC_LOAD_ERROR (A or B)"

   wbstring RC_LOG; 

   #define CG_SKIP_DEPS1 1e-12 
   #define CG_SKIP_DEPS2 1e-14

   #define CG_SKIP_REPS  1e-17 

#ifdef __WB_MPFR_HH__
   #define CG_SKIP_EPS1 1e-16 
   #define CG_SKIP_EPS2 1e-18 

   #define CG_EPS  1e-38      
   #define CG_EPS0 1e-36      
   #define CG_EPS1 1e-16      
   #define CG_EPS2 1e-18      

#else

   #define CG_EPS  2.22e-16    
   #define CG_EPS0 1e-15       

   #define CG_SKIP_EPS1 1e-14  
   #define CG_SKIP_EPS2 1e-16

   #define CG_EPS1 1e-10
   #define CG_EPS2 1e-12

#endif

   double cg_eps1=CG_EPS1;
   double cg_eps2=CG_EPS2;

enum QT_QSPACE {
   QT_UNKNOWN,
      QT_P,   
      QT_ZN,  
      QT_U1,  
      QT_A4,  
      QT_SUN, 
      QT_SpN, 
      QT_SON, 
      QT_SEN, 
   QT_NUM_TYPES
};

#define QT_ABELIAN  QT_U1 
#define QT_RANK0    QT_A4 

const char* QT_STR[QT_NUM_TYPES]=  {
  "QT:-",   
   "P",     
   "ZN",    
   "A",     
   "A4",
   "SUN",   
   "Sp2N",  
   "SO2N",  
   "SO2N_1" 
};

#define CGC_ALL_ABELIAN "A*"

enum CD_TYPE {
     CD_UNKNOWN,    
     CD_ABELIAN,    
     CD_IDENTITY,   

     CD_REF_INIT,   

     CD_BSZ_INIT,   

     CD_FROM_CGR,   

     CD_FROM_CTR,   
     CD_FROM_DEC,   

     CD_STD3,       

     CD_1JSY_ST3,   
     CD_1JSY_GEN,   

     CD_STD3_X,     
     CD_GEN3_X,     
     CD_IMPLICIT,   
     CD_EXPLICIT,   

     CD_OTHER_15,   

     CD_COMPLETE,   
     CD_COMPLETE_U, 
     CD_NP_ZERO,    
  CD_NUM_TYPES_X    
};

#define CD_NUM_TYPES    CD_COMPLETE  

#define CID_RANK1_Q0  60   

const char* CD_TYPE_STR[CD_NUM_TYPES_X+2] = {
   "C:--",   
   "C:iA",   
   "C:Id",   
   "C:iR",   
   "C:iS",   
   "C:cR",   
   "C:cX",   
   "C:cD",   
   "C:s3",   
   "C:r1J",  
   "C:c1J",  
   "C:c3x",  
   "C:cgx",  
   "C:iI",   
   "C:iE",   
   "C?15",   
   "C:#1",   
   "C:f",    
   "C:U",    
   "C:Z",    
"C:#X"};     

enum CR_TYPE {
     CR_DEFAULT,    
     CR_ABELIAN,    

     CR_CTR_SCALAR, 

     CR_CTR_ZERO,   

CR_NUM_TYPES };

const char* CR_TYPE_STR[CR_NUM_TYPES]=  {
    "R:--", 
    "R:iA", 
    "R:cS", 
    "R:c0", 
};

namespace Mx { 

template <>
size_t Mx::Array<CD_TYPE>::copy_to(CD_TYPE* b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

template <>
size_t Mx::Array<CR_TYPE>::copy_to(CR_TYPE* b, const wbperm *P, char tcheck)
const { return rcopy_to(b,P,tcheck); };

};

template <class T> class qset;
template <class T> class CStore;

template <class T> class QMap;
template <class TQ, class TD> class blockSpaceQS;

#include <map>

namespace CG {

class FOM__ : public wbvector<int> { 
public:

   int& me() { int i=my_tid(FL); return data[i]; };
   int tid() { return my_tid(FL); };

   int& operator++()     { int i=my_tid(FL); return ++data[i];     };
   int& operator--()     { int i=my_tid(FL); return --data[i];     };
   int& operator=(int x) { int i=my_tid(FL); return (data[i]=x);   };

  void status(const char *F=nullptr, int L=0, const char *fct=nullptr) {
     wblog(F_L,"TST %s() %s= [ %s ] (len=%d)",
     fct? fct:FCT, fct? "thread_fOM ":"", STR_(this),len);
  };

protected:
private:

   unsigned my_tid(const char *F, int L) { 
      unsigned l=0, tid=Wb::get_omp_tid_nn(&l);

      if (!len) { 
         #pragma omp critical (CGT_FOM__)
         { int n=omp_get_num_threads();
               n=MAX(n,QSP_NUM_THREADS);
               n=MAX(n,OMP_NUM_THREADS); n=MAX(n,int(tid+1));
               n=Wb::pow2_ceil(n); 
           init(n);
         }
      }
      if (unsigned(tid)>=len) wblog(FL, 
        "ERR CG::FOM__::%s() index out of bounds (%d/%d; l=%d)",tid,len,l);
      return tid;
   };
}; 

   FOM__ thread_fOM;

   int MAX_OM=9; 

   template <class TD> inline
   int signFirstVal(const char *F, int L,
      const TD *d, SPIDX_T n,
      double eps1=CG_EPS1, double eps2=CG_EPS2
   );

   template <class TD> inline
   int rangeSignConvention(const char *F, int L,
      TD *d, SPIDX_T n,
      double eps1=CG_EPS1, double eps2=CG_EPS2
   );

   template<class T>
   double FixRational(const char *F, int L,
      T* d, SPIDX_T n, unsigned niter=0,
      double eps1=-1, double eps2=-1
   );
}; 

namespace DY { 

   template<class TQ>
   size_t wdim_A(unsigned n, const TQ *qs); 

   template<class TQ>
   size_t wdim_B(unsigned n, const TQ *qs); 

   template<class TQ>
   size_t wdim_C(unsigned n, const TQ *qs); 

   template<class TQ>
   size_t wdim_D(unsigned n, const TQ *qs); 

};

class QType { 

  public:

    QType(QT_QSPACE t=QT_UNKNOWN, unsigned s=0)
     : type(t), sub(s) { if (t || s) { validType(FL); }};

    QType(const char *s ) : type(QT_UNKNOWN), sub(0) { init(0,0,s); };
    QType(const QType &q) : type(q.type), sub(q.sub) { };

    QType(const char *F, int L, const mxArray *a)
      : type(QT_UNKNOWN), sub(0) { init(F,L,a); };

   ~QType() { init(); };

    QType& init() { type=QT_UNKNOWN; sub=0; return *this; };

    QType& initx(QT_QSPACE t, unsigned m=0) {
       if (!t && !m) { init(); return *this; }
       type=t; sub=m; validType(FL);
       return *this;
    };

    QType& init(QT_QSPACE t, unsigned m=0) {

       if (!t && !m) { init(); return *this; }

       type=t;
       if (type==QT_SpN) { sub=m/2;
          if (m%2) wblog(FL,"ERR %s() invalid sub=%d with SpN",FCT,m);
       }
       else { sub=m; }

       validType(FL); return *this;
    };

    int init_s(const char *s);

    QType& init(const char *s) { return init(FL,s); }

    QType& init(const char *F, int L, const char *s) {
       if (!s || !s[0]) { type=QT_UNKNOWN; sub=0; }
       else {
          if (init_s(s)) wblog(FL,"ERR invalid QType '%s'",s);
          else validType(F_L);
       }
       return *this;
    };

    QType& init(const char *F, int L, const mxArray *a) {
       wbstring s(F_L,a);
       return init(F_L,s.data);
    };

    QType& init(const QType &q) {
       type=q.type; sub=q.sub; return *this;
    };

    void swap(QType &q) { 
       SWAP(type,q.type);
       SWAP(sub, q.sub );
    };

    inline char validType(const char *F=0, int L=0) const;

    bool isUnknown() const { return validType(FL)==0; };
    bool isKnown()   const { return validType(FL)>0; };        
    explicit operator bool() const { return validType(FL)>0; }; 

    bool isAbelian(const char *F=NULL, int L=0) const {
       if (type<=QT_UNKNOWN || type>=QT_NUM_TYPES) wblog(F_L,
          "ERR invalid type %d, `%s'",type,STR(*this));
       return ((type && type<=QT_ABELIAN)? 1 : 0);
    };

    bool isNonAbelian(const char *F=NULL, int L=0, char lflag=0) const {
       if (!lflag && !isUnknown()) { 
       if (type<=QT_UNKNOWN || type>=QT_NUM_TYPES) wblog(F_L,
          "ERR invalid type %d, `%s'",type,STR(*this)); }
       return (type>QT_ABELIAN? 1 : 0);
    };

    bool isSU2() const { return (type==QT_SUN && sub==1); };
    bool isSUN() const { return (type==QT_SUN); };
    bool isSU(unsigned N) const { return (type==QT_SUN && sub+1==N); };

    bool isSpN() const { return (type==QT_SpN); };
    bool isSON() const { return (type==QT_SON); };
    bool isSEN() const { return (type==QT_SEN); };

    bool isU1()  const { return (type==QT_U1 && sub==0); };
    bool isParity()  const { return (type==QT_P && sub==0); };

    bool isZN()  const { return (type==QT_ZN); }; 
    bool isZ2()  const { return (type==QT_ZN && sub==2); };
    bool isZN_even() const { return (type==QT_ZN && (sub%2)==0); };
    bool isZN_odd()  const { return (type==QT_ZN && (sub%2)==1); };
    bool isA4()  const { return (type==QT_A4); }; 

    bool isAdditive() const { 
       if (type==QT_U1) {
          if (sub) wblog(FL,
            "WRN invalid type %s (%d,%d)",STR(*this),type,sub);
          return 1;
       }
       return 0;
    };

    bool permitsFerm() const { 
       return (isU1() || isZN_even() || isParity() || isSU(2));
    };

    int getZ2(int q) const { 
       if (isU1() || isZN_even() || isSU2()) { q %= 2; } else 
       if (isParity()) { 
          if (q==+1) { q=0; } else
          if (q==-1) { q=1; } else wblog(FL,
          "ERR %s() invalid symmetry label %d for %s",FCT,q,STR(*this));
       }
       else wblog(FL,"ERR %s() not yet implemented for %s",FCT,STR(*this));
       return q;
    };

     bool signedLabels() const { 
        return (type>QT_UNKNOWN && type<=QT_RANK0 && type!=QT_ZN);
     };

    int permitsOM() const { 
       if (isRank0()) {
          if (isAbelian()) { return 0; }    
          if (isA4()     ) { return 2; }    
          else wblog(FL,"ERR %s() invalid %s",FCT,STR(*this));
       }
       return (sub<=1 ? 1:2);
    };

    int permitsOM(unsigned r) const { 
       if (r<3) { return 0; } 
       if (isRank0()) {
          if (isAbelian()) { return 0; } 
          if (isA4()) { return 1; } 
          else wblog(FL,"ERR %s() invalid %s",FCT,STR(*this));
       }
       return ((r==3 && sub<=1) ? 0:1);
    };

    unsigned qlen() const {
       if (type && type<=QT_RANK0) { return 1; }
       if (type>=QT_NUM_TYPES) wblog(FL,
          "ERR invalid type `%s' (%d)",STR(*this),type);
       return sub; 
    };

    unsigned qrank() const {
       if (!type || type>=QT_NUM_TYPES) wblog(FL,
          "ERR invalid type `%s' (%d)",STR(*this),type);
       if (type<=QT_RANK0) { return 0; } 
       if (!sub) wblog(FL,   
          "ERR %s() got %s @ sub=%d",FCT,STR(*this),sub); 
       return sub; 
    };

    bool isRank0(const char *F=NULL, int L=0) const {
       if (!type || type>=QT_NUM_TYPES) wblog(F_L,
          "ERR invalid type %d, `%s'",type,STR(*this));
       return (type<=QT_RANK0? 1 : 0);
    };

    bool isLargeD(unsigned d) const {
       if (type && type<=QT_ABELIAN) { return (d>1); }
       switch (type) {
          case QT_SUN: return (d>12000 || d>(sub!=1 ? exp10(qlen()) : 20));
          case QT_SON:
          case QT_SEN:
          case QT_SpN: return (d>12000 || d>exp10(qlen()));

          case QT_UNKNOWN : return (d>100);
          default:
          wblog(FL,"ERR invalid type `%s' (%d)",STR(*this),type);
       }
       return 1; 
    };

    bool isLargeD(const wbvector<unsigned> &dd) const {
       for (unsigned i=0; i<dd.len; ++i) {
          if (isLargeD(dd[i])) return 1;
       }
       return 0;
    };

    unsigned getDimDef() const { 
       unsigned d=1; 
       if (type>QT_RANK0) {
          switch (type) {
             case QT_SUN: d=sub+1;   return d;
             case QT_SpN: d=2*sub;   return d;
             case QT_SON: d=2*sub+1; return d;
             case QT_SEN: d=2*sub;   return d;
             default:
             wblog(FL,"ERR invalid type `%s' (%d)",STR(*this),type);
          }
       }
       return d;
    };

    template<class TQ>
    size_t wdim(const TQ *qs) const {
       unsigned d=1; 
       if (type>QT_ABELIAN) {
          switch (type) { 
             case QT_SUN: d=DY::wdim_A(sub,qs); break; 
             case QT_SpN: d=DY::wdim_C(sub,qs); break; 
             case QT_SON: d=DY::wdim_B(sub,qs); break; 
             case QT_SEN: d=DY::wdim_D(sub,qs); break; 
             default:
             wblog(FL,"ERR invalid type `%s' (%d)",STR(*this),type);
          }
       }
       return d;
    };

    template<class TQ>
    size_t wdim(const qset<TQ> &qs) const { return wdim(qs.data); }

    unsigned maxDimLocal() const {
       unsigned Dloc=1; 
       if (type) {
          Dloc=maxDimLocal0();
          if (Dloc<20) { Dloc=20; } 
       }
       return Dloc;
    };

    unsigned maxDimLocal0() const {
       unsigned Dloc=2; 
       if (type>QT_ABELIAN) {

          Dloc=round(pow(double(getDimDef()),2.25+double(sub)/4.));

       }
       return Dloc;
    };

    void printDimInfo(const char *F=NULL, int L=0) const {
       wblog(F_L," *  %s: sub=%d, dim_def=%d, dloc<=%d",
       STR(*this), sub, getDimDef(), maxDimLocal());
    };

    template <class T>
    unsigned qdim(const T* q) const; 

    template <class T>
    wbvector<unsigned>& QDim(
       const T* q, unsigned r, unsigned stride,
       wbvector<unsigned> &S) const;

    bool isSelfDual() const { 
       switch (type) {
          case QT_P   : 
          case QT_SpN :
          case QT_SON : return 1; 

          case QT_SEN : return (qlen()%2 ? 0 : 1);

          case QT_U1  : 
          case QT_ZN  :
          case QT_SUN : return 0;

          default:
          wblog(FL,"ERR dual not yet defined for '%s'",STR(*this));
       }
       return 0;
    };

    template<class TQ>
    char setDual(TQ* q) const {
       switch (type) {
          case QT_P  : { return 0; } 
          case QT_U1 : { q[0]=-q[0]; return 1; }
          case QT_ZN : { q[0]=(sub-q[0])%sub; return 1; }
          case QT_A4 : {
             if (abs(q[0])==1) { q[0]=-q[0]; }
          }

          case QT_SUN: { 
             unsigned n=qlen(); if (n<2) return 0; 
             else {
                unsigned i=0, l=n-1; TQ x;  
                for (n/=2; i<n; ++i) { x=q[i]; q[i]=q[l-i]; q[l-i]=x; }
                return 1;
             }
          }
          case QT_SpN: 
          case QT_SON: return 0;

          case QT_SEN: {
             unsigned n=qlen();
             if (n<3) wblog(FL,"ERR %s() got %s",FCT,STR(*this));
             if ((n%2)==0) { return 0; } 
             else {
                unsigned l=n-2; n-=1;
                if (q[l]==q[n]) { return 0; } 
                else {
                   TQ x=q[l]; q[l]=q[n]; q[n]=x;
                   return 1;
                }
             }
          }

          default:
          wblog(FL,"ERR dual not yet defined for '%s'",STR(*this));
       }
       return 0;
    };

    template<class TQ>
    char getDual(const TQ* q0, TQ* q) const {
       if (q!=q0) { memcpy(q,q0,qlen()*sizeof(TQ)); }
       return setDual(q);
    };

    template<class TQ>
    bool isDual(const TQ *q, const TQ *x) const { 
       switch (type) {
          case QT_P  : { return x[0]== q[0]; } 
          case QT_U1 : { return x[0]==-q[0]; }
          case QT_ZN : { return x[0]==TQ((sub-q[0])%sub); }
          case QT_A4 : {
             if (abs(q[0])==1)
                  { return x[0]==-q[0]; }
             else { return x[0]== q[0]; }
          }

          case QT_SUN: { 
             unsigned i=0, l=qlen();
             if (!l) wblog(FL,"ERR %s() got n=%d",FCT,l);
             for (--l; i<=l; ++i) { if (x[i]!=q[l-i]) return 0; }
             return 1;
          }
          case QT_SpN: 
          case QT_SON: {
             if (x!=q) {
                unsigned i=0, n=qlen();
                for (; i<n; ++i) { if (x[i]!=q[i]) return 0; }
             }; return 1;
          }
          case QT_SEN: {
             unsigned n=qlen();
             if (n<3) wblog(FL,"ERR %s() got %s",FCT,STR(*this));

             if (x==q) {
                return (n%2 ? x[n-2]==x[n-1] : 1);
             }
             else {
                unsigned i=0, l=n-(n%2 ? 2 : 0);
                for (; i<l; ++i) { if (x[i]!=q[i]) { return 0; }}
                if (i<n) {
                   if (x[i]!=q[i+1] || x[i+1]!=q[i]) { return 0; }
                }; return 1;
             }
          }

          default:
          wblog(FL,"ERR dual not yet defined for '%s'",STR(*this));
       }
       return 0;
    };

    QType& operator=(const char *s) { return init(FL,s); };
    QType& operator=(const QType &q) { return init(q); };

    QType& operator=(QT_QSPACE t) {
       type=t; sub=0; if (t) { validType(FL); }
       return *this;
    };

    bool operator<(const QType &q) const {
       if (type!=q.type)
            return (type<q.type);
       else return (sub<q.sub);
    };

    char cmp(const QType &q) const { 
       if (type<q.type) return -1;
       if (type>q.type) return +1;
       if (sub <q.sub ) return -1;
       if (sub >q.sub ) return +1; else return 0;
    };

    bool operator> (const QType &q) const { return !((*this)<=q); };
    bool operator<=(const QType &q) const {
       if (type!=q.type) return (type<q.type);
       else return (sub<=q.sub);
    };

    bool operator== (const QType &q) const {
        return (type==q.type && sub==q.sub);
    };
    bool operator!= (const QType &q) const {
        return (type!=q.type || sub!=q.sub);
    };

    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char tflag) const {
       wbstring s(8); 
       char e=0; s[0]=0;

       if (type==QT_U1  ) { s=(tflag ? QT_STR[type] : "U(1)"  ); } else
       if (type==QT_P   ) { s=(tflag ? QT_STR[type] : "Parity"); } else
       if (type==QT_A4  ) { s=QT_STR[type]; } else 
       if (!type && !sub) { s=(tflag ? "???" : "unknown");}

       if (s[0]) {
          if (sub) wblog(FL,
             "ERR invalid QType %s with sub=%d",s.data,sub);
          return s;
       }

       if (sub<1 || sub>99) { e|=1; } 
       else if (type==QT_ZN) { if (sub<2) e|=2;
          if (tflag)
               { snprintf(s.data,s.len,"Z%d",  sub); }
          else { snprintf(s.data,s.len,"Z(%d)",sub); }
       }
       else if (type==QT_SUN) {
          if (tflag) 
               { snprintf(s.data,s.len,"SU%d",  sub+1); }
          else { snprintf(s.data,s.len,"SU(%d)",sub+1); }
       }
       else if (type==QT_SpN) { if (sub<2) e|=4;
          if (!tflag)
               { snprintf(s.data,s.len,"Sp(%d)",2*sub); }
          else { snprintf(s.data,s.len,"Sp%d",  2*sub);
             if (tflag=='u') { s[1]='P'; } 
          }
       }
       else if (type==QT_SON) {
          if (!tflag)
               { snprintf(s.data,s.len,"SO(%d)",2*sub+1); }
          else { snprintf(s.data,s.len,"SO%d",  2*sub+1); }
       }
       else if (type==QT_SEN) { if (sub<2) e|=8;
          if (!tflag)
               { snprintf(s.data,s.len,"SO(%d)",2*sub); }
          else { snprintf(s.data,s.len,"SO%d",  2*sub); }
       }
       else {
          if (type) wblog(FL,"ERR %s() invalid type=%d",FCT,type);
          s=(tflag ? "???" : "unknown");
       }

       if (e) wblog(FL,
          "ERR invalid QType %s with sub=%d", type < QT_NUM_TYPES ?
           QT_STR[type] : "(type out of bounds)", sub
       );

       return s;
    };

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    mxArray* toMx(const char tflag=0) const {
       const size_t n=32; wbstring s(n); 
       size_t l=snprintf(s.data,n,"%s",toStr(tflag).data);
       if (l>=n) wblog(FL,
          "WRN %s() QType too long (%s; %d)",FCT,s.data,l);
       return s.toMx();
    };

    QT_QSPACE type;
    unsigned sub; 

  private:

    int atoi(const char *s, int &k, unsigned n) {
       if (s[0]!='(') { return Wb::atoi(s,k,n); } 
       else {
          unsigned l=strlen(++s);
          if (l && s[l-1]==')') {
             wbvec<char> x(l,s); x[l-1]=0;
             return Wb::atoi(x.data,k,n);
          }
       }
       return -6;
    }
};

QType SU2("SU2");

bool operator!(const QType &q) {
   if (q.type == QT_UNKNOWN) {
      if (q.sub) wblog(FL,
         "ERR %s() invalid QType %s",FCT,STR(q));
      return 1;
   }
   return 0;
};

namespace CG {

template <class TQ> inline
unsigned gotQAlpha(const TQ* q, unsigned n, const QType *t=NULL) {

   if (!WbUtil<TQ>().isInt()) { return 0; }
   if (t && !t->qrank()) { return 0; }

   if (int(n)<0) wblog(FL,"ERR %s() got n=%d",FCT,n);

   TQ qmax=61;
   #if __APPLE__
      qmax=35; 
   #endif

   for (unsigned i=0; i<n; ++i) { if (q[i]<0 || q[i]>qmax) return 0; }

   return 1;
};

template <> inline
unsigned gotQAlpha(const char* s, unsigned n, const QType *t) {

   unsigned i=0;

   #if __APPLE__
      for (; i<n; ++i) { if (!isalnum(s[i]) || s[i]>90) break; }
   #else
      for (; i<n; ++i) { if (!isalnum(s[i])) break; } 
   #endif

   if (int(n)<0) { return i; }        
   else { return (i<n ? (i+1) : 0); } 

   return 1;
};

template <class TQ> inline 
unsigned cstr2qset(const char *s, TQ *q, unsigned n=-1) { 

   unsigned i=0;
   for (; i<n; ++i) {
      if (s[i]>='0' && s[i]<='9') { q[i]=s[i]-'0';    } else
      if (s[i]>='A' && s[i]<='Z') { q[i]=s[i]-'A'+10; } else
   #if !(defined __APPLE__)
      if (s[i]>='a' && s[i]<='z') { q[i]=s[i]-'a'+36; } else
   #endif
      break;
   }

   if (int(n)<0) { return i; } 
   else { return (i<n ? (i+1) : 0); } 
};

template <class TQ> inline 
int qset2cstr(const TQ *q, char *s, unsigned n) { 

   for (unsigned i=0; i<n; ++i) {
      if (q[i]<0) { return (i+1); }
      if (q[i]<10) { s[i]=('0'   )+q[i]; } else
      if (q[i]<36) { s[i]=('A'-10)+q[i]; } else
   #if !(defined __APPLE__)
      if (q[i]<62) { s[i]=('a'-36)+q[i]; } else
   #endif
      { return (i+1); }
   }
   return 0;
};

}; 

   map<QType, time_t> load_cstore;

QType load_store_qtype(const char *F, int L, const Wb::matFile &M);
void load_RCStore(const char *F, int L, const QVec &qvec); 

class QVec : public wbvector<QType> { 
  public:

    QVec() : wbvector<QType>() {};

    QVec(unsigned l, const QType* d0, char ref=0)
     : wbvector<QType>(l,d0,ref) {};

#ifndef NOMEX
    QVec(const char *F, int L, const mxArray *a, unsigned d=0)
     : wbvector<QType>() { init(F,L,a,d); };

     QVec& init(const char *F, int L, const mxArray *a, unsigned d=0){
        if (!a) { init(); return *this; }
        if (mxIsChar(a)) { wbstring s(a); return init(F,L,s.data,d); }
        else wblog(FL,"ERR %s() invalid QVec (%s)",FCT,mxGetClassName(a));
        return *this;
     };
#endif

     QVec& init(unsigned l=0)  {
        wbvector<QType>::init(l); return *this;
     };

     QVec& init(const QType &t)  { 
        wbvector<QType>::init(1,&t);
        return *this;
     };

     QVec& init(const char *F, int L, const char *s, unsigned d=0);

     int checkInit() const; 

     unsigned Qlen(unsigned n=-1) const;

     unsigned Qrank() const;

     unsigned Qlen(wbvector<unsigned> &dd) const; 
     unsigned Qlen(
        wbvector<unsigned> &dd, 
        wbvector<unsigned> &dz  
     ) const;

     unsigned Qpos(wbvector<unsigned> &dc) const;

     unsigned Qlenz() const {
        unsigned i=0, n=0;
        for (; i<len; ++i) { n+=(data[i].qlen() + data[i].qrank()); }
        return n;
     }

     wbvector<gTQ>& getScalar(wbvector<gTQ> &qs, char zflag=0) {
        unsigned l=0, i=0;
        for (; i<len; ++i) {
           if (data[i].isParity()) {
              if (data[i].qlen()!=1 || data[i].qrank()!=0) wblog(FL,
                 "ERR %s() %s @ %d/%d", 
                 FCT,STR(data[i]),data[i].qlen(),data[i].qrank());
              break;
           }
        }
        if (i==len) { 
           qs.init(zflag ? Qlenz() : Qlen());
        }
        else if (!zflag) {
           wbvector<unsigned> dd; qs.init(Qlen(dd));
           for (i=0; i<len; ++i) {
              if (data[i].isParity()) { qs[l]=1; }; l+=dd[i];
           }
        }
        else {
           wbvector<unsigned> dd,dz;
           unsigned l1=Qlen(dd,dz), l2=dz.sum(); qs.init(l1+l2);

           for (i=0; i<dd.len; ++i) {
              if (data[i].isParity()) { qs[l]=1; }
              l+=dd[i]; if (zflag!='z') l+=dz[i];
           }
        }
        return qs;
     };

     bool hasCG() const { 
        for (unsigned i=0; i<len; ++i) {
           if (!data[i].isAbelian()) { return 1; }}
        return 0;
     };

     unsigned qrank() const {
        unsigned r=0, i=0, q=0;
        for (; i<len; ++i) { q=data[i].qrank(); if (r<q) { r=q; }}
        return r;
     };

     template <class TQ>
     unsigned QDim(const TQ *q) const;

     template <class TQ>
     unsigned QDim(const TQ *qq, unsigned k) const;

     template <class TQ>
     wbvector<unsigned>& QDim(
        const TQ *qq, wbvector<unsigned> &S
     ) const;

     template <class TQ>
     wbvector<unsigned>& QDim(
        const TQ *qq, unsigned k, unsigned r, wbvector<unsigned> &S
     ) const;

     template <class TQ>
     wbMatrix<TQ>& getQsub(
        const wbMatrix<TQ>& Q, const wbindex &I,
        wbMatrix<TQ>& QI,
        wbMatrix<TQ> *Qx=NULL  
     ) const;

     bool allU1() const { 
        for (unsigned i=0; i<len; ++i) { if (!data[i].isU1()) return 0; }
        return 1; 
     };

     void ReduceU1() { 
        if (len) {
           for (unsigned i=0; i<len; ++i) { if (!data[i].isU1())
              return;
           }; init();
        }
     };

     char allAbelian() const {  
        char q=1; 
        for (unsigned i=0; i<len; ++i) {
           if (!data[i].isAbelian()) { q=0; break; }
           if (data[i]!=QT_U1) q|=2; 
        }
        return q;
     };

     bool anyAbelian() const {
        if (!len) return 1;
        for (unsigned i=0; i<len; ++i) {
           if (data[i].isAbelian()) return 1; }
        return 0;
     };

     bool anyParity() const {
        if (len) { 
           for (unsigned i=0; i<len; ++i) {
           if (data[i].isParity()) return 1; }
        }
        return 0;
     };

     bool allNonAbelian() const { return !anyAbelian(); };

     char isNonAbelian() const {
        if (!len) return 0;
        for (unsigned i=0; i<len; ++i) {
           if (data[i].isNonAbelian()) { return data[i].permitsOM(); }}
        return 0;
     };

     bool allAdditive() const { 
        for (unsigned i=0; i<len; ++i) {
           if (!data[i].isAdditive()) return 0;
        }
        return 1; 
     };

     int permitsOM() const {
        int q=0; if (len) { unsigned i=0; int qi;
           for (; i<len; ++i) { if (q<(qi=data[i].permitsOM())) { q=qi; }}
        }
        return q;
     };

     int permitsOM(unsigned r) const {
        for (unsigned i=0; i<len; ++i) { if (data[i].permitsOM(r)) return 1; }
        return 0; 
     };

     bool sameAs(const QVec &b, char lflag=1) const { 
        if (len==b.len || !lflag) { return (*this)==b; }
        if (len && !b.len) return   allU1();
        if (b.len && !len) return b.allU1();
        return 0;
     };

     bool sameType (const QVec &b) const {
        if (len!=b.len) return 0;
        for (unsigned i=0; i<len; i++)
        if (data[i].type!=b.data[i].type) return 0;
        return 1;
     };

     wbstring toStr() const { return toStr(0); }; 
     wbstring toStr(const char vflag) const;
     int print_qset(
        const char *F, int L, const gTQ *qs, wbvec<char> &s) const;

};

namespace DY {

template <class TQ>
class weight_info { 

  public:

     weight_info() : m(0) {};

     weight_info(const weight_info& w) { init(w); };

     weight_info(unsigned m_, unsigned n) : m(m_) {
        p.init(n); }; 

     weight_info& init(unsigned m_, unsigned n) {
        m=m_; p.init(n); return *this;
     };

     weight_info& init(const weight_info& w) {
        m=w.m; p.init(w.p);
        return *this;
     };

     weight_info& operator=(const weight_info& w) { return init(w); }

     weight_info& save2(weight_info& w) {
        m=w.m; w.m=0; w.p.save2(p);
        return *this;
     };

     int m; 

     qset<TQ> p;  

  protected:
  private:
};

template <class TQ>
class Weights { 

  public:
     Weights() : n(0) {};

     Weights(const qset<TQ> &qs) : n(qs.len), qm(qs) {
        if (!n) wblog(FL,"ERR %s() got empty qset",FCT);
        weight_info<TQ> &w=W[qs]; 
        w.init(1,n); 
        X[qs]=&w;    
     };

     Weights& init(const qset<TQ> &qs) {
        if (!(n=qs.len)) wblog(FL,"ERR %s() got empty qset",FCT);
        qm.init(qs); W.clear(); X.clear();
        weight_info<TQ> &w=W[qs]; 
        w.init(1,n); 
        X[qs]=&w;    
        return *this;
     };

     size_t dim() const {
        size_t mtot=0, merr=0;
        for (auto I=W.begin(); I!=W.end(); ++I) {
           if (I->second.m>0) { mtot+=I->second.m; }
           else { ++merr; }
        }
        if (merr) wblog(FL,
           "WRN %s() got %d undetermined weight spaces",FCT,merr);
        return mtot;
     };

     void print(
        const char *F, int L, const QType &q,
        char vflag=0) const;

     mxArray* toMx() const;

     unsigned n;         
     qset<TQ> qm;        

     map <
        qset<TQ>,        
        weight_info<TQ>  
     > W;

     map <
        qset<TQ>,  
        const weight_info<TQ>* 
     > X;

  protected:
  private:
};

class Symmetry { 

  public:

    Symmetry(const QType &q_) : q(q_) { init(q); };
    Symmetry(const mxArray *a) { q.init(FL,a); init(q); };

    Symmetry& init(const QType &q);

    template <class TQ>
    size_t getWeightsFT(
       const char *F, int L, const qset<TQ> &qs,
       Weights<TQ> &R  
    );

    template <class TQ>
    int getInnerMultFT(
       const qset<TQ> &qs, const weight_info<TQ> &w,
       Weights<TQ> &R, char vflag=0
    );

    template <class TQ>
    int Get(
       const qset<TQ> &qk, map < qset<TQ>, weight_info<TQ> > &W,
       char vflag=0
    );

    template <class TQ>
    int wExpand(const char *F, int L, Weights<TQ> &R) const;

    mxArray* toMx() const;

    QType q;          
    unsigned n;       

    wbarray<double> R; 

    wbMatrix<int> A; 

    wbvector<int> n2; 

    wbarray<int> M; 

  protected:
  private:

};

}; 

class QDir: public wbvector<char> { 

  public:

    QDir() : wbvector<char>() {};

    QDir(unsigned l, const char* d=NULL, char ref=0)
     : wbvector<char>(l,d,ref) {};

    QDir(const iTags& it) { init(it); };
    QDir(const char *F, int L,const char *s) { init(F,L,s); };
    QDir(const char *s) { init(FL,s); };

#ifndef NOMEX
    QDir(const char *F, int L, const mxArray *a) { init(F,L,a); };
#endif

    QDir& init(unsigned l=0, const char* d=NULL, char ref=0)  {
       wbvector<char>::init(l,d,ref);
       return *this;
    };

    QDir& operator=(const char *s) { return init(FL,s); };

    QDir& init(const char *F, int L,const char *s);
    QDir& init(const char *s) { return init(FL,s); };
    QDir& init(const iTags& it);

    QDir& init(const char *F, int L,const mxArray *a) {
       wbstring s(F,L,a); init(F,L,s.data);
       return *this;
    };

    QDir& init_iout(const char *F, int L,
       unsigned r, unsigned iout 
    );

    unsigned isSorted() const {
       unsigned i=0, j=0;

     #ifndef WB_SKIP_ASSERT
       for (; i<len; ++i) if (!data[i]) {
          wblog(FL,"ERR %s() invalid qdir '%s'",FCT,STR(*this));
       }; i=0;
     #endif

       for (; i<len; ++i) { if (data[i]>0) break; };
       if (i==len) return 1; 
       if (i) return 0;

       for (; i<len; ++i) { if (data[i]<0) break; }
       if (i==len) return (len+1); 

       for (j=i+1; i<len; ++i) { if (data[i]>0) break; }
       if (i==len) return j; 

       return 0;
    };

    bool operator==(const char *s) const;
    bool operator!=(const char *s) const { return !((*this)==s); };

    bool operator==(const QDir &q) const {
       return this->wbvector<char>::operator==(q); };
    bool operator!=(const QDir &q) const { return !((*this)==q); };

    bool operator<(const QDir &q) const { 
       return wbvector<char>::operator<(q); };

    QDir& Conj() {
       for (unsigned i=0; i<len; ++i) { data[i]=-data[i]; }
       return *this;
    };

    char& Conj(const unsigned &i) {
       if (i>=len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d)",FCT,i,len);
       return (data[i]=-data[i]);
    };

    unsigned nconj() {
       unsigned i=0, n=0;
       for (; i<len; ++i) { if (data[i]<0) ++n; }
       return n;
    };

    unsigned nconj(const ctrIdx &Ix) {
       unsigned n=0, i=0;
       wbvector<char> q(*this); const unsigned* ix=Ix.data;

       for (; i<Ix.len; ++i) {
          if (ix[i]>=len) wblog(FL,"ERR %s() index out of bounds"
             "(%d/%d: %d/%d)",FCT,i,Ix.len, ix[i],len);
          q[ix[i]]=0;
       }
       if (Ix.conj)
            { for (i=0; i<len; ++i) { if (q[i]>0) ++n; }}
       else { for (i=0; i<len; ++i) { if (q[i]<0) ++n; }}

       return n;
    };

    int find1() const {
       int rval=-1024; unsigned k=0, i=0, o=0, ni=0, no=0;

       for (; k<len; ++k) { if (data[k]>0) { ++ni; i=k; } else { ++no; o=k; }}
       if (no==1) { rval=-o; } else 
       if (ni==1) { rval= i; } 

       return rval;
    };

    wbstring toStr() const;
    wbstring toTag() const;

    mxArray* toMx() const { return toStr().toMx(); };

  protected:
  private:

};

   #define USR_CBUF_IDX 2

   #define usr_cbuf_ACTIVE 1

class cgdStatus { 

  public:

    cgdStatus() : ctype(0), ctime(0), mtime(0), cID(0) {
       memset(user,0,4); };

    explicit cgdStatus(CD_TYPE t) : cgdStatus() {
       init_type_(t); 
       init_time();
    };

    explicit cgdStatus(double f_) : cgdStatus() {
       init_type_(f_); 
       init_time();
    };

    cgdStatus(const char *F, int L, mxArray *a)
     : cgdStatus() { init(F,L,a); };

    template <class TQ>
    cgdStatus(const CDATA_TQ& C) { init(C.cstat); }

    cgdStatus& init() { 
       ctype=0; ctime=mtime=0; cID=0;
       memset(user,0,2); 
       return *this;
    };

    cgdStatus& init(unsigned ct, char flag=0) {
       ctype=ct;
       return init_time(flag);
    };

    cgdStatus& init(CD_TYPE t, char flag=0) {
       init_type_(t); 
       return init_time(flag);
    };

    cgdStatus& init(CD_TYPE t1, CD_TYPE t2) { 
       init_type_(t1); set(t2);
       return init_time();
    };

    cgdStatus& init_(CD_TYPE t) { 
       ctime=mtime=0; cID=0;
       memset(user,0,2); 
       init_type_(t);   
       return *this;
    };

    cgdStatus& update_m() { return init_time('m'); };
    cgdStatus& update_m(CD_TYPE t) {
       init_type_(t);   
       return init_time('m');
    };

    void swap(cgdStatus &b) {
       SWAP(ctype, b.ctype);
       SWAP(ctime, b.ctime);
       SWAP(mtime, b.mtime);
       SWAP(cID,   b.cID  );

       SWAP(user[0],b.user[0]); 
       SWAP(user[1],b.user[1]); 
    };

    cgdStatus& operator=(const cgdStatus &b) { return init(b); };

    cgdStatus& init(const cgdStatus &b) {
       if (cID && b.cID && (cID!=b.cID || ctime!=b.ctime)) wblog(FL,
          "ERR got initialized cstat with cID mismatch\n%s\n%s",
           STR2(*this,'V'),STR2(b,'V'));
       ctype=b.ctype; ctime=b.ctime; mtime=b.mtime; cID=b.cID;

       user[0]=b.user[0]; 
       user[1]=b.user[1]; 

       return *this;
    };

    CR_TYPE init(const char *F, int L, const mxArray *a);

    int sameID(const cgdStatus &b) const { 
       return (cID==b.cID && ctime==b.ctime);
    };

    inline void init_type_(CD_TYPE t) { 
       if (!t) { ctype=0; }
       else { int i=t; 
          if (i) {
             if (i<=0 || i>=CD_NUM_TYPES) wblog(FL,
               "ERR %s() invalid CData flags i=%d",FCT,i);
             ctype=(1<<(i-1));
          }
       }
    };

    inline void init_type_(double f_) { 
       unsigned f=f_; 
       if (double(f)!=f_) { wblog(FL,
          "ERR %s() invalid CData flags f=%g",FCT,f_); }
       ctype=f;
    };

    inline void init_type_( 
       const char *F, int L, double c0, double cx) {

       unsigned c2=cx;
       if (c0<0 || c0>=(1<<CD_NUM_TYPES) || cx<0 || cx>=(1<<12)) wblog(F_L,
          "ERR %s() invalid CData flags f=(%g,%g)",FCT,c0,cx);

       c2<<=(CD_NUM_TYPES);
       if (double(c2>>CD_NUM_TYPES)!=cx) wblog(F_L,
          "ERR %s() invalid CData flags f=(%g,%g)",FCT,c0,cx);

       init_type_(c0);
       ctype|=c2;
    };

    cgdStatus& set(CD_TYPE t) {
       if (t==CD_UNKNOWN || t>=CD_NUM_TYPES) wblog(FL, 
          "ERR %s() got %s",FCT,CD_TYPE_STR[t]);
       ctype|=(1<<(t-1)); 
       return *this;
    };

    cgdStatus& set(unsigned i, unsigned q=1) {
       if (i>=CD_NUM_TYPES_X)
          wblog(FL,"ERR %s() invalid CGD flag i=%d",FCT,i);
       ctype|=(q<<i);
       return *this;
    };

    cgdStatus& setComplete(unsigned q=1) { 
       if (q>7) { wblog(FL,"ERR %s() invalid q=%d",FCT,q); }
       ctype|=(q<<CD_COMPLETE);
       return *this;
    };

    cgdStatus& cpyComplete(const cgdStatus &S) { 
       return setComplete( S.isComplete(7) );
    };

    int isComplete_() const {
       int r=0; 
       if (any(CD_IDENTITY, CD_1JSY_ST3, CD_1JSY_GEN)) { r=2; } else
       if (any(CD_STD3, CD_STD3_X)) { r=3; } else
       if ((*this)==CD_ABELIAN) { r=-1; }
       return r;
    };

    unsigned isComplete(unsigned q=7) const { 
       return ((ctype>>CD_COMPLETE) & q);
    };

    char setuser_BUF(unsigned char q=1) {
       unsigned char &u=user[USR_CBUF_IDX];
       if (q) { u|=q; } else { u=0; }
       return u;
    };

    char setuser_BUF(unsigned char q, unsigned l) {
       unsigned char &u=user[USR_CBUF_IDX];
       if (l>7 || (!l && (q&1))) wblog(FL,"WRN got %s(%d,%d) u=%d",FCT,q,l,u);
       if (q) { u |=  (q<<l); }  
       else   { u &= ~(1<<l); }  
       return u;
    };

    char setuser_BUF_active(unsigned char q=1) {
       unsigned char &u=user[USR_CBUF_IDX];
       if (!(u&1)) wblog(FL,"WRN %s() got u=%d (q=%d)",FCT,u,q);
       if (q) 
            { u |= ( q<<usr_cbuf_ACTIVE   ); }
       else { u &= ((1<<usr_cbuf_ACTIVE)-1); };
       return u;   
    };

    char setuser_BUF_passive() { 
       unsigned char &u=user[USR_CBUF_IDX];
       if (!(u&1)) wblog(FL, 
          "WRN %s() got u=%s",FCT,BITS(u));
       if (u>1) { setuser_BUF_active(0); }
       if (!u ) { setuser_BUF(); }
       return u;
    };

    char gotuser_BUF() const { return user[USR_CBUF_IDX]; };
    char gotuser_BUF_active() const {
       return (user[USR_CBUF_IDX] >> usr_cbuf_ACTIVE); };

    int cmp(const char *F, int L, const cgdStatus &b) const;

    bool olderThan(const cgdStatus &b) const { return (cmp(FL,b)<0); };
    bool newerThan(const cgdStatus &b) const { return (cmp(FL,b)>0); };

    bool hasID(char lflag=0) const {
       if (!lflag)
            return (ctime && mtime && cID); 
       else return (ctime || mtime || cID);
    };

    bool operator<=(CD_TYPE t) const {
       if (t)
            { return ((ctype&((1<<CD_NUM_TYPES)-1)) <= (1U<<(t-1)) ? 1 : 0); }
       else { return (ctype ? 0 : 1); }
    };

    inline bool operator==(CD_TYPE t) const {
       if (t)
            { return (ctype & (1<<(t-1)) ? 1 : 0); }
       else { return (ctype? 0 : 1); }
    };

    inline bool operator!=(CD_TYPE t) const { return !(*this==t); };

    bool any(CD_TYPE t1, CD_TYPE t2) const {
       return ((*this)==t1 || (*this)==t2);
    };
    bool any(CD_TYPE t1, CD_TYPE t2, CD_TYPE t3) const {
       return ((*this)==t1 || (*this)==t2 || (*this)==t3);
    };

    bool operator==(const cgdStatus &b) const {
        return (ctime==b.ctime && mtime==b.mtime && cID==b.cID);
    };

    bool operator!=(const cgdStatus &b) const { return !((*this)==b); };

    bool isSet() const { return (ctime || mtime || cID); }; 
    bool isEmpty() const { return (!ctime && !mtime && !cID && !ctype); };

    int inValid() const { 
       int e=0; 
       if ((ctime!=0) ^ (mtime!=0)) { e|=1; }
       if ((*this)==CD_IMPLICIT) { 
          if (ctime || cID>99) { e|=8; }
       }
       else {
          if ((ctime!=0) ^ (cID!=0)) { e|=2; }
          if ((ctime!=0) ^ (ctype!=0)) { e|=4; }
       }
       return e;
    };

    int sameAs(const cgdStatus &b, char lflag=0) const;

    wbstring tstr() const;
    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char vflag) const;
    wbstring u2Str(const char *istr=NULL) const;

    mxArray* toMx(CR_TYPE rt=CR_DEFAULT) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    mxArray* add2MxStruct(mxArray *S, unsigned k) const;

    unsigned ctype; 

    double ctime;  
    double mtime;  
    unsigned cID;  

    unsigned char user[4]; 

  protected:
  private:

    cgdStatus& init_time(char flag=0);

    wbstring to_vstr(double t) const;

    void setID() { 
       unsigned i=0, n=3; 
       for (; i<n; ++i) {
          cID=(::rand() & ((1<<20)-1)); if (cID) break;
       }
       if (!cID) wblog(FL,"ERR %s() got cID=%d (%d/%d)",FCT,cID,i,n);
    };
};

bool operator!(const cgdStatus &q) {
   return (!q.ctype && !q.ctime && !q.mtime && !q.cID);
};

class cgrType { 

  public:

    cgrType(CR_TYPE t_=CR_DEFAULT) : t(t_) {};

    cgrType(double t_) { init(FL,t_); };

    template<class T>
    cgrType& init(const char *F, int L, T t_) {
        t=CR_TYPE(t_);
        if (T(t)!=t_ || t>=CR_NUM_TYPES) wblog(F_L,
           "ERR %s() invalid CRef type %g",FCT,double(t_)
        );
        return *this;
    };

    cgrType& operator=(const cgrType  &q) { t=q.t; return *this; };
    cgrType& operator=(const CR_TYPE &t_) { t=t_; return *this; };

    explicit operator bool() const { return (t!=0); };

    bool operator< (const CR_TYPE t_) const { return (t< t_); };
    bool operator<=(const CR_TYPE t_) const { return (t<=t_); };
    bool operator> (const CR_TYPE t_) const { return (t> t_); };
    bool operator>=(const CR_TYPE t_) const { return (t>=t_); };
    bool operator==(const CR_TYPE t_) const { return (t==t_); };
    bool operator!=(const CR_TYPE t_) const { return (t!=t_); };

    bool operator< (const cgrType &q) const { return (q.t< t); };
    bool operator<=(const cgrType &q) const { return (q.t<=t); };
    bool operator> (const cgrType &q) const { return (q.t> t); };
    bool operator>=(const cgrType &q) const { return (q.t>=t); };
    bool operator==(const cgrType &q) const { return (q.t==t); };
    bool operator!=(const cgrType &q) const { return (q.t!=t); };

    const char* tostr() const {
       if (t>=CR_NUM_TYPES) wblog(FL,"ERR %s() "
          "cgrType out of bounds (%d/%d)",FCT,t,CR_NUM_TYPES);
       return CR_TYPE_STR[t];
    };

    wbstring toStr() const { return tostr(); }

    CR_TYPE t;

  protected:
  private:
};

template <class TQ>
class qset : public wbvector<TQ> { 

  public:
    qset() : wbvector<TQ>() {};

    qset(unsigned l, const TQ* d0=NULL, char ref=0)
     : wbvector<TQ>(l,d0,ref) {};

    qset(const qset<TQ> &q)
     : wbvector<TQ>(q.len,q.data) {};

    qset(const QType &Q, const TQ *qq)
     : wbvector<TQ>() { init(Q,qq); };

    qset(const QVec &Q, const TQ *qq)
     : wbvector<TQ>() { init(Q,qq); };

    qset(const char *F, int L, const mxArray *a) : wbvector<TQ>() {
       if (mxIsChar(a)) { wbstring s(a); init(F_L,s.data); }
       else init(F,L,a,NULL,1);
    };

    qset(const char *F, int L, const mxArray *a,
       const QVec *Q, char check_type=1) : wbvector<TQ>() {
       init(F,L,a,Q,check_type);
    };

    qset(const qset<TQ> &q1, const qset<TQ> &q2)
     : wbvector<TQ>() { init(q1,q2); };

    qset(const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q3)
     : wbvector<TQ>() { init(q1,q2,q3); };

    qset(
       unsigned l1, const TQ *d1,
       unsigned l2, const TQ *d2)
     : wbvector<TQ>() { init(l1,d1,l2,d2); };

    qset(const qset<TQ> &q1, char op, const qset<TQ> &q2);

    qset& init(const qset<TQ> &q1) {
       if (this!=&q1) { wbvector<TQ>::init(q1); }
       return *this;
    };

    qset& operator=(const qset<TQ> &q1) { 
       if (this!=&q1) { wbvector<TQ>::init(q1); }
       return *this;
    };

    qset& init(const char *F, int L, const char *s, QType *t=NULL); 
    qset& init(
       const char *F, int L, const mxArray *a,
       const QVec *Q=NULL, char check_type=1
    ){
       if (mxIsChar(a)) { wbstring s(a); init(F_L,s.data); }
       else wbvector<TQ>::init(F,L,a,NULL,check_type);

       if (Q) {
          unsigned n=Q->Qlen();
          if ((this->len)%n) wblog(FL,
             "ERR %s() got invalid qset.len=%d @ %d",FCT,this->len,n
          );
       }
       return *this;
    };

    qset& init(const QType &Q, const TQ *qq) {
       wbvector<TQ>::init(Q.qlen(),qq);
       return *this;
    };

    qset& init(const QVec &Q, const TQ *qq) {
       init(Q.Qlen(),qq); 
       return *this;
    };

    qset& init(const qset<TQ> &q1, const qset<TQ> &q2) {
       init (q1.len, q1.data, q2.len, q2.data);
       return *this;
    };

    qset& init(const qset<TQ> &q1, const qset<TQ> &q2, const qset<TQ> &q3) {
       init(q1.len, q1.data, q2.len, q2.data, q3.len, q3.data);
       return *this;
    };

    qset& init() { wbvector<TQ>::init(); return *this; };

    qset& init(unsigned l, const TQ *d0=NULL) {
       wbvector<TQ>::RENEW(l,d0);
       return *this;
    };

    qset& init(
       unsigned l1, const TQ *d1,
       unsigned l2, const TQ *d2
    );

    qset& init(
       unsigned l1, const TQ* d1,
       unsigned l2, const TQ *d2,
       unsigned l3, const TQ *d3
    );

    qset<TQ>& init(const TQ q1) { return init(1U,&q1); }

    qset<TQ>& init(const TQ q1, const TQ q2) {
       return init(1,&q1,1,&q2);
    };
    qset<TQ>& init(const TQ q1, const TQ q2, const TQ q3) {
       return init(1,&q1,1,&q2,1,&q3);
    };

    qset<TQ>& initStride(
       const TQ* d_, unsigned n, unsigned r, unsigned QDIM
    ){
       init(n*r); if (d_ && this->len) {
          Wb::cpyStride(this->data,d_,n,NULL,r,-1,QDIM);
       }
       return *this;
    };

    bool operator<(const qset<TQ> &b) const;

    int checkNQs(const QVec &qq, const char *F=NULL, int L=0);

    bool isConsistent(
    const QVec &qq, int rank=-1, const char *F=NULL, int L=0) const;

    unsigned QDim(const QVec &q, int rank=-1) const;

    qset<TQ>& times(TQ fac, qset<TQ> &x) const {
        x=*this; x*=fac; return x;
    };

    wbstring toStr() const { return wbvector<TQ>::toStr(); };

    wbstring toStr(const QType &t QS_UNUSED_VAR) const {
       if (CG::gotQAlpha(this->data, this->len)) {  
          wbstring sout(this->len+1);
          if (CG::qset2cstr(this->data,sout.data,this->len)) wblog(FL,
             "ERR %s() failed to obtain compact qset string",FCT);
          sout.data[this->len]=0;
          return sout;
       }
       return wbvector<TQ>::toStr();
    };

};

template <class TQ>
bool qset<TQ>::operator<(const qset<TQ> &b) const {

   if (!this->len || !b.len) wblog(FL,
      "ERR %s() got empty object (%d/%d)",FCT,this->len,b.len);

   return wbvector<TQ>::operator<(b);
};

template <>
bool qset<double>::operator<(const qset<double> &b) const {

   if (!len || !b.len) wblog(FL,
      "ERR %s() got empty object (%d/%d)",FCT,len,b.len);
   if (len!=b.len) wblog(FL,
      "ERR %s() got length mismatch (%d/%d)",FCT,len,b.len);

   double *d=b.data;
   for (unsigned i=0; i<len; i++) {
      if (float(data[i])<float(d[i])) return 1;
      if (float(data[i])>float(d[i])) return 0; 
   }
   return 0;
};

template <>
bool qset<float>::operator<(const qset<float> &b) const {
   wblog(FL,"ERR %s() double precision prefered over '%s' (0x%lX)",
   TSTR(float),&b); 
   return 0;
};

template <class TQ>
class QSet { 

  public:

    QSet() : t(QT_UNKNOWN) {}; 

    QSet(const QSet &B) { t=B.t; qs=B.qs; qdir=B.qdir; };

    QSet(const QType &t_, const qset<TQ> &qs_,
       unsigned iout=0, 
       char isref=0
    ){
       init(t_,qs_,iout,isref);
    };

    QSet(const QType &t_, const iTags &it, unsigned r=-1){
       if (int(r)<0) { r=it.len; } else
       if (r!=it.len) wblog(FL,"ERR %s() rank mismatch (%d/%d)",FCT,r,it.len);
       t=t_; qdir.init(it); qs.init(r*t.qlen());
    };

    QSet(const CRef<TQ> &A) { init(A); };

    template <class TD>
    QSet(const char *F, int L,
       const QSpace<TQ,TD> &A, unsigned i, unsigned isym=-1
    ){ init(F,L,A,i,isym); }

    QSet(const char *F, int L, const mxArray *a) { init(F,L,a); };
    QSet(const mxArray *a) { init(0,0,a); };

    QSet(const char *F, int L,
       const CRef<TQ> &A, const ctrIdx &ica,
       const CRef<TQ> &B, const ctrIdx &icb, wbperm *P=NULL
    ){ init(F,L,A,ica,B,icb,P); };

    QSet(const char *F, int L, const char *s) { init_str(F,L,s); };

    QSet& init(const CRef<TQ> &A) {
       if (A.cgb)
            { init(*A.cgb); Permute(A.cgp); }
       else { init(); }
       return *this;
    };

    QSet& init(const char *F, int L, const mxArray *a, unsigned k=0);

    template <class TD>
    QSet& init(const char *F, int L,
       const QSpace<TQ,TD> &A, unsigned i, unsigned isym=-1);

    QSet& init(const QType t_ = QT_UNKNOWN) { t=QT_UNKNOWN;
       t=t_; if (qs.len  ) qs.init();
             if (qdir.len) qdir.init();
       return *this;
    };

    QSet& init(const QSet &B) {
       t=B.t; qs=B.qs; qdir=B.qdir;
       return *this;
    };

    QSet& init2ref(const QSet &B) { 
       t=B.t; qs.init2ref(B.qs); qdir.init2ref(B.qdir);
       return *this;
    };

    QSet& init(const QType &t_, const TQ* q, const char *qd);

    QSet& init(
       const QType &t_, const qset<gTQ>&qs_,
       unsigned iout=0, 
       char isref=0
    );

    QSet& init(
       const QType t_, const TQ* qs_, unsigned r, unsigned N, 
       unsigned iout=0  
    );

    QSet& init(
       const QType t_, const TQ* qs_, unsigned r, unsigned N, 
       const iTags &it 
    );

    QSet& init1J(const QType &t_, const TQ *q, char iflag=0) {
       unsigned m=t_.qlen();
       t=t_; qdir.init2val(2, iflag ? -1 : +1);

       qs.init(2*m); memcpy(qs.data, q, m*sizeof(TQ));
       t.getDual(q, qs.data+m);

       return *this;
    };

    QSet& init2(const QType &t_, const TQ *q=NULL) {
       unsigned m=t_.qlen();
       t=t_; qdir.init_iout(FL,2,2); qs.init(2*m);
       if (q) { unsigned l=m*sizeof(TQ);
          memcpy(qs.data,   q, l);
          memcpy(qs.data+m, q, l);
       }
       return *this;
    };

    QSet& init2(const QType &t_,
       const qset<TQ> &q1, const qset<TQ> &q2, const char *qds="++") {

       unsigned m=t_.qlen(), l=m*sizeof(TQ);
       if (q1.len<m || q2.len<m) wblog(FL,
          "ERR %s() invalid qset length (%s: %d,%d/%d)",
          FCT,STR(t_),q1.len,q2.len,m
       );
       t=t_; qdir=qds; qs.init(2*m);
       memcpy(qs.data,   q1.data, l);
       memcpy(qs.data+m, q2.data, l);
       return *this;
    };

    QSet& init1(const QType &t_, TQ *q_=NULL) {
       t=t_; qs.init(t_.qlen(),q_);
       qdir.init2val(1,+1); 
       return *this;
    };

    QSet& init3(const QType &t_) {
       t=t_; qs.init(3*t_.qlen());
       qdir.init_iout(FL,3,3); return *this;
    };

    QSet& init3(const QType &t_,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J
    ){
       unsigned n=t_.qlen();
       if (J1.len!=n || J1.len!=J2.len || J1.len!=J.len) {
          if (J1.len!=n || J2.len!=n || J.len!=n) wblog(FL,
             "ERR %s() severe length inconsistency [%d %d %d; %d]", FCT,
             J1.len,J2.len,J.len,n
          );
       }
       qs.init(J1,J2,J); t=t_;
       qdir.init_iout(FL,3,3);

       return *this;
    };

    QSet& init3(const QType &t_, 
       const TQ &J1, const TQ &J2, const TQ &J
    ){
       unsigned n=t_.qlen();
       if (n!=1) wblog(FL,
          "ERR %s() got QType %s (%d/1)",FCT,STR(t_),n);
       qs.init(3U); qs[0]=J1; qs[1]=J2;  qs[2]=J;
       qdir.init_iout(FL,3,3); t=t_;
       return *this;
    };

    QSet& init(const char *F, int L,
       const CRef<TQ> &A, const ctrIdx &ica,
       const CRef<TQ> &B, const ctrIdx &icb, wbperm *Pcgd=NULL
    );

    QSet& init(const char *F, int L,
       const QType &t_, const QDir &qd, std::initializer_list<TQ> ql);

    QSet& init_str(const char *F, int L, const char *s);

    QSet& reduceTo1J(const char *F, int L, QSet &B) const;

    QSet& Conj() { qdir.Conj(); return *this; };
    QSet& Conj(unsigned i); 

    unsigned rank(const char *F=NULL, int L=0) const {
       unsigned n=t.qlen();
       if (n*qdir.len!=qs.len) wblog(F_L, 
          "ERR %s() invalid qs (%d*%d == %d?)",FCT,n,qdir.len,qs.len);
       return qdir.len;
    };

    QSet& permute(QSet &B, const wbperm &P) const;

    QSet& Permute(const wbperm &P) {
       QSet<TQ> X(*this); X.permute(*this,P);
       return *this;
    };

    inline TQ* rec(unsigned i) const {
       unsigned r=rank(FL), n=t.qlen();
       if (i>=r || r!=qdir.len || !qdir[i]) wblog(FL,
          "ERR %s() invalid input (i=%d/%d/%d)",FCT,i,r,qdir.len);
       if (qs.len!=r*n) wblog(FL,
          "ERR %s() invalid input (len=%d = %d * %d ?)",FCT,qs.len,r,n);
       return (qs.data+i*n);
    };

    QSet& Revert(unsigned i) {
       t.setDual(rec(i)); qdir[i]=-qdir[i];
       return *this;
    };

    unsigned qdim(unsigned i) const { return t.qdim(rec(i)); };

    void swap(QSet &X) {
       if (this!=&X) {
          SWAP(t,X.t); qdir.swap(X.qdir); qs.swap(X.qs);
       }
    };

    QSet& save2(QSet &X) {
       if (this!=&X) { this->swap(X); init(); }
       return X;
    };

    QSet& operator=(const QSet &B) { return init(B); };

    bool operator==(const QSet &B) const {
       return (t==B.t && qs==B.qs && qdir==B.qdir);
    };

    bool operator!=(const QSet &B) const {
       return (t!=B.t || qs!=B.qs || qdir!=B.qdir);
    };

    char cmp(const QSet<TQ> &B) const { 
       char i=0;
       if (!(i=t.cmp(B.t))) {
       if (!(i=qdir.cmp(B.qdir))) {
            (i=qs.cmp(B.qs)); }}
       return i;
    };

    bool operator<(const QSet<TQ> &B) const { return cmp(B)<0; }

    explicit operator bool() const { return !isEmpty(); } 
    bool     operator!    () const { return  isEmpty(); }

    bool isEmpty(char check=3) const {
       if (check>3) {
          if (check=='c'              ) { check=2; } else
          if (check=='C' || check=='!') { check=3; } else
          wblog(FL,"ERR %s() invalid flag=%s",FCT,cSTR(check));
       }

       if (qs.len || qdir.len) {
          if ( ((check&1) && t==QT_UNKNOWN) ||
               ((check&2) && qs.len!=qdir.len*t.qlen()) ) {
             wblog(FL,"ERR %s() mismatch itag / QSet %s (%d*%d / %d)",
             FCT,STR(*this), qdir.len, t.qlen(), qs.len);
          }
          return 0;
       }
       return (t==QT_UNKNOWN ? 1 : 0);
    };

    bool gotRCData() const;

    int isZero() const; 

    int permissible(const char *F=NULL, int L=0) const {
       int q=0; 

       if (!qdir.len) { return q; }
       if ((q=isZero())>=0) { return (q=(q>0 ? -2 : 1)); }

       if (qdir.len<3) wblog(FL, 
          "ERR isZero() returned q=%d for %s",FCT,q,STR(*this));
       q=generateFullOM(F_L,'t'); 

       if (q>=0) { q+=1; } 

       return q;
    };

    char generateFullOM(const char *F=NULL, int L=0, char test=0) const;

    unsigned isStd3() const;

    unsigned checkStd3() const;

    int init2opt3(
       const QSet<TQ> &Q, wbperm &P, wbvector<char> &cflags);

    QSet& Convert(const char *s, wbvector<char> &cflags);

    wbstring sizeStr() const;

    bool isScalar() const; 
    bool is1J(unsigned r=2) const;

    bool isAbelian() const { return t.isAbelian(); };
    int  permitsOM() const { return t.permitsOM(rank(FL)); }; 

    int checkQ_abelian(const char *F, int L) const;

    int checkQ_SU2(const char *F, int L) const;

    QSet& Sort(wbperm *cgp=NULL, char iflag=0);

    QSet& sort(QSet &Qs, wbperm *cgp=NULL, char iflag=0) const {
       Qs.init(*this);
       return Qs.Sort(cgp,iflag);
    };

    bool isSorted() const;

    size_t memSize() const {
       return ( sizeof(t)
         + sizeof(qs  ) + qs  .len*sizeof(TQ)
         + sizeof(qdir) + qdir.len*sizeof(char)
       );
    };

    wbstring QStrS(const wbperm *cgp=NULL, char sep='|') const;
    wbstring QStr(char compact=1) const;

    wbstring toStr(const char *astr, char compact=0) const;

    wbstring toStr() const { return toStr(NULL,0); };

    wbstring toTag() const;

    mxArray* toMx() const;

    QType t;     
    qset<TQ> qs; 

    QDir qdir;   

 protected:
 private:

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;
};

template <class TQ>
unsigned QSet<TQ>::isStd3() const {

   unsigned q=0; { if (qdir!="++-") return q; } 

   unsigned m=t.qlen(), d1,d2,d3, l=t.maxDimLocal(), D=MAX(1U<<6,l), D2=D*D;

   d1=t.qdim(qs.data    );
   d2=t.qdim(qs.data+  m);
   d3=t.qdim(qs.data+2*m);

   if (qdir.len*m!=qs.len) wblog(FL,
      "ERR %s %d*%d =? %d",FCT,STR(*this),qdir.len,m,qs.len);

   q=(d1*d2<=D2 || d1<=l || d2<=l) ? 1 : 2;

   if (d1>2*d3 && (d2*d3<=D2 || d2<=l || d3<=l)) q|=4; 
   if (d2>2*d3 && (d1*d3<=D2 || d1<=l || d3<=l)) q|=8; 

   return q;
};

template <class TQ>
unsigned QSet<TQ>::checkStd3() const {

   unsigned q=0; { if (qdir.len!=3) return q; } 
   unsigned m=t.qlen(), d1,d2,d3, l=t.maxDimLocal(), D=MAX(1U<<6,l), D2=D*D;

   d1=t.qdim(qs.data    );
   d2=t.qdim(qs.data+  m);
   d3=t.qdim(qs.data+2*m);

   if (qdir.len*m!=qs.len) wblog(FL,
      "ERR %s having %d*%d =? %d",FCT,STR(*this),qdir.len,m,qs.len);

      if (d1*d2<=D2 || d1*d3<=D2 || d2*d3<=D2) q|=(1<<0);

   if (!q) { D2=1U<<16; 
      if (d1*d2<=D2 || d1*d3<=D2 || d2*d3<=D2) q|=(1<<1);
   }
   if (!q) { D2=1U<<19; 
      if (d1*d2<=D2 || d1*d3<=D2 || d2*d3<=D2) q|=(1<<2);
   }

   if (d1<=l) q|=(1<<3); 
   if (d2<=l) q|=(1<<4); 
   if (d3<=l) q|=(1<<5); 

   if (d1<d2) {
      if (d2<=d3) q|=(3<<6); 
      else        q|=(2<<6); 
   } else {
      if (d1<=d3) q|=(3<<6); 
      else        q|=(1<<6); 
   }

   return q;
};

template <class TQ>
int QSet<TQ>::init2opt3(
   const QSet<TQ> &Q, wbperm &iP, wbvector<char> &cflags) {

   unsigned i3=Q.checkStd3(), w=(i3>>6);
   unsigned m=Q.t.qlen();
   wbperm P;

   if (!(i3&63)) { return -1; }

   if (!w || w>3 || m*Q.qdir.len!=Q.qs.len) wblog(FL,
      "ERR %s() got %s\n(%d -> %d)",FCT,STR(Q),BITS(i3),BITS(w));
   if (w==1) P.init("231"); else 
   if (w==2) P.init("132"); else 
   if (w==3) P.init();           

   Q.permute(*this,P);    
   Convert("++-",cflags); 

   if (Wb::cmpRange(qs.data,qs.data+m,m)>0) {
      wbperm p2("213"); Permute(p2);
      P.Permute(p2); cflags.Permute(p2);
   }

   P.invert(iP);

   return 0;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::Convert(const char *s, wbvector<char> &cflags) {

   unsigned i=0, r=(!s || !s[0] ? 0 : strlen(s));

   if (r!=qdir.len) wblog(FL,
      "ERR %s() size mismatch %s [0/%d]",FCT,STR(*this),qdir.len);
   if (!s || !s[0]) { cflags.init(); return *this; }

   cflags.init(r);

   for (; i<r; ++i) {
       if ((s[i]=='+' && qdir[i]<0) || (s[i]=='-' && qdir[i]>0)) {
          Conj(i); cflags[i]=1;
       }
       else if (s[i]!='+' && s[i]!='-')
       wblog(FL,"ERR %s() invalid qdir='%s'",FCT,s);
   }

   return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::Conj(unsigned i) {

   unsigned m=t.qlen();
   if (i>=qdir.len || !qdir[i]) wblog(FL,"ERR %s() "
      "index out of bounds %s (%d/%d)",FCT,STR(*this),i+1,qdir.len);
   if (m*qdir.len!=qs.len) wblog(FL,"ERR %s() "
      "invalid  %s (%d*%d =? %d)",FCT,STR(*this),m,qdir.len,qs.len);

   qdir[i]=-qdir[i];
   t.setDual(qs.data+i*m);

   return *this;
};

template <class TQ>
wbstring QSet<TQ>::sizeStr() const { 
   unsigned m=t.qlen(), r=qdir.len;

   if (r*m!=qs.len) wblog(FL,
      "ERR %s: %d != %d*%d",STR(*this),qs.len,r,m);
   try {
      if (qs.len) { wbvector<size_t> dd(r);
         for (unsigned i=0; i<r; ++i) { dd[i]=t.qdim(qs.data+i*m); }
         return SSTR(dd);
      }
   }
   catch (...) {
      wblog(FL,"ERR %s() failed to obtain size(Q)\n%s",FCT,STR(*this));
   }
   return "";
};

template <class TQ>
class QHash { 

  public:

    size_t operator()(const QSet<TQ> &S, const char *istr=0) const {

       unsigned i, j=1, r=S.qdir.len, nq=S.qs.len,
          l=1+(1+nq+r)/sizeof(long), 
          lb=l*sizeof(long); 
       unsigned long  h=5381U; wbvec<unsigned long> x(l);
       char *s=(char*)x.data;

       const TQ *qs=S.qs.data;

       if ((!r && nq) || (r && nq%r)) wblog(FL,
          "ERR %s() got QSet inconsistency (%d/%d)",FCT,nq,r);
       x[l-1]=0;

       s[0]=S.t.type;
       s[1]=S.t.sub;

       for (i=0; i<nq; ++i) {
           s[++j]=char(qs[i]);
       }; ++j;

       if ((j+r)>lb) { 
          wblog(FL,"ERR %s() string/id out of bounds (%d/(%d*%d))",
          FCT,j,l,sizeof(long));
       }

       memcpy(s+j,S.qdir.data,r); 

       if (istr && *istr) {
          unsigned n=strlen(istr); if (n>lb) {
             wblog(FL,"ERR %s() istr `%s' out of bounds (%d/%d*%d)",
             FCT,istr,n,l,sizeof(long));
          }
          for (j=lb-n, i=0; i<n; ++i) { s[j+i] ^= istr[i]; }
       }

       for (i=0; i<l; ++i) { h ^= ((h<<6) + (h>>2)) + x[i]; }

       return h;
    };
};

template <class T,  ENABLE_IF_isINT(T)>
unsigned hash_unsigned(T x, unsigned m=20) {

   unsigned q=0, i=0, n=1+(8*sizeof(x)-1)/(m? m:1);

   if (m<8 || m>32) wblog(FL,"ERR %s() got m=%d",FCT,m);
   for (; i<n; ++i) { q^=x; x>>=m; }

   q &= ((1<<m)-1);
   return q;
};

template <class TD>
class cdata : public wbsparray<TD> { 

 public:

    cdata& init(SPIDX_T l=0) {
       wbsparray<TD>::init(l);
       return *this;
    };

    cdata& init(
       SPIDX_T D1, SPIDX_T D2, SPIDX_T D,
       unsigned M=1 
    ){
       if (M==1)
            wbsparray<TD>::init(D1,D2,D);
       else wbsparray<TD>::init(D1,D2,D,M);
       return *this;
    };

    cdata& init(const wbvector<SPIDX_T> &S) {
       wbsparray<TD>::init(S); return *this;
    };

    cdata& init(const cdata &C) {
       wbsparray<TD>::init(C);
       return *this;
    };

    template <class TQ>
    cdata& init(const CRef<TQ> &R, unsigned l=-1, char adapt=1); 

    template <class T2>
    cdata& init(const wbsparray<T2> &C) {
       wbsparray<TD>::init(C);
       return *this;
    };

    cdata& init( 
       const char *F, int L, const mxArray *a, unsigned k,
       const QType &q, unsigned r=-1
    );

    cdata& initScalar(double x=1, unsigned r=-1) {
       wbsparray<TD>::initScalar(x,r);
       return *this;
    };

    cdata& initIdentity(const QType &q, SPIDX_T d, TD dval=1);
    cdata& initIdentity(const wbvector<SPIDX_T> &S);

    bool isScalar() const {
       return wbsparray<TD>::isScalar(FL,'z');
    };

    unsigned len() const { return this->SIZE.len; };

    unsigned getOM(
       const char *F=NULL, int L=0, const char *istr=NULL) const;

    cdata& permute(cdata &B, wbperm P) const;    
    cdata& Permute(const wbperm &P) {            
       cdata X; this->save2(X);
       return X.permute(*this,P);
    };

    cdata& Kron(const cdata &B){
        wbsparray<TD>::Kron(B,'N','N','k');
        return *this;
    };

    TD contract( 
       const char *F, int L, const ctrIdx &ica,
       const cdata &B, const ctrIdx &icb, cdata &C,
       cdata *Cx=NULL, const wbperm *P=NULL, char normalize=1
    ) const;

    double SkipTiny(const char *F=NULL, int L=0); 

    TD NormSignC(
       const char *F=NULL, int L=0, unsigned m=0,
       TD eps=CG_EPS1,   
       TD eps2=CG_EPS2   
    );

    void info() const {
       printf("\n data: %s\n",SSTR(*this));
    };

 protected:
 private:
};

template <class TQ, class TD>
class CData : public QSet<TQ> { 

 public:

    CData() : cstat(CD_UNKNOWN) { init(); };

    CData(const CData &C)
     : QSet<TQ>((QSet<TQ>&)C), cgd(C.cgd), cstat(C.cstat) {};

    explicit CData(const QSet<TQ> &Q) : QSet<TQ>(Q) {};

    explicit CData(const CRef<TQ> &R) { init(R); }

    explicit CData(const CRef<TQ> &R, unsigned l) { init(R,l); };

    CData& init() {
       this->t.init(); this->qs.init(); this->qdir.init();
       cgd.init(); cstat.init();
       return *this;
    };

    CData& init(const QSet<TQ> &Q, bool full=true) {
       this->t=Q.t; this->qs=Q.qs; this->qdir=Q.qdir;
       if (full) {
          cgd.init(); cstat.init();
       }
       return *this;
    };

    CData& init(const CRef<TQ> &R, unsigned l=-1); 

    CData& init(const CData &B) {
       this->t=B.t; this->qs=B.qs; this->qdir=B.qdir;
       cgd=B.cgd; cstat=B.cstat;
       return *this;
    };

    template <class T2>
    CData& init(const CData<TQ,T2> &B) {
       this->t=B.t; this->qs.init(B.qs); this->qdir=B.qdir;
       cgd.wbsparray<TD>::init(B.cgd); cstat=B.cstat;
       return *this;
    };

    CData& operator=(const CData &B) { return(init(B)); }

    CData& operator=(const QSet<TQ> &Q) { return init(Q,0); };

    void swap(CData &B); 

    CData& save2(CData &B);

    CData& initAbelian(const QSet<TQ> &Q) {
       this->t=Q.t; this->qs=Q.qs; this->qdir=Q.qdir;
       cstat.init(CD_ABELIAN); cgd.init();
       return *this;
    };

    CData& initScalar(const QSet<TQ> &Q) {     
       if (!Q.t.isAbelian() && Q.qs.norm2()) { 
          wblog(FL,"ERR %s() got QSet %s",FCT,STR(Q)); }
       this->t=Q.t; this->qs=Q.qs; this->qdir=Q.qdir;
       cstat.init(CD_IMPLICIT); cgd.initScalar(1,Q.qdir.len);
       return *this;
    };

    unsigned getIQ(char flag=0) const; 

    cgdStatus& cstat_init_(CD_TYPE t, char flag=0) {
       unsigned iQ=getIQ(flag), id_=cstat.cID;
       cstat.init(t,flag); 
       if (cstat.cID!=id_) { cstat.cID=iQ; }

       return cstat; 
    };

    cgdStatus& cstat_init_(unsigned ct, char flag=0) {
       unsigned iQ=getIQ(flag), id_=cstat.cID;
       cstat.init(ct,flag);
       if (cstat.cID!=id_) { cstat.cID=iQ; }

       return cstat; 
    };

    CData& init3(const QType &t_,
       SPIDX_T D1, SPIDX_T D2, SPIDX_T D,
       unsigned M=1 
    ){
       QSet<TQ>::init3(t_);
       cgd.init(D1,D2,D, M); cstat.init(CD_UNKNOWN);
       return *this;
    };

    CData& init3(const QType &t_,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J,
       unsigned M=1 
    ){
       SPIDX_T D1,D2,D;
       QSet<TQ>::init3(t_,J1,J2,J); cstat.init(CD_UNKNOWN);

       D1=t_.qdim(J1.data);
       D2=t_.qdim(J2.data);
       D =t_.qdim(J .data); cgd.init(D1,D2,D,M);

       return *this;
    };

    CData& init3_abelian(const QType &t_, 
       const TQ &J1, const TQ &J2, const TQ &J
    ){
       QSet<TQ>::init3(t_,J1,J2,J);
       cgd.init();
       cstat.init(CD_ABELIAN); 
       return *this;
    };

    CData& init(const char *F, int L, 
       const QType &t_, const QDir &qd, std::initializer_list<TQ> ql,
       unsigned M=0); 

    int init3FT(const char *F, int L, 
       const QType &t_, const QDir &qd,
       std::initializer_list<TQ> ql,
       std::initializer_list<double> Dfull,  
       char flags=0, char mp3=0);

    int RefInit(
       const char *F, int L, const CData<TQ,TD> &B,
       char lflag=0, char bare=0);

    CData& RefInit_auxtr( 
       const char *F, int L,
       const wbvector<SPIDX_T> &S, const wbvector<RTD> &cgt
    );

    int Reduce2Ref(const char *F=0, int L=0, char force=0);
    int LoadRef   (const char *F=0, int L=0, char force=1);

    int Load_CRef(
       const char *F, int L, CData &Cr, const char *cgd_file) const;

    CData& initX3(
       const char *F, int L, const QSet<TQ> &Q,
       wbvector<char> & cflags, wbperm &P,
       unsigned loadRC=0, unsigned char flag=0); 

    int init_mxCRef( 
       const char *F, int L, const mxArray *S, unsigned k
    );

    TQ* qptr(unsigned k) const { 
       unsigned n=this->t.qlen();
       if (this->qs.len%n || k>=this->qs.len/n) wblog(FL,
          "ERR %s() index out of bounds (%d/%d <> %d)",
          FCT,this->qs.len,n,k);
       return (this->qs.data+k*n);
    };

    bool olderThan(const CData &B) const {
       if ((QSet<TQ>&)*this!=(QSet<TQ>&)B) wblog(FL,
          "ERR %s() got incompatible CData\n   %s\n<> %s",
          FCT, STR(*this), STR(B)
       );
       checkSameStat(FL,B);
       return cstat.olderThan(B.cstat);
    };

    bool operator!=(const CData &B) const { return !((*this)==B); };
    bool operator==(const CData &B) const { 
       return (
          this->t==B.t && this->qs==B.qs && this->qdir==B.qdir &&
          cgd==B.cgd && cstat.sameAs(B.cstat,0)
       );
    };

    int gotDiff(const CData &B, char lflag=0) const;
    bool sameAs(const CData &B, char lflag=0) const {
       return !gotDiff(B,lflag);
    };

    bool operator!=(const QSet<TQ> &B) const {
       return (this->QSet<TQ>::operator!=(B));
    };
    bool operator==(const QSet<TQ> &B) const {
       return (this->QSet<TQ>::operator==(B));
    };

    bool operator!=(const cgdStatus &b) const { return (cstat!=b); };
    bool operator==(const cgdStatus &b) const { return (cstat==b); };

    explicit operator bool() const { return !isEmpty(); } 

    bool isEmpty(char check=3) const {
       return (QSet<TQ>::isEmpty(check) && cgd.isEmpty() && cstat.isEmpty());
    };

    int isComplete(int q=7) const {
        q=cstat.isComplete(q);
        if (!q && !isEmpty() && isAbelian()) { q=3; }
        return q;
    };

    bool valid() const { 
       if (!cgd || cstat.isComplete(4)) { return 0; } 
       unsigned r=this->qdir.len, l=cgd.SIZE.len;
       return (l==r || (l==r+1 && r>2) ? 1 : 0);
    };

    bool NP_zero() const { 
       if (cstat==CD_BSZ_INIT) {
          if (cstat.isComplete(7)) wblog(FL,
             "WRN %s() got %s",FCT,STR2(*this,2));
          return 0;
       }
       return (!cgd && cstat.isComplete(4)); 
    };

    bool isAbelian(const char *F=NULL, int L=0) const {
       if (this->t.isAbelian()) {
          if (cstat==CD_ABELIAN) {
             if (cgd.D.len || cgd.SIZE.len) wblog(F_L,
                "ERR %s() invalid abelian CData\n%s",FCT,STR(*this));
             return 1;
          }
          if (!isScalar() || (cgd.D.len && (cgd.D.len>1 || cgd.D[0]!=1)))
          wblog(F_L,"ERR %s() invalid abelian CData\n%s",STR(*this));
          return 1;
       }
       return 0;
    };

    bool sizeable(unsigned Dmin=1024) const {
       return (cgd.D.len>Dmin);
    };

    bool isScalar(char dflag=0) const;

    TD getScalar(const char *F=NULL, int L=0) const {
       if (!isScalar()) wblog(F_L,
          "ERR %s() got non-scalar CData\n%s",FCT,STR(*this));
       return (cgd.D.len ? cgd.D.data[0] : TD(1));
    };

    unsigned rank(const char *F=NULL, int L=0) const;

    unsigned rankS() const { return cgd.rank(); };

    bool isrank(unsigned r, unsigned *r_=NULL) const {
       unsigned l=rank(FL); if (r_) { (*r_)=l; }
       return (l!=r);
    };

    template <class T>
    wbvector<T>& getSize(
       wbvector<T> &S, char bare=0, const wbperm *cgp=NULL) const;

    SPIDX_T dim() const;

    bool isRefInit(char check=1) const {
       if (cstat!=CD_REF_INIT) { return 0; }
       if (check) { 
          const unsigned r=cgd.SIZE.len, m=cgd.D.len;
          if (!r || (m>1 && m!=cgd.SIZE[r-1])) { wblog(FL, 
             "ERR %s() invalid CData size ref (%s @ D.len=%d)",
             FCT, SSTR(cgd), cgd.D.len);
          }
       }
       return 1;
    };

    int reportRefInit( 
       const char *F, int L, const char *fct, const char *istr) const {
       if (isRefInit()) {
          wblog(F_L,"WRN %s() got RefInit %s: %s",
             fct?fct:"(null)", istr?istr:"(null)", STR(*this));
          return 1;
       }
       return 0;
    };

    bool checkValidStat(const char *F=0, int L=0) const {
       if (this->t.isAbelian()) {
          if (cstat.hasID('l') || cstat!=CD_ABELIAN) {
             if (F) wblog(F,L,"ERR %s() "
                "abelian\n%s\n%s",FCT,STR(*this),STR2(cstat,'V'));
             return 0;
          }
       }
       else if (cstat.inValid() || cstat<=CD_ABELIAN) {
          if (F) wblog(F,L,"ERR %s() "
             "non-abelian\n%s\n%s",FCT,STR(*this),STR2(cstat,'V'));
          return 0;
       }
       return 1;
    };

    int checkSameStat(
       const char *F, int L, const CData &B, char lflag=1) const {

       int e=0;
       if (!lflag) { if (cstat.ctype!=B.cstat.ctype) e=1; }
       if (cstat.cID!=B.cstat.cID || cstat.ctime!=B.cstat.ctime) e|=2;

       if (e && F) {
          wblog(F_L,"WRN %s() cstat mismatch (e=%d)",FCT,e);
          wblog(FL,"  > %s\n  > %s\n... having\n  > %s\n  > %s",
             STR(*this), STR(B), STR2(cstat,'V'), STR2(B.cstat,'V'));
          if (e>1) wblog(FL,"ERR %s() e=%d",FCT,e);
       }
       return e;
    };

    bool sameType(const CData &S, const char *F=NULL, int L=0) const;

    int  cmpOM(const CData &B) const; 

    template <class DB>
    bool sameSizeR( 
       const CData<TQ,DB> &B,
       unsigned *r=NULL, const char *F=NULL, int L=0) const;

    template <class DB>
    bool sameSizeR(
       const cdata<DB> &B,
       unsigned *r=NULL, const char *F=NULL, int L=0) const;

    template <class T2>
    int sameSizeR(const char *F, int L, const wbvector<T2> &S) const;

    bool isSymmmetric() const;

    CData& initOM(const char *F, int L);

    unsigned getOM(const char *F=NULL, int L=0) const {
       unsigned m=0; 

       if (cgd.SIZE) { if (!(m=numOM(F,L))) { m=1; }} else
       if (this->qdir.len==2 && cgd.isDiag()) { m=1; } else
       if (cgd.isEmpty()) {
          if (cstat==CD_ABELIAN) { m=1; } else
          if (cstat.isComplete()) { m=0; } 
       }
       else if (!isEmpty()) { wblog(FL, 
          "WRN %s() got %s -> OM=%d",FCT,STR(*this),m); 
          this->rank(F_L); 
          cgd.wbsparray<TD>::checkSize(F_LF);
       }
       return m;
    };

    unsigned numOM(const char *F=NULL, int L=0) const {
       unsigned M=0, l=cgd.SIZE.len; 

       if (this->qdir.len && l) {
          unsigned r=this->rank(F_L); 
          if (l==r+1) { M=cgd.SIZE[r];
             if (!M || M>9999 || (M>1 && !this->t.permitsOM(r))) {
                wblog(FL,"ERR %s() invalid OM (i=%d, r=%d)\n%s",
                FCT,M,r,STR(*this));
             }
          }
          else if (l!=r) { 
             wblog(FL,"ERR %s() invalid OM setting\n%s",FCT,STR(*this));
          }
       }
       return M;
    };

    unsigned Mdims() const { 
       unsigned l=cgd.SIZE.len, r=this->qdir.len; 
       if (l!=r) { if (l<r || l>r+1 || r<=2)
           wblog(FL,"ERR %s() got r=%d with %s",FCT,r,SSTR(cgd)); }
       return l-r;
    };

    unsigned checkOM(const char *F=NULL, int L=0) const {
       if (this->qdir.len) {
          unsigned r=this->rank(F_L); 
          if (this->t.permitsOM(r)) {
             if (this->cgd.SIZE.len>r) {
                SPIDX_T i=this->cgd.SIZE[r]; if (!i)
                   wblog(FL,"ERR %s() got OM=%d",FCT,i);
                return i;
             }
             return 1;
          }
          else if (this->cgd.SIZE.len>r) wblog(FL,
          "ERR %s() CData with invalid OM\n%s",FCT,STR(*this));
       }
       return 0;
    };

    int hasFullOM(const char *F=NULL, int L=0) const;

    int completeOM_DegQ(
       const char *F, int L, const wbperm &p,
       unsigned level=0 
    );

    CData& Project(
       const char *F, int L, cdata<TD> X, 
       wbvector<TD> &w
    );

    double normDiff(
       const char *F, int L, const CData &S,
       double eps=1E14
    ) const;

    double norm2(unsigned k=-1) const; 

    wbvector<TD>& norm2(wbvector<TD> &x2, wbvector<char> *sgn=NULL) const;

    CData& NormSignC(wbvector<TD> &x);

    double NormSignC() { 
       double x=1; 
       wbvector<TD> nrm; NormSignC(nrm);

       if (nrm.len) {
          TD Dx=0, dx=0;
          for (unsigned i=1; i<nrm.len; ++i) {
             dx=Wb::abs(Wb::abs(nrm[i])-Wb::abs(nrm[i-1]));
             if (Dx<dx) { Dx=dx; }
          }
          if ((x=double(Dx))>CG_EPS1) wblog(FL,
             "ERR %s() got varying OM normalization @ %.3g",FCT,x);
          x=double(Wb::abs(nrm[0]));
       }
       return x;
    };

    wbvector<RTD>& trace(
       const char *F, int L, wbvector<RTD> &cgt) const;

    wbvector<RTD> trace(const char *F=NULL, int L=0) const {
       wbvector<RTD> cgt; return trace(F,L,cgt); 
    };

    CData& trace(const char *F, int L,
       ctrIdx &i1, ctrIdx &i2, CData<TQ,TD> &X) const;

    CData& AddMultiplicity(const char *F, int L, const cdata<TD> &c);

    CData& AddMultiplicity(const char *F, int L, const CData &C) {
       if (this->t!=C.t || this->qs!=C.qs) wblog(F_L,
          "ERR %s() incompatible CData\n[%s] <> [%s] (%s, %s)",
          FCT, STR(*this), STR(C), sizeStr().data, C.sizeStr().data
       );
       return AddMultiplicity(F_L,C.cgd); 
    };

    cdata<TD>& getMultiplicity(unsigned im, cdata<TD> &c) const;

    bool isBasicCG(const char *F=NULL, int L=0) const;

    SPIDX_T getBasicCGSize(unsigned *m=NULL) const;

    wbsparray<TD>&
    getBasicCG(unsigned k, wbsparray<TD> &a) const;

    int getCG_set(
       const char *F, int L, wbIndex &Idx, wbsparray<TD> &a) const;

    CData& reduceTo1J(const char *F, int L, CData<TQ,TD> &C) const;

    int checkConsistency(const char *F=NULL, int L=0) const;
    int checkNormSign(const char *F=NULL, int L=0, char xflag=0) const;

    int checkQ(const char *F, int L,
       const QSet<TQ> &Q, const char* istr=NULL) const;

    bool checkAdditivityZ(const char *F, int L,
       const wbMatrix<double> &Z1, 
       const wbMatrix<double> &Z2, const wbMatrix<double> &Z,
       TD eps=CG_EPS2
    ) const;

    double SkipTiny(const char *F=NULL, int L=0); 

    CData& Permute(const wbperm &P) {
       CData<TQ,TD> X; this->save2(X);
       return X.permute(*this,P);
    };

    CData& permute(CData &B, wbperm P) const { 
       QSet<TQ>::permute((QSet<TQ>&)B,P);

       cgd.permute(B.cgd,P);
       B.cstat=cstat;

       return B;
    };

    CData& Conj() { QSet<TQ>::Conj(); return *this; }

    CData& Conj(unsigned i); 

    size_t memSize() const {
       return (
          QSet<TQ>::memSize() + cgd.memSize() + sizeof(cstat)
       );
    };

    void info(const char *istr=NULL, const char *F=NULL, int L=0) const;

    wbstring dStr() const { return this->qdir.toStr(); };
    wbstring qStr() const { return this->t.toStr(); };

    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char vflag) const;

    wbstring sizeStr() const;

    mxArray* toMx() const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    mxArray* mxCreateCell(unsigned m, unsigned n) const {
       wblog(FL,"ERR %s(%d,%d)",FCT,m,n); return 0; };
    void add2MxCell(mxArray *S, unsigned i, char tst=0) const {
       wblog(FL,"ERR %s(%lX,%d,%d)",FCT,S,i,tst); };

    void put(const char *vname, const char *ws="caller"
    ) const { put(0,0,vname,ws); };

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx();
       int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F_L,"I/O putting '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };

    inline  char setuser_BUF(unsigned char q=1) {
    return cstat.setuser_BUF(q); }

    inline  char setuser_BUF(unsigned char q, unsigned l) {
    return cstat.setuser_BUF(q,l); }

    inline  char setuser_BUF_active(unsigned char q=1) {
    return cstat.setuser_BUF_active(q); }

    inline  char setuser_BUF_passive() {
    return cstat.setuser_BUF_passive(); }

    inline  char gotuser_BUF() const {
    return cstat.gotuser_BUF(); };

    inline  char gotuser_BUF_active() const {
    return cstat.gotuser_BUF_active(); }

    cdata<TD> cgd;   

    cgdStatus cstat; 

 private:
};

template <class TQ, class TD>
unsigned CData<TQ,TD>::getIQ(char flag) const {

   if (QSet<TQ>::rank(FL)<2) wblog(FL,
      "ERR %s() unexpected (empty?) QSet %s",FCT,STR(*this));

   if (flag!='f') {
      if (flag && flag!='c' && flag!='m') wblog(FL,
      "ERR %s() unexpected flag %s",FCT,cSTR(flag));
      else if (cstat.cID) wblog(FL,
      "WRN %s() already got initialized cID\n-> %s",FCT,STR(cstat));
   }

   size_t l=QHash<TQ>()( (QSet<TQ>&) *this );
   unsigned q=hash_unsigned(l,20);

   if (!q) wblog(FL,"ERR %s() got cID=%d (%lX)",FCT,q,l); 
   return q;
};

template <class TQ>
class CRef { 

 public:

    CRef( CR_TYPE rt_=CR_DEFAULT ) : cgb(NULL), rtype(rt_) {};

    CRef(const CRef &R) { init(R); };
    CRef(const CDATA_TQ &C, unsigned m=-1) : CRef() { initBase(C,m); };

    CRef& init( CR_TYPE rt_=CR_DEFAULT ) {
       cgb=NULL; rtype=rt_; cgp.init(); cgw.init();
       return *this;
    };

    CRef& init(const CRef &B, char full=1) {
       if (this!=&B) {
          cgb=B.cgb; rtype=B.rtype; cgp=B.cgp;
          if (full) { cgw=B.cgw; }
       }
       return *this;
    };

    CRef& operator=(const CRef &B) { return init(B); }

    void swap(CRef &B) { 
       if (this!=&B) {
          cgp.swap(B.cgp); SWAP(cgb,  B.cgb  );
          cgw.swap(B.cgw); SWAP(rtype,B.rtype);
       }
    };

    CRef& save2(CRef &B) {
       if (this!=&B) {
          cgp.swap(B.cgp); B.cgb=cgb;     cgb=0;
          cgw.swap(B.cgw); B.rtype=rtype; rtype=0;
          cgp.init(); cgw.init(); 
       }; return B;
    };

    double safeCpy(const char *F, int L, const CRef &B, CRef &C) const;

    CRef& init_1(const CRef &B, unsigned j) {
       unsigned m=B.wdim2(); if (j>=m) wblog(FL, 
          "ERR %s() index out of bounds (%d/%d)\nhaving %s",FCT,j,m,STR(B));
       init(B,0); B.cgw.getCol(j,cgw);
       return *this;
    };

    CRef& init_wId(const CRef &B,
        unsigned wid=1, unsigned l1=-1, const unsigned *d_=NULL
    ){
       double w=1;
       unsigned m = (wid==1 ? B.wdim1() : B.wdim2());
       if (!wid || wid>2) { wblog(FL,"ERR %s() invalid wid=%d",FCT,wid); }

       init(B,0); 

       if (int(l1)>=0) { unsigned d=qdim(l1);
          if (d_ && d!=(*d_)) wblog(FL, 
             "ERR %s() qdim inconsistency (%d/%d)",FCT,d,*d_);
          w=::sqrt(double(d));
       }
       else if (d_) {
          if (*d_<1) wblog(FL,"ERR %s() invalid qdim=%d",FCT,*d_);
          w=::sqrt(double(*d_));
       }

       cgw.initIdentity(m,0,w); 
       return *this;
    };

    CRef& init( 
       const char *F, int L, const mxArray *a, unsigned k,
       char refC=0, const QSet<TQ> *Q=NULL,
       char xflag=0 
    );

    CRef<TQ>& initIdentityR(const char *F, int L,
       const QType &q, const TQ *qs, unsigned dim=-1, char xflag=0);

    CRef<TQ>& initIdentity1J(const char *F, int L,
       const QType &q, const TQ *qs, unsigned dim=-1, char xflag=0);

    CRef<TQ>& Reduce2Identity(char xflag=0);

    wbvector<RTD>& trace(const char *F, int L, wbvector<RTD> &cgt) const;

#ifndef QS_SKIP_MPFR
    wbvector<double>& trace(const char *F, int L, wbvector<double> &cgt) const {
       wbvector<RTD> cgt_; trace(F,L,cgt_);
       return cgt.initT(cgt_);
    };
#endif

    wbvector<double> trace(const char *F, int L) const {
       wbvector<double> cgt; 
       return trace(F,L,cgt);
    };

    double trace() const;

    wbarray<double> trace(const char *F, int L,
       ctrIdx i1, ctrIdx i2, CRef *Rt=NULL) const;

    CRef& initAbelian() {
       cgb=NULL; rtype=CR_ABELIAN; 
       cgp.init(); cgw.init();
       return *this;
    };

    CRef& initAbelian(double w) {
       cgb=NULL; rtype=CR_ABELIAN; 
       cgp.init(); cgw.init(1,1); cgw[0]=w;
       return *this;
    };

    CRef& initCtrScalar(double x=1) {
       cgb=NULL; rtype=CR_CTR_SCALAR;
       cgp.init(); cgw.init(1,1,&x);
       return *this;
    };

    CRef& initBase(const CDATA_TQ *cgr, unsigned m=-1, const CDATA_TQ* X=NULL);
    CRef& initBase(const CDATA_TQ &C,   unsigned m=-1);

    int LoadRef(const char *F=0, int L=0, char force=1) const {
       if (!cgb) {
          wblog(FL,"WRN %s() got cgb=null\n%s",FCT,STR(*this));
          return 0;
       }
       return ((CDATA_TQ*)cgb)->LoadRef(F,L,force); 
    };

    int Reduce2Ref(const char *F=0, int L=0, char force=0) const {
       if (isRefInit(F,L,force)) { return 0; }
       if (!cgb) wblog(FL,"ERR %s() got cgb=null\n%s",FCT,STR(*this));
       return ((CDATA_TQ*)cgb)->Reduce2Ref(F,L,force);
    };

    explicit operator bool() const { return !isEmpty(); } 
    bool operator! () const { return isEmpty(); }

    bool isEmpty() const {
        if (cgb && !cgb->isEmpty()) { return 0; }

        size_t n=cgw.numel();           
        if (n>1 || cgp.len) wblog(FL,   
           "ERR %s() invalid abelian cref (%d)",FCT,STR(*this));
        return (rtype==CR_DEFAULT && !n); 
    };

    int isComplete() const { int q=0; 
        if (cgb) { return cgb->cstat.isComplete(); }
        else {
           if (!isAbelian()) wblog(FL,"ERR %s() got %s",FCT);
           return q=3;
        }
    };

    bool isSymmmetric(const char *F=0, int L=0) const;

    bool isRefInit(const char *F=0, int L=0, char check=1) const {
       if (cgb && cgb->isRefInit(check)) {
          if (rtype!=CR_DEFAULT) wblog(F_L,
             "ERR %s() REF_INIT mismatch\n%s",FCT,STR(*this));
          return 1;
       }
       return 0;
    };

    bool isAbelian(const char *F=NULL, int L=0) const;

    bool isScalar(char dflag=0) const;

    size_t isw2(const char *F=NULL, int L=0) const;

    size_t isw3() const {
       return (cgw.SIZE.len==3 ? cgw.numel() : 0);
    };

    unsigned Reduce_w3Id(unsigned d=-1, char lflag=0);

    int cgw_check(const char *F, int L) const;

    bool cgw_exists(unsigned i, unsigned j);  
    int  cgw_check_std3(
       const char *F=NULL, int L=-1, unsigned d3=0, unsigned im=0);

    int checkQ(const char *F, int L, const QSet<TQ> &Q) const;

    unsigned  numel() const;  
    unsigned wnumel() const { return cgw.numel(); }

    bool wscalar(const char *F=NULL, int L=0) const;
    bool wscalar1(const char *F=NULL, int L=0) const { 
       return (wscalar(F,L) && cgw[0]==1);
    };

    double wel(unsigned k) const; 

    bool wOM(const char *F=NULL, int L=0) const { 
       return (isw2(F_L)>1); };

    bool wdim_is(unsigned m) const {
       return (isw2() ? cgw.isSMatrix(m) : 0); };

    unsigned wdim(const char *F=NULL, int L=1) const;

    unsigned wdim1(const char *F=NULL, int L=0) const {
       return isw2(F_L); }   

    unsigned wdim2(const char *F=NULL, int L=0) const { isw2(F_L);
       return cgw.SIZE[1]; }

    unsigned wdim12_(const char *F, int L, unsigned &d2) const;

    unsigned wdim12(
       const char *F, int L, unsigned &d2, unsigned d0=1) const;

    unsigned wdim3(const char *F=NULL, int L=0) const;

    double wget0() const;
    double wget1() const;

    bool wSame( 
       const wbarray<double> &cgw_, char lenient=0, double eps=1e-14) const;

    template <class T>
    wbvector<T>& getSize(wbvector<T> &S, char bare=0) const;

    unsigned Size( 
       unsigned k, 
       unsigned r=-1) const;

    unsigned qdim( 
       unsigned k, 
       const QType *qt=NULL) const;

    int checkAbelian(const char *F=NULL, int L=0) const;

    bool isDiagCSC(RTD eps=1e-14) const;
    bool isIdentityCG(wbvector<double> *nrm=NULL, double eps=1e-14) const;

    int wisId() const;

    double norm2(char checks=1) const;

    double NormSignW(
       const char *F=NULL, int L=0,
       char useExt=1, 
       double eps =CG_SKIP_DEPS1, 
       double eps2=CG_SKIP_DEPS2  
    );

    double getNormSignW( 
       const char *F=NULL, int L=0,
       char useExt=1, 
       double eps =CG_SKIP_DEPS1  
    ) const;

    char sameUptoFac(const CRef &B, double *fac=NULL, double eps=1e-12) const;
    char sameAs(const CRef &B, double eps=1e-12) const;
    char sameAs_fix(const char *F, int L, CRef &B, double eps=1e-12);

    bool sameQSet(const QSet<TQ> &Q) const; 
    bool sameQSet(const CRef &B) const;
    bool sameQDir(const char *F, int L, const QDir &qd) const;

    bool sameQDir(const char *F, int L, const CRef& B) const {
       QDir qdb; B.get_qdir(qdb);
       return sameQDir(F,L,qdb);
    };

    bool sameQDir(const iTags& b) const; 

    int gotSameCData(
      const char *F, int L, const CRef& B, char xflag=0
    ) const;

    double NormStd(const char *F, int L, unsigned r=-1);

    double normExt(const char *F=NULL, int L=0) const;

    int SortDegQ(const char *F=NULL, int L=0, QSet<TQ> *Q=NULL);

    bool isSortedDegQ(wbperm *pxt, QSet<TQ> *Q=NULL) const;

    char got3(const char *F, int L,
       const QType &q, const qset<TQ> &J12, const qset<TQ> &J3
     ) const;

    bool operator==(const CRef &B) { return (sameAs(B)==0 ? 1 : 0); };
    bool operator!=(const CRef &B) { return (sameAs(B)==0 ? 0 : 1); };

    double normDiff2(const char *F, int L, const CRef &B) const;
    double normDiff (const char *F, int L, const CRef &B) const {
       return Wb::sqrt(normDiff2(F,L,B)); };

    wbarray<double>& wProd(const CRef &B, wbarray<double> &x) const;

    int check(const char *F, int L) const;
    int check(const char *F, int L,
      const QType &q, const QDir &qdir) const;

    unsigned getOM(const char *F=NULL, int L=0, char wflag=0) const;

    unsigned checkOM(const char *F=NULL, int L=0) const {
       return (cgb ? cgb->checkOM(F,L) : 0);
    };

    unsigned rankS(char lflag=0) const; 
    unsigned rank(const char *F=NULL, int L=0, char lflag=0) const;

    CRef wget(unsigned j) const { 
       CRef<TQ> B; 
       return B.init(*this,j);
    };

    CRef& Permute(wbperm P, char isnew=0); 
    CRef& permute(CRef &B, const wbperm &P, char isnew=0) const {
        B=(*this); return B.Permute(P,isnew);
    };

    CRef& Conj() {
       if (cgb) { cgp.Conj(); } else
       if (cgp) { wblog(FL,"WRN %s() got %s",FCT,STR(cgp)); }
       return *this;
    };

    int HConjOpScalar();  

    bool gotConj() const { return Wb::conj2bool(cgp.conj); };

    bool gotPerm(const char *F=NULL, int L=0) const {
       if (cgp.relevant()) {
          if (!cgb) wblog(F_L,
            "ERR %s() invalid perm (len=%d; cgb=NULL)",FCT, cgp.len);
          if (cgp.len!=cgb->qdir.len) wblog(F_L,
            "ERR %s() invalid perm (len=%d/%d)",FCT, cgp.len, cgb->qdir.len);
          return 1;
       }
       return 0;
    };

    bool anyTrafo() const { return (cgp ? 1 : 0); };
    bool  noTrafo() const { return (cgp ? 0 : 1); };

    unsigned getP(unsigned k) const {
       if (cgp.len) {
          if (k>=cgp.len) wblog(FL,
             "ERR %s() index out of bounds (%d/%d)",FCT,k,cgp.len);
          return cgp.data[k];
       }
       return k;
    };

    wbperm& getP(wbperm &p, char iflag=0) const {
       if (!cgp.len)
            return p.init(rank(FL));
       else return p.init(cgp,iflag);
    };

    char stat() const { char i=0; 
       if (rtype!=CR_DEFAULT) { i|=1; }
       if (cgp.conj) { i|=2; }
       if (cgp.len && !cgp.isIdentityPerm()) { i|=4; }
       return i;
    };

    bool affectsQDir() const;

    qset<TQ>& adapt(qset<TQ> &qs, char iflag=0) const {
       if (cgp.len) {
          wbperm P(cgp, iflag? 0:'i');
          qs.BlockPermute(P);
       }
       return qs;
    };

    QSet<TQ>& adapt(QSet<TQ> &Q, char iflag=0) const;

    ctrIdx& adapt(ctrIdx &I, char iflag=0) const;

    QDir& get_qdir(QDir &qd) const {
       if (cgb)
            { cgb->qdir.permute(qd,cgp); }
       else { qd.init(); }
       return qd;
    };

    QDir get_qdir() const { QDir qd; 
       return get_qdir(qd);
    };

    wbstring qdir2Str(char vflag=0) const {

       if (!cgb || !cgb->qdir.len) { return ""; }

       if (vflag<=0) { return get_qdir().toStr(); }
       if (vflag==1 || vflag=='v') {
          wbvec<char> s(2*(cgb->qdir.len)+5);
          s.catf(FL,"%s => %s",STR(cgb->qdir), STR(get_qdir()));
          return s.data;
       }
       else {
          wbvec<char> s(2*(cgb->qdir.len)+5);
          s.catf(FL,"%s %s", STR(cgb->qdir), STR(cgp));
          return s.data;
       }
    };

    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char lflag) const;

    wbstring sizeStr() const;

    wbstring statStr(char vflag=0) const;

    wbstring qStr() const { return cgb ? cgb->qStr() : ""; }; 
    wbstring QStr() const {
       return cgb ? cgb->QStrS(&cgp) : "";
    };

    mxArray* toMx(char flag=0) const;  

    mxArray* toMX() const { return toMx('f'); }; 

    mxArray* mxCreateStruct(unsigned m, unsigned n, char flag=0) const;
    void add2MxStruct(mxArray *S, unsigned i, char flag=0) const;

    const CDATA_TQ *cgb; 

    wbarray<double> cgw;    

    wbperm cgp; 

    cgrType rtype; 
};

   #define MAP32  std::map < qset<TQ>, CRef<TQ> >
   #define MAP31  std::map < qset<TQ>, MAP32 >

   #define MP3_LDC   1  
   #define MP3_LDR   2  

   #define TP3_ITER  3  
   #define TP3_LDM   4  
   #define TP3_LDR   8  
   #define TP3_TST  16  
   #define TP3_DBG  64  
   #define TP3_LOAD  (TP3_LDM|TP3_LDR)

   #define LB_LOAD   4  
   #define LB_UPD    8  

   #define LB_CALC  16  
   #define LB_REF   32  

   #define LB_UPD__  (LB_UPD  | LB_LOAD)          
   #define LB_CALC__ (LB_CALC | LB_LOAD)          
   #define LB_GEN    (LB_LOAD | LB_UPD | LB_CALC) 

   #define LB_ANY    (31<<2)  

   #define QF_LOADC 128

   #define CX3_SAVE   1
   #define CX3_FORCE  2

template <class TQ>
class CStore { 

  public:

    CStore() {
       BUF.max_load_factor(0.2); 
    };

    void clear() { 
       map3.clear();
       BUF.clear();
    };

    typedef typename MAP32::iterator iMAP32;
    typedef typename MAP31::iterator iMAP31;

    void init() {};

    const CDATA_TQ& getIdentityC(const char *F, int L,
       const QType &t, const TQ *qs, unsigned dim=-1,
       unsigned loadRC=LB_CALC__ 
    );
    const CDATA_TQ& getIdentity1J(
       const char *F, int L, CRef<TQ> &C,
       const QType &t, const TQ *qs, unsigned dim=-1,
       unsigned loadRC=LB_CALC__ 
    );

    void add_CData_abelian( 
       const QType &q, const TQ &J1, const TQ &J2);

    void add_CData(
       const QType &q, const qset<TQ> &J1,const qset<TQ> &J2);

    void add2BUF(const char *F, int L, const CDATA_TQ &S, char nflag=0);

    unsigned add3(const QType &q, const mxArray* S, char loadC=0);
    unsigned add2BUF(const QType &q, const mxArray* S, unsigned k=-1);

    CDATA_TQ& getBUF(
       const char *F, int L, const QSet<TQ> &Q, unsigned loadRC=0
    );

    size_t reduceMemUsage(const char *F=NULL, int L=0, double mfac=-1);

    void getQfinal_1(const char *F, int L, const QType &type,
       const qset<TQ> &J1, const qset<TQ> &J2, wbMatrix<TQ> &J,
       wbMatrix< const CRef<TQ>* > *ss=NULL, 
       char loadC=0 
    );

    void getQfinal_v(const char *F, int L, const QVec &type,
       const qset<TQ> &q1, const qset<TQ> &q2, wbMatrix<TQ> &QQ,
       wbMatrix< const CRef<TQ>* > *SS=NULL, 
       char loadC=0 
    );

    int getQfinal_zdim(
       const char *F, int L, const QVec &qvec,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q,
       wbvector<widx_t> &s1, wbvector<widx_t> &s2, wbvector<widx_t> &s,
       wbvector<unsigned> *m=NULL, wbMatrix< const CRef<TQ>* > *S3=NULL
    );

    void getQfinal(const char *F, int L, const QVec &type,
       const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, QMap<TQ> &M,
       unsigned char cgflag=0 
    );

    void Info(const char *F=0, int L=0, char vflag=0) const;

    int checkInit(const char *F, int L, const QType &t) const;

    unsigned LoadStore(const char *F, int L, const char *file);

    mxArray* toMx(const QType &q=QT_UNKNOWN) const;

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' to %s (%d)",vname,i,ws);
       if (F) wblog(F_L,"I/O putting '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };

    void put(const char *vname, const char *ws="caller") const {
       put(0,0,vname,ws);
    };

    int find_map3(const QType &q, const qset<TQ> &J1, const qset<TQ> &J2
     ) const {
       qset<TQ> J12(J1,J2);
       if (J1.len!=J2.len || J1.len!=q.qlen()) wblog(FL,
          "ERR %s() qset length inconsistency (%d/%d/%d)",
          FCT,J1.len,J2.len,q.qlen());
       return find_map3(q,J12);
    };

    int find_map3(const QType &q, const qset<TQ> &J12) const {
       if (J12.len!=2*q.qlen()) {
          if (!q.isKnown()) wblog(FL,
             "ERR %s() got uninitialized symmetry",FCT); else
          wblog(FL,"ERR %s() "
             "qset length inconsistency (%d/2*%d)",FCT,J12.len,q.qlen()
          );
       }
       auto it1 = map3.find(q);
       if (it1!=map3.end()) {
          auto it2 = it1->second.find(J12);
          if (it2!=it1->second.end()) { return it2->second.size(); }
       }
       return -1;
    };

    int valid_mp3_data(const QType &t,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J3,
       char rflag=0  
    ) const;

    CDATA_TQ* BUF_find(const QSet<TQ> &Q);

    int Reduce2Ref(
       const char *F, int L, const QSet<TQ> &Q,
       char force=0, unsigned Dmin=1024);

    map <QType,         
      map <qset<TQ>,    
        map <qset<TQ>,  
            CRef<TQ>
        >
      >
    > map3; 

    unordered_map <
       QSet<TQ>,  
       CDATA_TQ,  
       QHash<TQ>
    > BUF;

  protected:
  private:

    iMAP31 get_mp3_data(const char *F, int L,
       const QType &q, const qset<TQ> &J1, const qset<TQ> &J2,
       char force=0
    ){
       if (J1.len!=J2.len) wblog(FL, 
          "ERR qset length inconsistency (%d/%d)",J1.len,J2.len);

       iMAP31 iJ12;
       qset<TQ> J12(J1,J2);

       auto itype = map3.find(q); if (itype==map3.end()) { checkInit(FL,q);
            itype = map3.find(q); } 

       if (itype!=map3.end()) {
          iJ12 = itype->second.find(J12);
          if (iJ12!=itype->second.end()) {
             if (!iJ12->second.size()) wblog(FL,
                "ERR %s() got empty map3 data (%s)",FCT,STR(q));
             return iJ12;
          }
       }

       if (!force) {
          add_CData(q,J1,J2);
          return get_mp3_data(F,L,q,J1,J2,'f');
       }

       wblog(F,L,"ERR failed to find or generate map3 target\n"
         "%s [%s] x [%s] => [%s]",STR(q), STR(J1), STR(J2), STR(J12));
       return iJ12; 
    };

    void print_Coeff(const char *F, int L, const QType &type,
       const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J,
       double c
    ){
       wblog(F,L," *  CG coeff. %-8s [%s, %s; %s] = %.6g",
       type.toStr(), STR(J1), STR(J2), STR(J),c);
    };
};

   CStore<gTQ> gCS;

   CData<gTQ,double> EMPTY_CGSTORE;

   mxArray* map3_CreateStructMatrix(unsigned m, unsigned n);

   template <class TQ>
   void map3_add2MxStruct(mxArray *S, unsigned k,
      const QType &q, const qset<TQ> &J12);

   gTQ get_qtot_abelian(
      const QType &t, const gTQ *q0, const widx_t n,
      const widx_t *idx=NULL, const widx_t stride=1
   );

template <class TM> 
class cgc_contract_id : public wbvector<TM> { 

  public:

    cgc_contract_id() { wbvector<TM>::init(); };

    template <class TQ, class TD>
    cgc_contract_id(
       const CData<TQ,TD> &A, const ctrIdx &ica,
       const CData<TQ,TD> &B, const ctrIdx &icb
    );

    template <class TQ, class TD>
    void extract(
       CData<TQ,TD> &A, ctrIdx &ica,
       CData<TQ,TD> &B, ctrIdx &icb,
       char xflag=1 
    ) const;

    size_t memSize() const {
       return ( sizeof(*this) + this->len*sizeof(TM) );
    };

    void wblog_(const char* F=0, int L=0, const char *istr=0) const {
       #ifdef QS_USING_OMP
         Wb::ompGuard myLK(wblog_lk,1);
       #endif

       CData<gTQ,RTD> a,b; ctrIdx ica,icb;
       this->extract(a,ica,b,icb,'l');

       wblog(F_L,"TST cgc_contract_id() %s",istr ? istr : "");
       wblog(F_L,"\b a: %-32s @ %s\n b: %-32s @ %s%N\n    %s\n    %s%N",
          STR(a), STR(ica), STR(b), STR(icb),
          STR2(a.cstat,'V'),STR2(b.cstat,'V')
       );
    };

    wbstring toStr(char vflag=7) const { 
       wbstring sout; 

       CData<gTQ,RTD> a,b; ctrIdx ica,icb;
       this->extract(a,ica,b,icb, vflag & 6 ? 'l' : 0);   

       if (vflag&4) { 
          sout.init(256); snprintf(sout.data,sout.len,
            "%-20s @%s #%05X | %-20s @%s #%05X",
             STR2(a,1), STR(ica), a.cstat.cID,
             STR2(b,1), STR(icb), b.cstat.cID
          );
       }
       else if (vflag&2) {
          sout.init(128); snprintf(sout.data,sout.len,
            "%-20s @%s | %-20s @%s",
             STR2(a,1), STR(ica), STR2(b,1), STR(icb)
          );
       }
       else {
          sout.init(64); snprintf(sout.data,sout.len,
          "%s_%s * %s_%s",STR2(a,0), STR(ica), STR2(b,0), STR(icb));
       }
       return sout;
    };

  protected:
  private:

     template <class T1, class T2>
     unsigned set_val(const char *F, int L, T1* &s, const T2 &x) {
         s[0]=x;
         if (T2(s[0])!=x) { wblog(F_L,
            "ERR got non-%s data entry (%g)",TSTR(T1), double(x)); }
         s+=1; return 1;
     };

     template <class T1, class T2>
     void get_val(
        const char *F, int L, const T1* &s, T2 &x) const {
        get_range(F,L,s,1,&x);
     };

     template <class T1, class T2>
     unsigned set_vec(const char *F, int L, T1* &s, const wbvector<T2> &x) {
         s[0]=x.len;
         if (unsigned(s[0])!=x.len) { wblog(F_L,
            "ERR unexpected size (out of bounds %s @ %d)",TSTR(T1),x.len); }
         set_range(F_L,++s,x.len,x.data);
         return (1+x.len);
     };

     template <class T1, class T2>
     void get_vec(
        const char *F, int L, const T1* &s, wbvector<T2> &x) const {
        if (!s || int(s[0])<0) wblog(FL,
           "ERR %s() got invalid vec.len=%d",FCT,s?s[0]:-1);
        x.init(s[0]); ++s;
        get_range(F,L,s,x.len,x.data);
     };

     template <class T1, class T2>
     unsigned set_range(
        const char *F, int L, T1* &s, unsigned n, const T2 *x) {
        for (unsigned i=0; i<n; ++i) {
           s[i]=x[i];
           if (T2(s[i])!=x[i]) { wblog(F_L,
              "ERR got non-%s data entry (%g)",TSTR(T1), double(x[i])); }
        }
        s+=n; return n;
     };

     template <class T1, class T2>
     void get_range(
        const char *F, int L, const T1* &s, unsigned n, T2 *x) const {
        for (unsigned i=0; i<n; ++i) {
           x[i]=s[i];
           if (T1(x[i])!=s[i]) { wblog(F_L,
              "ERR got non-%s data entry (%g)",TSTR(T1), double(x[i]));
           }
        }
        s+=n;
     };

     template <class TQ, class TD>
     unsigned set_CData(const char *F, int L,
        char* &s, unsigned n, const CData<TQ,TD> &A,
        char incl_type=1 
     ){
        char *s0=s;
        unsigned r=A.rank(), l=4+A.qs.len+r+(4+r)*sizeof(unsigned);

        wbvector<size_t> S;

        if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);

        if (incl_type) {
           set_val(FL,s, A.t.type);
           set_val(FL,s, A.t.sub );
        }

        set_vec(FL,s, A.qs);
        set_vec(FL,s, A.qdir);  

        set_vec  (FL,(unsigned*&)s,A.getSize(S,'b')); 

        if ((l=s-s0)>=n) wblog(F_L,
           "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
        return l;
     };

     template <class TQ, class TD>
     unsigned get_CData(
        const char *F, int L, const char* &s, CData<TQ,TD> &A,
        char incl_type=1 
      ) const {

        A.init();
        if (incl_type) {
           A.t.initx(QT_QSPACE(s[0]),unsigned(s[1]));
           s+=2;
        }

        get_vec(F_L,s, A.qs);
        get_vec(F_L,s, A.qdir);  

        get_vec(F_L,(const unsigned*&)s, A.cgd.SIZE);
        A.cstat.init_(CD_BSZ_INIT);

        return 0;
     };

};

template <class TM>
class MHash { 

  public:

    size_t operator()(const cgc_contract_id<TM> &S) const {

       MTI h=5381U, *const x=S.data;
       for (unsigned i=0; i<S.len; ++i) {
          h ^= ((h<<6) + (h>>2)) + x[i];
          h ^= ((h<<6) + (h>>2)) + x[i];
       }
       return h;
    };
};

template <class TQ, class TD>
class x3map { 

  public:

    x3map() : cgb(NULL), zflag(0), rtype(CR_DEFAULT)
    #ifndef QS_USING_MPFR
    , X3(x3)
    #endif
    {};

    x3map(const char *F, int L,
       const CRef<TQ> &A_, const ctrIdx &ica,
       const CRef<TQ> &B_, const ctrIdx &icb, wbperm &cgp,
       char xCGR=1) 
     : x3map() { contract_x3(F_L,A_,ica,B_,icb, cgp, xCGR); }; 

    x3map& init() {
       a.init(); b.init(); c.init(); pab.init(); x3.init();
       #ifdef QS_USING_MPFR
          X3.init();
       #endif
       cgb=NULL; zflag=0; rtype=CR_DEFAULT;
       return *this;
    };

    int contract_x3(const char *F, int L, 
       const CRef<TQ> &A_, const ctrIdx &ica_,
       const CRef<TQ> &B_, const ctrIdx &icb_, wbperm &cgp,
       char xCGR=0); 

    x3map& init(
       const char *F, int L, const mxArray* S, unsigned k=0,
       ctrIdx *ica=NULL, ctrIdx *icb=NULL, char recalc=1); 

    bool isEmpty() const {
       return (
          a.isEmpty() && b.isEmpty() && c.isEmpty() &&
          !pab.len && x3.isEmpty()
       );
    };

    wbstring toStr() const { return toStr(0); }; 
    wbstring toStr(char vflag) const;

    wbstring x3Str(char vflag=1) const;

    void print(
        const char *F=NULL, int L=0, const char *istr="", char vflag=0,
        const ctrIdx *ica=NULL, const ctrIdx *icb=NULL
      ) const;

    size_t memSize() const {
       return ( sizeof(cgb) + sizeof(char)
         + a.memSize() + b.memSize() + c.memSize()
         + sizeof(pab) +  pab.len*sizeof(widx_t)
         + sizeof(x3) + x3.numel()*sizeof(TD)
         + sizeof(rtype)
       );
    };

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;

    void add2MxStruct(
       mxArray *S, unsigned l,
       const cgc_contract_id<MTI> *idc=NULL) const;

    mxArray* toMx(const cgc_contract_id<MTI> *idc=NULL) const { 
       mxArray *S=mxCreateStruct(1,1);
       add2MxStruct(S,0,idc);
       return S;
    };

    CData<TQ,RTD> a, b, c;

    const CDATA_TQ* cgb; 

    wbperm pab; 

    int zflag; 

    cgrType rtype;

    wbarray<TD> x3;  

#ifdef QS_USING_MPFR
    wbarray<RTD> X3; 
#else
    wbarray<TD> &X3; 
#endif

  protected:
  private:
};

template <class TQ, class TD>
class X3Map { 

  public:

    X3Map() {
       XBUF.max_load_factor(0.2); 
    };

    void clear() {
       XBUF.clear(); 
    };

    int init(const char*F, int L, const mxArray *S, unsigned k=0);

    x3map<TQ,TD>& getXBUF(  
       const char *F, int L, 
       const cgc_contract_id<MTI> &idc
    );

    int contains(const char *F, int L,
       const x3map<TQ,TD> &x, const cgc_contract_id<MTI> &idc);

    int contractCGR(const char *F, int L, 
       const CRef<TQ> &A, const ctrIdx &ica,
       const CRef<TQ> &B, const ctrIdx &icb, CRef<TQ> &C,
       char xCGR=0 
    );

    CRef<TQ>& contractDegQ(const char *F, int L,
       const CRef<TQ> &A, CRef<TQ> &B);

    unsigned add(
       const char *F, int L, const cgc_contract_id<MTI> &idc,
       const mxArray* S, unsigned k
    );

    mxArray* toMx() const;

    void Info(const char *F=0, int L=0) const;

    unordered_map <
       cgc_contract_id<MTI>,
       x3map<TQ,TD>,  
       MHash<MTI>
    > XBUF;

  protected:
  private:
};

   X3Map<gTQ,double> gXS; 

template <class TQ>
class QMap { 

  public:

    void init() { qvec.init();
       Q1.init(); Q2.init(); Q.init();
       I1.init(); I2.init(); II.init(); D.init();
    };

    QMap<TQ>& init(const char *F, int L,
       const QVec &qv, const wbMatrix<TQ> &q1, const wbMatrix<TQ> &q2,
       const wbMatrix< wbMatrix<TQ> > &QQ,
       const wbMatrix< wbMatrix< const CRef<TQ>* > > *SS=NULL
    );

    QMap<TQ>& init(const char *F, int L, const mxArray *a);

    template<class TD>
    int getCGZlist(const char *F, int L,
       wbMatrix<TQ> &qc,
       wbvector<TD> &dc,
       wbvector<widx_t> *IC=NULL,
       const wbvector<widx_t>* Ix=NULL, 
       const qset<TQ> *q3=NULL,     
       wbIndex *m3=NULL,  
       const wbvector<TD>* cc=NULL,
       double eps1=1e-10,
       double eps2=1e-14  
    ) const;

    void Skip(const wbindex &I);

    template<class TD>
    void getIdentityQ(const char *F, int L,
       const wbvector<widx_t> &Sa, const wbvector<widx_t> &Sb,
       QSpace<TQ,TD> &A, char vflag=1
     ) const;

    bool isConsistent(
       const char *F=NULL, int L=0, const QVec *qv=NULL
     ) const;

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i, char tst=0) const;

    mxArray* toMx() const {
       mxArray *S=mxCreateStruct(1,1);
       add2MxStruct(S,0); return S;
    };

    void put(const char *vname, const char *ws="caller") const
    {  put(0,0,vname,ws); }

    void put(const char *F, int L,
       const char *vname, const char *ws="caller"
    ) const {
       mxArray *a=toMx(); int i=mexPutVariable(ws,vname,a);

       if (i) wblog(F_L,
          "ERR failed to write variable `%s' (%d)",vname,i);
       if (F) wblog(F_L,"I/O putting '%s' to %s",vname,ws);
       mxDestroyArray(a);
    };

    QVec qvec; 

    wbMatrix<TQ> Q1,Q2,Q;

    wbvector<widx_t> D; 

    wbvector<widx_t> I1, I2;
    wbMatrix<widx_t> II;

    wbMatrix< const CRef<TQ>* > cg3;

};

template <class TQ>
bool QMap<TQ>::isConsistent(const char *F, int L, const QVec *qv) const {

   unsigned d=qvec.Qlen();

   if (Q1.dim2!=d || Q2.dim2!=d) {
      if (!F) return 0; else wblog(F,L,
      "ERR QMap: invalid qset record length (%d,%d/%d)",
       Q1.dim2, Q1.dim2, d);
   }
   if (qv && (!qvec.sameType(*qv) || d!=(qv->Qlen()))) {
      if (!F) return 0; else wblog(F,L,
      "ERR QMap: Q type mismatch: %s; %s (%d/%d)",
       STR_(qv), STR(qvec), qv->Qlen(), d);
   }
   return 1;
};

template <class T>
class pairPatternCG { 
  public:

    pairPatternCG& init(unsigned i0=0, unsigned j0=0, unsigned n=0) {
       i=i0; j=j0; k.init(n); fac.init(n);
       return *this;
    };

    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    unsigned i,j;
    wbvector<unsigned> k;
    wbvector<T> fac;

  protected:
  private:
};

template <class TQ, class TD>
class genRG_base{ 

  public:

     genRG_base() : err(0) {};
     genRG_base(const genRG_base& S) { init(S); }

     bool isEmpty() const {
        return (!q && !Sp.len && !Sz.len && !J.len && Z.isEmpty());
     };

     double normDiff(const char *F, int L,
        const genRG_base& B,
        TD eps=1E12 
      ) const;

     void checkOrthoCommRel(
        const char *F, int L, const genRG_base &B) const;

     genRG_base& Sort();
     double SkipTiny(const char *F=0, int L=0); 

     genRG_base& ApplyQFac(
        const wbvector<double> &qfac, const wbarray<double> *JM=NULL);

     genRG_base& init() {
        q=0; J.init(); Z.init(); Sp.init(); Sz.init();
        err=0; istr.init(); 
        return *this;
     };

     genRG_base& init(const genRG_base& S) {
        q=S.q; J=S.J; Z=S.Z; Sp=S.Sp; Sz=S.Sz;
        err=S.err; istr=S.istr; 
        return *this;
     };

     genRG_base& operator=(const genRG_base& S) { return init(S); }

     SPIDX_T dim() const {
        if (!Sz.len || Sz.len!=Sp.len) wblog(FL,
           "ERR %s() got empty generators for %s [%s]",
           FCT,STR(q),STR(J),Sz.len,Sp.len);
        SPIDX_T d=Sz[0].dim();
        if (d!=Sp[0].dim()) wblog(FL,"ERR %s() got "
           "inconsistent Sp/z dimensions (%d/%d)",FCT,Sp[0].dim(),d);
        return d;
     };

     void compareStdSU2(const char *F=NULL, int L=0) const;

     unsigned block_decompose_Sm(const char *F, int L, char w,
        wbMatrix< wbarray<TD> > &SM,
        wbMatrix<unsigned> &IJ, WBINDEX *D=NULL) const;

    int get1J_gen_aux(const char *F, int L,
       const genRG_base<TQ,TD> &G2, const CData<TQ,TD> &Z,
       wbsparray<TD> &C, double &dx, double &dn, char flag=0) const;

     genRG_base& save2(genRG_base& S) {
        S.q=q; J.save2(S.J); Z.save2(S.Z); q=0;
        Sp.save2(S.Sp); Sz.save2(S.Sz);
#ifdef CG_CHECK_MW_PERM
        P0.save2(S.P0);
#endif
        S.err=err; istr.save2(S.istr); 
        return S;
     };

     wbstring info() const { 
        wbvec<char> s(128); char sep[2]=" ";
        if (!q.isAbelian() && J.wbvector<TQ>::allIn(0,9)) { sep[0]=0; }
        s.catf(FL,"%s [%s]",STR(q), J.toStrf("",sep,q.qlen(),";").data);
        return s.data;
     };

     mxArray* toMx() const {
        mxArray *a=mxCreateStruct(1,1);
        add2MxStruct(a,0);
        return a;
     };

     size_t memSize() const {
        size_t mtot=sizeof(QType) + J.len*sizeof(TQ);
        unsigned i=0;
        for (   ; i<Sp.len; ++i) mtot+=Sp[i].memSize();
        for (i=0; i<Sz.len; ++i) mtot+=Sz[i].memSize();
        return mtot;
     };

     mxArray* mxCreateStruct(unsigned m, unsigned n) const;
     void add2MxStruct(mxArray *S, unsigned i) const;

     QType q;

     wbvector< SPARRAY_TD > Sp; 
     wbvector< SPARRAY_TD > Sz; 

     qset<TQ> J;

     wbMatrix<double> Z; 

#ifdef CG_CHECK_MW_PERM
     wbperm P0;
#endif

     double err;    
     wbstring istr; 

  protected:
  private:
};

template <class TQ, class TD>
class genRG_struct { 
  public:

    genRG_struct() {};

    genRG_struct& checkInit(const char *F, int L, const QType &q0) {
       if (q0.isNonAbelian()) {
          unsigned n=RSet.size();
          if (n>0) {
             if (!q0.type) wblog(FL,"ERR %s() got %d elements "
                "in buf for '%s'",FCT,n,STR(q));
             return *this;
          }
          SetupSym(F,L,q0);
       }
       return *this;
    };

    genRG_struct& SetupSym(
       const char *F, int L, const QType &q0, qset<TQ> *qs_=NULL);

    genRG_struct& Setup_A4 (const char *F, int L, qset<TQ> *qs=NULL);

    genRG_struct& Setup_SUN(const char *F, int L, qset<TQ> *qs=NULL);
    genRG_struct& Setup_SpN(const char *F, int L, qset<TQ> *qs=NULL);
    genRG_struct& Setup_SON(const char *F, int L, qset<TQ> *qs=NULL);
    genRG_struct& Setup_SEN(const char *F, int L, qset<TQ> *qs=NULL);

    void initSU2(const qset<TQ> &J);

    void initCommRel(const char *F, int L,
       const genRG_base<TQ,TD> &R, char vflag=0);

    double checkCommRel(
      const char *F, int L, const genRG_base<TQ,TD> &R
    ) const;

    void printStatus(const char *F, int L, const char *istr=NULL);

    int get1J_gen_old(const char *F, int L, const qset<TQ> &q, char flag=0);

    int get1J_gen(const char *F, int L, const qset<TQ> &q, char flag=0);

    unsigned block_decompose_Sm(const char *F, int L,
      const qset<TQ> &J, char w, wbMatrix< wbarray<TD> > &SM);

    int getTensorProdReps_gen(
       const qset<TQ> &q1, const qset<TQ> &q2,
       unsigned flag=TP3_LDM); 

    void getTensorProdReps(
       const qset<TQ> &q1, const qset<TQ> &q2
    ){
       q.validType(FL);
       int n=getTensorProdReps_gen(q1,q2,TP3_TST);
       if (n>=0) return;

       wblog(FL,"ERR %s() "
         "failed to access RSets of input spaces (%d)",FCT,n);

    };

    unsigned genTensorProds(const qset<TQ> &q, unsigned flag=0);
    unsigned genTensorProds(unsigned dmax=-1, char sdig=0, char vflag=0);

    void put(const char *vname, const char *ws="caller") const
    {  put(0,0,vname,ws); }

    mxArray* toMx(char bflag=1) const;
    void put(const char *F, int L,
      const char *vname, const char *ws="caller"
    ) const;

    void wblog_RC_1st_build(
       const char *F, int L, const char *fct, int nrep) const {

       if (CG_VERBOSE && !q.isSU2()) {
          if (q.sub>1) { PRINTF("\n"
          "   RC_STORE %s() first-time build for %s\n"
          "   %d passes, please wait ...\n\n",fct,STR(q),nrep); }
          else wblog(F_L,
             "NB! %s() RC_STORE first-time build for %s",fct,STR(q));
       }
    };

    QType q;

    map <qset<TQ>,          
         genRG_base<TQ,TD>  
    > RSet;

    qset<TQ> qdef;

    wbvector< pairPatternCG<double> > CR;

    wbMatrix<double> DZ;

  protected:
  private:
};

template <class TQ>
TQ getTensorProdReps_abelian(const QType &q, const TQ &J1, const TQ &J2) {

   TQ J=0; 

   switch (q.type) {
      case QT_U1: J=J1+J2; break;

      case QT_ZN: J=int(J1+J2)%q.sub;

         if (TQ(int(J1)%q.sub)!=J1 || TQ(int(J2)%q.sub)!=J2)
            wblog(FL,"ERR %s() invalid %s labels (%g,%g)",
            FCT, STR(q), double(J1), double(J2));
         break;

      case QT_P: J=J1*J2;

         if (fabs(double(J1))!=1 || fabs(double(J2))!=1)
            wblog(FL,"ERR %s() invalid %s labels (%g,%g)",
            FCT,STR(q),double(J1),double(J2));
         break;

      default: wblog(FL,"ERR %s() got invalid symmetry %s",FCT,STR(q));
   }

   return J;
};

template <class TQ, class TD>
class RStore { 

 public:

   RStore(const char *s=NULL) {
      if (s) {
         if (!strcmp(s,"info")) {
            #ifdef __WB_MPFR_HH__
               wblog(FL,"<i> %s() uses type %s (=mpfr@%d) ",
               FCT, TSTR(TD), TD().prec());
            #else
               wblog(FL,"<i> %s() uses type %s",FCT,TSTR(TD));
            #endif
         }
         else wblog(FL,"ERR %s() invalid usage (%s)",FCT,TSTR(TD));
      }
   };

   size_t bufsize(const QType &q) {
      auto it=buf.find(q);
      return (it==buf.end() ? 0 : it->second.RSet.size());
   };

   bool isEmpty(const QType &q) {
      return (bufsize(q)>0 ? 1 : 0);
   };

   genRG_struct<TQ,TD>& Buf(const QType &q) {
      auto it=buf.find(q);
      if (it!=buf.end()) {
         genRG_struct<TQ,TD> &B=(it->second);
         if (!B.q.isKnown()) {
            if (B.RSet.size()) wblog(FL,
               "ERR %s() got RSet with unknown qtype (len=%d; %s)",
               FCT,B.RSet.size(),STR(B.q));
            B.checkInit(FL,q);
         }
         return B;
      }
      else {
         genRG_struct<TQ,TD> &B=buf[q]; 
         B.checkInit(FL,q);
         return B;
      }
   };

   genRG_base<TQ,TD>* find_RSet(const QType &q, const TQ *qptr){
      genRG_struct<TQ,TD> &B=Buf(q);
      qset<TQ> qs(q.qlen(),qptr,'r');

      auto it = B.RSet.find(qs);
      return (it!=B.RSet.end() ? &(it->second) : NULL);
   };

   genRG_struct<TQ,TD>& getRSet(
      const QType &q,
      qset<TQ> *qvec=NULL 
   ){
      genRG_struct<TQ,TD> &B=Buf(q);
      if (qvec) { (*qvec)=B.qdef; }
      return B;
   };

   genRG_base<TQ,TD>& getR(
      const char *F, int L, const QType &qtype, const TQ *q) const;

   genRG_base<TQ,TD>& getR(
      const char *F, int L, const CDATA_TQ &S,
      unsigned k 
   ) const { return getR(F_L, S.t, S.qptr(k)); };

   wbMatrix<TQ> getZ( 
      const char *F, int L, const QType &qtype, const TQ *q
   ) const {
      const genRG_base<TQ,TD> &R = getR(F_L,qtype,q);
      return R.Z;
   };

   unsigned qdim(
      const char *F, int L, const QType &qtype, const TQ *q
    ) const {

      if (qtype.isAbelian()) return 1;
      const genRG_base<TQ,TD> &R = getR(F_L,qtype,q);

      if (qtype.isSU2() && TQ(R.Z.dim1)!=q[0]+1) wblog(FL,"ERR %s() "
         "got multiplet dim mismatch (%d/%d)",FCT,R.Z.dim1,int(q[0]+1));

      return R.Z.dim1;
   };

   mxArray* toMx(const QType &q=QT_UNKNOWN, char bflag=1) const;

   unsigned LoadStore(const char *F, int L, const char *file);
   unsigned add(const QType &q, const mxArray* S, unsigned k);

   void Info(const char *F=0, int L=0) const; 

   map <QType,              
        genRG_struct<TQ,TD> 
   > buf;

   protected:
   private:
};

   RStore<gTQ,RTD> gRS; 

#ifdef QS_USING_OMP
   Wb::ompNLock RCS_lock("RCStore"); 
#endif

   map<QType,int> g_sym; 

class RCStore { 

 public:

   RCStore() { };

   int setupDirs(const char *F=NULL, int L=0) {

      if (root.idx.len) { return -root.idx.len; } 

     #ifdef QS_USING_OMP
      Wb::ompGuard dLK(RCS_lock); 
      if (root.idx.len) { return -root.idx.len; }  
     #endif

      int i1,i2=0,e=0;

      const char rcs[]="RC_STORE", rsy[]="RC_SYNC";

      i1=setup_root_directory(F_L,rcs,root);
      if (!i1) wblog(FL,"ERR missing env variable '%s'",rcs);

      i2=setup_root_directory(F_L,rsy,root,'+');

      if (i1<0 || i2<0) { root.init(); wblog(FL, 
         "ERR invalid environemental variables RC_STORE%s",
          i2<0 ? " or RC_SYNC":""); }
      if (i2>1) { root.init(); wblog(FL, 
         "ERR %s() only one RC_SYNC directory allowed",FCT); }

      for (auto it=root.begin(); it.end(); ++it) {
         if (!strncmp(it.at(),"/home/",7)) { e=1;
            wblog(FL,"WRN %s() RC_%s points to HOME (%d/%d)\n"
              "(recommended to use /data or similar, instead)",
              FCT, !it.i || it.i+1<it.n?"STORE":"SYNC", it.i+1,it.n);
            root.print("RC_STORE+SYNC"); break;
         }
      }
      if (!e && CG_VERBOSE>6) { root.print("RC_STORE+SYNC"); }

      if (!i2 && i1) {
         root.append(FL,root.last(),":"); 
      }
      return root.idx.len;
   };

   template <class TQ, class TD>
   int save_1J(const char *F, int L,
      const CData<TQ,TD> &A, const char *istr=0);

   template <class TQ, class TD>
   int save_CData(const char *F, int L, const CData<TQ,TD> &A);

   template <class TQ, class TD>
   int load_CData(
      const char *F, int L, const QSet<TQ> &Q, CData<TQ,TD> &A,
      char ref=0); 

   template <class TQ>
   int save_mp3( 
      const char *F, int L, const QType &t, const qset<TQ> &J12,
      const wbvector<double> *c2eps=NULL) const;

   template <class TQ>
   int save_mp3_extended(
      const char *F, int L, const QType &t, const qset<TQ> &J12,
      const wbvector<double> *c2eps=NULL) const;

   template <class TQ>
   int load_mp3(const char *F, int L, 
       const QType &t, const qset<TQ> &J1, const qset<TQ> &J2,
       unsigned loadRC=0) const;

   template <class TQ, class TD>
   int save_RSet(
      const char *F, int L, const genRG_base<TQ,TD> &R, char xflag=0);

   template <class TQ>
   int load_RSet(
      const char *F, int L, const QType &t, const qset<TQ> &J,
      char vflag=0);

   int save_XMap( 
      const char *F, int L, const cgc_contract_id<MTI> &idc) const;

   int load_XMap(
      const char *F, int L, const cgc_contract_id<MTI> &idc) const;

   template <class TQ>
   int get_file_name(
      const char *F, int L, wbstring &file,
      const QSet<TQ> &Q, const char *ext=NULL, char mflag=0) const;

   template <class TQ, class TD>
   int get_file_name(
      const char *F, int L, wbstring &file,
      const genRG_base<TQ,TD> &R, const char *ext=NULL, char mflag=0
    ) const;

   int get_file_name(
      const char *F, int L, wbstring &file,
      const cgc_contract_id<MTI> &idc, const char *ext=NULL, char mflag=0,
      QType *t=NULL
    ) const;

   size_t rclog(const QType &t,
      const char *F, int L, const char *fmt, ...);

   size_t rclog(const QType &t,
      const char *F, int L, char vflag, const char *fmt, ...);

   size_t vrclog(const QType &t, 
      const char *F, int L, const char *fmt, va_list args);

   void Info(const char *F=0, int L=0) { 
      if (F) wblog(F,L," *  %s()\n",FCT);
      root.print("RCStore::root");
      fprintf(stdout,"\n\n");
   }

   int is_new(const QType &t) {
       return (g_sym[t] ? 0 : 1);
   };

   Wb::Path root; 

 protected:
 private:

   int get_rcs_path( 
       const char *F, int L, char *path, unsigned len, const QType &t,
       const char *sub1=NULL, const char *sub2=NULL, const char *file=NULL,
       char mflag=0
     ) const;

   int setup_root_directory(
       const char *F, int L, const char *env, Wb::Path &root,
       char append=0
    );
};

   RCStore gStore;  

   Wb::IOstat gstatR,gstatM,gstatC,gstatX;

   #define RC_MDIR 1 
   #define RC_SYNC 2
   #define RC_CURR 4 

   #define RC_SAVE   (RC_MDIR|RC_CURR)
   #define RC_SAVE_S (RC_MDIR|RC_SYNC)

   #define RCL_V1  1  
   #define RCL_V2  2  

   #define RCL_ALL (RCL_V1|RCL_V2)

int RCStore::setup_root_directory(
   const char *F, int L, const char *env, Wb::Path &root, char append
){
   int m=0, k=0, e=0, kmax; 
   const char *sval;

   if (!env || !env[0]) wblog(F_L,
      "ERR %s() no variable specified",FCT);

   sval=getenv(env); if (!sval || !sval[0]) {
      if (!append) { root.init(); }
      return m;
   }

   if (!append)
        { m=root.init(F_L,sval,":",1024); }
   else { k=root.numel(); m=root.append(F_L,sval,":"); }

   for (kmax=k+m; k<kmax; ++k) { char *dir=root[k];

      if (!dir || !dir[0]) {
         if (!e) { PRINTF("\n"); }; e|=1;
         PRINTF("   ERR invalid %s entry (empty)\n",env);
      }
      if (Wb::fexist(dir,'d')) { continue; }

      if (k+1<kmax) {
         if (!e) { PRINTF("\n"); }; e|=2;
         PRINTF("   ERR invalid %s entry (%s) ^1\n",env,dir);
      }

      unsigned n=strlen(dir), i=n-1, i0=0; int q;

      for (; i<n; --i) { if (dir[i]=='/') { i0=i; break; }}
      if (i0 && (
         !strncasecmp(dir+i0+1,"RCStore",7) ||
         !strncasecmp(dir+i0+1,"RCSync", 6))
      ){
         dir[i0]=0; q=Wb::fexist(dir,'d');
         dir[i0]='/';
         if (q) {
            mode_t p=(S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            if ((q=mkdir(dir,p))==0) { 
               PRINTF("   creating %-8s root directory (%s)\n",
               env, Wb::repHome(dir).data);
            }
            else {
               if (!e) { PRINTF("\n"); }; e|=4;
               PRINTF("   ERR invalid %s entry (%s) ^2\n",env,dir);
            }
         }
         else { 
            if (!e) { PRINTF("\n"); }; e|=8;
            PRINTF("   ERR invalid %s entry (%s) ^3\n",env,dir);
         }
      }
      else if (!Wb::fexist(dir,'d')) {
         if (!e) { PRINTF("\n"); }; e|=16;
         PRINTF("   ERR invalid %s entry (%s) ^4\n",env,dir);
      }
   }
   if (e) { PRINTF("\n");
      if (e& 2) { PRINTF("   ^1 all but last entry must exist\n"); }
      if (e& 4) { PRINTF("   ^2 failed to create directory\n"); }
      if (e& 8) { PRINTF("   ^3 invalid %sdirectory\n",kmax>1?"last ":""); }
      if (e&16) { PRINTF("   ^4 %sdirectory must also exist if not named %s*\n",
         strncmp(env,"RC_SYNC",8) ? "RCStore":"RCSync", kmax>1?"last ":""); }
      if (e> 1) { PRINTF("\n"); }
   }

   if (e) {
      root.init(); 
      wblog(FL,"ERR invalid env %s directory (e=%d)",env,e);
   }

   return (e ? -m : m);
};

namespace CG {

template <class TQ, class TD>
double getSymmetryStates(const char *F, int L, const QType &q, 
   const wbvector< SPARRAY_TD > &Sp,
   const wbvector< SPARRAY_TD > &Sz,
   wbvector< SPARRAY_TD > &U, 
   wbvector<unsigned> &dd, 
   wbvector<genRG_base<TQ,TD> > &RR,
   wbvector<double> &c2eps, 
   wbMatrix<unsigned> *iOM=NULL, 
   char vflag=1
);

template <class TQ>
SPIDX_T findMaxWeight( 
   const QType &q, const wbMatrix<double> &Z, qset<TQ> *J=NULL, wbperm *P_=NULL
);

template <class TQ, class TD>
void get_SU2mat(TQ s2,   
   wbvector<double> &sz, 
   wbsparray<TD> &Sp, wbsparray<TD> &Sz,
   wbsparray<TD> &S2, wbsparray<TD> &E
);

int isQSet (const mxArray *S, unsigned k=-1);
int isCData(const mxArray *S, unsigned k=-1);
int isCRef (const mxArray *S, unsigned k=-1);

int MemStat(const char *F=0, int L=0) { 
   int n=0; const char *hstr="RCStore"; const char *xstr="/full";

   wblog(F_L,"%N%s() %33R",FCT,"-");
      gRS.Info();
      gCS.Info();
      gXS.Info();
   Wb::MemStat(FL); 

   if (!gstatR.isEmpty()) gstatR.print("RData", (++n==1)? hstr:0, xstr);
   if (!gstatM.isEmpty()) gstatM.print("map3",  (++n==1)? hstr:0, xstr);
   if (!gstatC.isEmpty()) gstatC.print("CData", (++n==1)? hstr:0, xstr);
   if (!gstatX.isEmpty()) gstatX.print("XData", (++n==1)? hstr:0, xstr);

   if (n) { PRINTF("\n"); }
   else wblog(FL,"I/O no data read or written yet to RCStore");

   return 0;
};

#ifdef QS_USING_OMP

   map <size_t, Wb::ompNLock*> lock_map;    

class Guard { 

  public:

     Guard() : active(false) {}; 

     template<class TQ>
     Guard(const char *F, int L, const QSet<TQ> &Q, const char *tag=NULL,
        const char *cond=NULL)
      : active(false) { acquire(F_L,Q,tag,cond); };

     template<class TM>
     Guard(const char *F, int L, const cgc_contract_id<TM> &x,
        const char *cond=NULL)
      : active(false) { acquire(F_L,x,cond); }

     template<class TQ>
     Guard(const char *F, int L, const QSet<TQ> &QA, const QSet<TQ> &QB)
      : active(false) { acquire(F_L,QA,QB); };

     template<class TQ>
     Guard(const char *F, int L,
        const QSet<TQ> &QA, const QSet<TQ> &QB, const QSet<TQ> &QC)
      : active(false) { acquire(F_L,QA,QB,QC); };

    ~Guard() {
        if (hid.len) { release(); hid.init(); }
        active=false;
     };

     explicit operator bool() const { return (hid.len || active); }
     bool operator! () const { return (!hid.len && !active); }

     int acquire(const char *F, int L, 
        const FileLock &flk, const char *tag=NULL, const char *cond=NULL);

     template<class TQ> 
     int acquire(const char *F, int L,
        const QSet<TQ> &Q, const char *tag=NULL, const char *cond=NULL);

     template<class TM> 
     int acquire(const char *F, int L,
        const cgc_contract_id<TM> &x, const char *cond=NULL);

     template<class TQ> 
     void acquire(const char *F, int L, const QSet<TQ>* Q[], unsigned n);

     template<class TQ>
     void acquire(const char *F, int L,
        const QSet<TQ> &QA, const QSet<TQ> &QB)
      { const QSet<TQ>* QQ[2]={ &QA, &QB }; acquire(F_L,QQ,2); };

     template<class TQ>
     void acquire(const char *F, int L,
        const QSet<TQ> &QA, const QSet<TQ> &QB, const QSet<TQ> &QC)
      { const QSet<TQ>* QQ[3]={ &QA, &QB, &QC }; acquire(F_L,QQ,3); };

     inline int acquire_set(const char *ss[], int ntry=-1) { 
        if (!hid.len) { return -1; }

        wbvec<int> got(hid.len);
        wbvec<Wb::ompNLock*> lk(hid.len);

        int nlks=acquire_set_iter(got.data,lk.data,ss,0);
        if (nlks<=0) { if (ntry==1 || !ntry) return nlks;
           nlks=acquire_set_wait(got.data,lk.data,ss,ntry);
        }

        if (nlks<=0) {
           if (ntry<=0 || ntry>99) wblog(FL,
             "ERR %s() failed to qdquire locks (%d/%d)",FCT,nlks,hid.len);
           return nlks;
        }

        for (unsigned i=0; i<hid.len; ++i) {
            lk[i]->set_istr(ss[i]);
            if (got[0]<got[i]) { got[0]=got[i]; } 
        }
        active=true;

        return got[0];
     };

     int acquire_set_iter( 
         int *got, Wb::ompNLock *lk[], const char *ss[], unsigned iter=0);
     int acquire_set_wait( 
         int *got, Wb::ompNLock *lk[], const char *ss[], int ntry=-1);

     void release();

     int level() { int l=-1;
        for (unsigned i=0; i<hid.len; ++i) {
            const Wb::ompNLock *lp = Guard::find(hid[i]);
            if (lp && l<(lp->level)) { l=lp->level; }
        }
        return l;
     };

     int used() { int l=0; 
        for (unsigned i=0; i<hid.len; ++i) {
            const Wb::ompNLock *lp = Guard::find(hid[i]);
            if (lp && l<(lp->level)) {
               if (l< lp->level  ) { l=lp->level;   }
               if (l< lp->nactive) { l=lp->nactive; }
               if (l< lp->ntot   ) { l=lp->ntot;    }
            }
        }
        if (!l && active) { l=-1; }
        return l;
     };

     void print(const char *F, int L);

     wbstring toStr() const;

     static wbstring toStr(size_t h, const char *istr_=NULL);
     static int toStr(wbstring &s, size_t h, const char *istr_=NULL);

     template<class TQ>
     static wbstring Status(const QSet<TQ> &Q, const char *tag=NULL);

     template<class TM>
     static wbstring Status(const cgc_contract_id<TM> &x);

     static wbstring Status(size_t h1, const char *istr=NULL);

     static Wb::ompNLock* find(size_t h);

     static int print_CG_locks(const char *F=NULL, int L=0);

  protected:
  private:

     wbvector<size_t> hid;

     bool active;

     const char* sprintf_locks(
        wbstring &sout, const int* got, const char *istr[],
        const char *indent=NULL);

     int deal_with_cond(const char *F, int L, size_t h0, const char *cond);

}; 

#endif 

class FileLock { 

  public:

     FileLock() : fid(0) {}; 

    ~FileLock() {
        if (fid) { ReleaseLock(FL); }
     };

     template<class TQ>
     FileLock(
        const QType &q_, const qset<TQ> &J1, const qset<TQ> &J2,
        const char *ext="mp3" 
     ){ init(q_,J1,J2,ext); }

     template<class TQ>
     int init(
        const QType &q_, const qset<TQ> &J1, const qset<TQ> &J2,
        const char *ext="mp3" 
     );

     int GetLock(const char *F=NULL, int L=0, char flag='W') {
        return tryGetLock(F,L,flag); 
     };

     int ReleaseLock(const char *F=NULL, int L=0) {
        return tryGetLock(F,L,'u' ); 
     };

     int ensureFOpen(const char *F=NULL, int L=0) const {
        int q=0;
        if (!fname) { q|=1; }
        if (fid<3 ) { q|=2; }
        if (q && F) {
           if (q&1) wblog(F_L,"ERR %s() got empty fname (e=%d)",FCT,q);
           if (q&2) wblog(F_L,"ERR %s() file not yet open "
              "(fid=%d; e=%d)\n%s",FCT,fid, q,fname.data
           );
        }
        return q;
     };

     int fid;

     QType q; 
     wbstring fname;

#ifdef QS_USING_OMP
     CG::Guard fLK; 
#endif

  protected:
  private:

     int tryGetLock(const char *F, int L, char flag);
};

}; 

template<class TD>
int BuildKrylovH(
   const QType &q,
   wbarray<TD> &HK, unsigned nk, 
   SPARRAY_TD &X,        
   wbvector<TD> &EK,
   const SPARRAY_TD &HL,            
   const SPARRAY_TD &HR,            
   const wbvector< SPARRAY_TD > SL, 
   const wbvector< SPARRAY_TD > SR,
   double &xmin
);

template<class TD>
void GetHPsi(
   SPARRAY_TD &HX, const SPARRAY_TD &X,
   const SPARRAY_TD &HL, 
   const SPARRAY_TD &HR, 
   const wbvector< SPARRAY_TD > SL, 
   const wbvector< SPARRAY_TD > SR  
);

#endif

