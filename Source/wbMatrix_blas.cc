/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbMatrix (matrix class, row-major)
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

#ifndef __WB_BLAS_CC__
#define __WB_BLAS_CC__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //
// wrapper routine

template<class TA, class TB, class TC>
void MMDIAG(
    const wbMatrix<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);

template<class TA, class TB, class TC>
wbMatrix<TC>& MMDIAG(
    const wbvector<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);

template<class TA, class TB, class TC>
wbMatrix<TC>& MMDIAG(
    const wbMatrix<TA> &A,
    const wbvector<TB> &B, wbMatrix<TC> &C,
    char aflag='N', char bflag='N',
    const TA afac=1., const TC cfac=0.
);

template<class TA, class TB, class TC>
void DZGEMM_VEC(
    const wbMatrix<TA> &A,
    const wbvector<TB> &B,
    wbvector<TC> &C, int m, int n,
    char aflag='N', const TA afac=1., const TC cfac=0.
);

template<class TD> 
void DZGEMM(
    const wbMatrix<TD> &A,
    const wbMatrix<TD> &B, wbMatrix<TD> &C,
    const unsigned k, 
    char aflag='N', char bflag='N',
    const TD afac=1., const TD cfac=0.
);

template<>
inline void DZGEMM<double>(
    const wbMatrix<double> &A,
    const wbMatrix<double> &B, wbMatrix<double> &C,
    const unsigned k, 
    char aflag, char bflag,
    const double afac, const double cfac
){

#ifdef WB_CLOCK
   wbc_dgemm.resume();

   if (aflag=='N')
        stat_dgemm.account(A.dim1*A.dim2*(A.dim2==B.dim1 ? B.dim2 : B.dim1));
   else stat_dgemm.account(A.dim1*A.dim2*(A.dim1==B.dim1 ? B.dim2 : B.dim1));
#endif

    dgemm(bflag, aflag, int(C.dim2), int(C.dim1), int(k),
       afac, B.data, int(B.dim2),
             A.data, int(A.dim2),
       cfac, C.data, int(C.dim2)
    );

#ifdef WB_CLOCK
   wbc_dgemm.stop();
#endif
}

template<>
inline void DZGEMM<wbcomplex>(
    const wbMatrix<wbcomplex> &A,
    const wbMatrix<wbcomplex> &B, wbMatrix<wbcomplex> &C,
    const unsigned k, 
    char aflag, char bflag,
    const wbcomplex afac, const wbcomplex cfac
){
#ifdef WB_CLOCK
    wbc_zgemm.resume();

    if (bflag=='N')
         stat_dgemm.account(A.dim1*A.dim2*(A.dim2==B.dim1 ? B.dim2 : B.dim1));
    else stat_dgemm.account(A.dim1*A.dim2*(A.dim1==B.dim1 ? B.dim2 : B.dim1));
#endif

    zgemm (
      bflag, aflag, (int)C.dim2, (int)C.dim1, (int)k,
      afac, B.data, (int)B.dim2, A.data, (int)A.dim2,
      cfac, C.data, (int)C.dim2
    );

#ifdef WB_CLOCK
    if (wbc_zgemm.stop()) {
       wbc_zgeNX.flag+=2; 
       wbc_zgeNX.tcpu+= (long unsigned)
       (double(A.dim1*A.dim2*B.dim2)*(4E-9*CLOCKS_PER_SEC)+1.5);
    }
    else {
       wbc_zgeNN.flag+=2; 
       wbc_zgeNN.tcpu+= (long unsigned)
       (double(A.dim1*A.dim2*B.dim2)*(4E-9*CLOCKS_PER_SEC)+1.5);
    }
#endif
}

namespace Wb {

void MatProd( 
    const wbMatrix<double> &A,
    const wbMatrix<wbcomplex> &B, wbMatrix<wbcomplex> &C,
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.,
    const char i00flag=0
){
    wbMatrix<double> R,I, rC, iC;

    B.getReal(R);
    B.getImag(I);

    if (bflag=='C') I.flipSign();

    if (cfac!=0.) { C.getReal(rC); C.getImag(iC); }

    MatProd(A, R, rC, aflag, bflag, afac, cfac, i00flag);
    MatProd(A, I, iC, aflag, bflag, afac, cfac, i00flag);

    C.set(rC,iC);
};

void MatProd( 
    const wbMatrix<wbcomplex> &A,
    const wbMatrix<double> &B, wbMatrix<wbcomplex> &C,
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.,
    const char i00flag=0
){
    wbMatrix<double> R,I, rC, iC;

    A.getReal(R);
    A.getImag(I);

    if (aflag=='C') I.flipSign();

    if (cfac!=0.) { C.getReal(rC); C.getImag(iC); }

    MatProd(R, B, rC, aflag, bflag, afac, cfac, i00flag);
    MatProd(I, B, iC, aflag, bflag, afac, cfac, i00flag);

    C.set(rC,iC);
};

template<class T>  
void MatProd( 
    const wbMatrix<T> &A,
    const wbMatrix<T> &B, wbMatrix<T> &C,
    char aflag, char bflag,
    const T afac, const T cfac,
    const char i00flag
){
    unsigned a1=A.dim1, a2=A.dim2, b1=B.dim1, b2=B.dim2;
    const char flags[]="NTCntc";

#ifdef WB_CLOCK
    wbc_matprod.resume();
#endif

    if (!strchr(flags,aflag) || !strchr(flags,bflag)) wblog(FL,
       "ERR %s() invalid flags %c<%d>, %c<%d>",
       FCT,aflag,aflag,bflag,bflag);

    if (aflag!='N') SWAP(a1,a2);
    if (bflag!='N') SWAP(b1,b2);

    if (a2!=b1) {
        wblog(FL, "ERR %s() dimension mismatch: (%d,%d) * (%d,%d) ?",
        FCT,a1,a2,b1,b2); return;
    }

    if (cfac!=0.) {
        if (C.data==NULL && i00flag) {
           if (cfac!=1.) wblog(FL,
              "WRN C = A*B + c*[] with c=%s !?", toStr(cfac).data);
           C.init(a1,b2);
        }
        else if (C.dim1!=a1 || C.dim2!=b2) { wblog(FL,
           "ERR %s() dimension mismatch: C=(%d,%d) =? (%d,%d)",
            FCT,C.dim1,C.dim2,a1,b2); return;
        }
    }

    if (a1==0 || b2==0) {
       if (cfac==0.) C.init(a1,b2);
       return;
    }

    if (a2==0) { wblog(FL, 
       "ERR %s() cannot multiply (%dx%d)*(%dx%d)!",
       FCT,a1,a2,b1,b2); return; }

    if (&C==&A || &C==&B) { wblog(FL, 
       "ERR I/O spaces must be distinct!\n[%7lX %7lX %7lX]",
       &A, &B, &C); return; }

    if (cfac==0.)
    C.init(a1,b2);

    if (A.isdiag || B.isdiag)
         MMDIAG(A,B,C,   aflag,bflag,afac,cfac);
    else DZGEMM(A,B,C,a2,aflag,bflag,afac,cfac);

#ifdef WB_CLOCK
    wbc_matprod.stop();
#endif
};

template<class TA, class TB, class TC> 
void MatProd(
    const wbMatrix<TA> &A,
    const wbvector<TB> &B, wbvector<TC> &C,
    char aflag, const TC afac, const TC cfac,
    const char i00flag 
){
    unsigned m=A.dim1, n=A.dim2;
    if (aflag!='N') SWAP(m,n);

    if (B.len!=n) wblog(FL,
       "ERR severe size mismatch Ax=b (%dx%d, %d)",m,n,B.len);

    if (!C.data) {
       if (cfac!=0. && i00flag) wblog(FL,
          "WRN C = A*B + c*[] with c=%s !?",toStr(cfac).data);
       C.init(m);
    }
    else if (C.len!=m) {
       if (cfac==0.) C.init(m); else wblog(FL,
          "ERR severe size mismatch Ax=b (%d, %d)", m, C.len);
    }

    DZGEMM_VEC(A,B,C,m,n,aflag,afac,cfac);
};

}; 

template<>
inline void DZGEMM_VEC(
    const wbMatrix<double> &A,
    const wbvector<double> &B,
    wbvector<double> &C, int m, int n,
    char aflag, const double afac, const double cfac
){
    dgemm (
      'N', aflag, 1, m, n, afac,
       B.data, 1, A.data, A.dim2, cfac,
       C.data, 1
    );
}

template<>
inline void DZGEMM_VEC(
    const wbMatrix<wbcomplex> &A,
    const wbvector<wbcomplex> &B,
    wbvector<wbcomplex> &C, int m, int n,
    char aflag, const wbcomplex afac, const wbcomplex cfac
){
    zgemm (
      'N', aflag, 1, m, n, afac,
       B.data, 1, A.data, A.dim2, cfac,
       C.data, 1
    );
}

template<>
inline void DZGEMM_VEC(
    const wbMatrix<double> &A,
    const wbvector<wbcomplex> &B,
    wbvector<wbcomplex> &C, int m, int n,
    char aflag, const double afac, const wbcomplex cfac
){
    wbvector<double> b,cr,ci;

    if (cfac && cfac!=1.) { C*=cfac; }

    B.getReal(b); C.getReal(cr); DZGEMM_VEC(A,b,cr,m,n,aflag,afac,1.);
    B.getImag(b); C.getImag(ci); DZGEMM_VEC(A,b,ci,m,n,aflag,afac,1.);

    C.set(cr,ci);
};

inline void DZGEMM_aux(
    const wbMatrix<double> &A,
    const wbMatrix<double> &B, wbMatrix<double> &C,
    const unsigned k, 
    char aflag='N', char bflag='N',
    const double afac=1., const double cfac=0.
){
    if (A.isdiag || B.isdiag)
         MMDIAG(A, B, C,    aflag, bflag, afac, cfac);
    else DZGEMM(A, B, C, k, aflag, bflag, afac, cfac);
}

template<class TA, class TB, class TC>
inline void MMDIAG(
    const wbMatrix<TA> &A,
    const wbMatrix<TB> &B, wbMatrix<TC> &C,
    char aflag, char bflag,
    const TA afac, const TC cfac
){
    unsigned i, s=C.dim1*C.dim2, dmax;

    if (!A.isdiag && !B.isdiag) wblog(FL,
    "ERR Wrong call - elements not diagonal.");

    if (A.isdiag && A.dim1>1 && A.dim2>1)  
    if (A(1,0)!=0. || A(0,1)!=0.) wblog(FL,
    "ERR Input A is NOT diag even though isdiag=%d", A.isdiag);

    if (B.isdiag && B.dim1>1 && B.dim2>1)  
    if (B(1,0)!=0. || B(0,1)!=0.) wblog(FL,
    "ERR Input B is NOT diag even though isdiag=%d", B.isdiag);

    if (cfac!=1.) for (i=0; i<s; i++) C.data[i]*=cfac;

    if (A.isdiag && !B.isdiag) {
        MMDIAG(A.getDiag(), B, C, aflag, bflag, afac);
    }
    else
    if ((!A.isdiag) && B.isdiag) {
        MMDIAG(A, B.getDiag(), C, aflag, bflag, afac);
    }
    else { 
        wbvector<TA> adiag;
        wbvector<TB> bdiag;
        A.getDiag(adiag); B.getDiag(bdiag);

        if (aflag=='C')
        for (i=0; i<adiag.len; i++) adiag[i]=conj(adiag[i]);

        if (bflag=='C')
        for (i=0; i<bdiag.len; i++) bdiag[i]=conj(bdiag[i]);

        dmax=MIN( MIN(A.dim1,A.dim2), MIN(B.dim1,B.dim2) );
        for (i=0; i<dmax; i++) C(i,i) += afac * adiag[i] * bdiag[i];
    }

    if (C.isdiag) if (!A.isdiag || !B.isdiag) C.isdiag=0;
}

template<class TA, class TB, class TC>
inline wbMatrix<TC>& MMDIAG(
   const wbvector<TA> &A,
   const wbMatrix<TB> &B, wbMatrix<TC> &C,
   char aflag, char bflag,
   const TA afac, const TC cfac
){
   unsigned i,j, m=B.dim1, n=B.dim2;
   char acc,btr;
   TA dbl;

   if (bflag=='C' && typeid(TB)!=typeid(wbcomplex)) bflag='T';

   if (!strchr("NTC",aflag) || !strchr("NTC",bflag)) wblog(FL,
   "ERR Invalid flag(s) %c<%d>, %c<%d>",aflag,aflag,bflag,bflag);

   acc=(aflag=='C');
   btr=(bflag!='N'); if (btr) SWAP(m,n);

   if (A.len!=m) wblog(FL,
   "ERR Size mismatch (diag(%d) * %dx%d)", A.len,m,n);
   if (&C==&B) wblog(FL,"ERR Output space C equal to input space.");

   if (C.isEmpty() || (cfac==0. && (C.dim1!=m || C.dim2!=n)))
      C.init(m,n);
   else {
      if (C.dim1!=m || C.dim2!=n) wblog(FL,
      "ERR Size mismatch (%dx%d; %dx%d)",C.dim1,C.dim2, m,n);
      if (cfac!=1.) {
         if (cfac==0.) C.reset();
         else C*=cfac;
      }
   }

   for (i=0; i<C.dim1; i++) {
      dbl=A[i];   if (acc) dbl=conj(dbl);
      dbl*=afac;  if (dbl==0.) continue;

      switch (bflag) {
        case 'N': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*B(i,j)); break;
        case 'T': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*B(j,i)); break;
        case 'C': for (j=0; j<C.dim2; j++) C(i,j)+=(dbl*conj(B(j,i)));
      }
   }

   return C;
}

template<class TA, class TB, class TC>
inline wbMatrix<TC>& MMDIAG(
   const wbMatrix<TA> &A,
   const wbvector<TB> &B, wbMatrix<TC> &C,
   char aflag, char bflag,
   const TA afac, const TC cfac
){
   unsigned i,j, m=A.dim1, n=A.dim2;
   char atr,bcc;
   TB dbl;

   if (aflag=='C' && typeid(TA)!=typeid(wbcomplex)) aflag='T';

   if (!strchr("NTC",aflag) || !strchr("NTC",bflag)) wblog(FL,
   "ERR Invalid flag(s) %c<%d>, %c<%d>",aflag,aflag,bflag,bflag);

   atr=(aflag!='N'); if (atr) SWAP(m,n);
   bcc=(bflag=='C');

   if (B.len!=n) wblog(FL,
      "ERR Size mismatch (%dx%d * diag(%d))", m,n, B.len);
   if (&C==&A) wblog(FL,"ERR Output space C equal to input space.");

   if (C.isEmpty() || (cfac==0. && (C.dim1!=m || C.dim2!=n)))
      C.init(m,n);
   else {
      if (C.dim1!=m || C.dim2!=n) wblog(FL,
      "ERR Size mismatch (%dx%d; %dx%d)",C.dim1,C.dim2, m,n);
      if (cfac!=1.) {
         if (cfac==0.) C.reset();
         else C*=cfac;
      }
   }

   for (j=0; j<C.dim2; j++) {
      dbl=B[j];   if (bcc) dbl=conj(dbl);
      dbl*=afac;  if (dbl==0.) continue;

      switch (aflag) {
        case 'N': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*A(i,j)); break;
        case 'T': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*A(j,i)); break;
        case 'C': for (i=0; i<C.dim1; i++) C(i,j)+=(dbl*conj(A(j,i)));
      }
   }

   return C;
}

inline void WbEigenSymmetric (
    const wbMatrix<double> &M, wbMatrix<double> &V, wbvector<double> &E
){
    unsigned n=M.dim1; pINT ni=0, q=0;
    double nd=0;

    if (M.dim1!=M.dim2) wblog(FL,
       "ERR WbEigenSymmetric requires square matrix (%d,%d).",
        M.dim1, M.dim2);

    V=M; E.init(n); if (n==0) return;

    dsyevd('V','U',n,V.data,n,E.data,&nd,-1,&ni,-1,q); 
    if (nd<1 || ni<1 || q) wblog(FL,
       "ERR DSYEVD returned lwork=%g/%d (e=%d) !?",nd,ni,q);

    wbvector<double> wd(nd);
    wbvector<pINT> wi(ni);

    dsyevd(
       'V',        
       'L',        
        n,         
        V.data,    
        n,         
        E.data,    
        wd.data, wd.len, 
        wi.data, wi.len, 
        q          
    );

    if (q) wblog(FL,"ERR DSYEV returned e=%d !?", q);
}

inline void WbEigenSymmetric (
    const wbMatrix<wbcomplex> &M, wbMatrix<wbcomplex> &V, wbvector<double> &E
){
    unsigned n=M.dim1; pINT q, ni=0;
    wbcomplex nz=0; double nr=0;

    if (M.dim1!=M.dim2) wblog(FL,
       "ERR %s() requires square matrix (%d,%d)",FCT,M.dim1,M.dim2);
    if (!M.isHConj(1E-12,&nr)) wblog(FL,
       "ERR %s() got non-hermitian matrix (%.3g)",FCT,nr);
    nr=0;

    V=M; E.init(n); if (!n) return;

    zheevd('V','U',n,V.data,n,E.data,&nz,-1,&nr,-1,&ni,-1,q); 
    if (nz.r<1 || nr<1 || ni<1 || q) wblog(FL,
       "ERR ZHEEVD returned lwork=%g/%g/%d (e=%d) !?",nz.r,nr,ni,q);

    wbvector<wbcomplex> wz(nz.r);
    wbvector<double> wr(nr);
    wbvector<pINT> wi(ni);

    zheevd(
       'V',        
       'L',        
        n,         
        V.data,    
        n,         
        E.data,    
        wz.data, wz.len, 
        wr.data, wr.len, 
        wi.data, wi.len, 
        q          
    );

    if (q<0) { wblog(FL,"ERR %s() ZHEEVD returned e=%d",FCT,q); }
    if (q>0) {
       V=M; E.init(n);

       wblog(FL,"WRN %s() ZHEEVD returned e=%d -> fall back to ZHEEV",FCT,q);
       MXPut(FL,0,"tmpfile").add(M,"M").add(q,"q");

       if (n<128) { ni=2*n-1; } 
       else { 
          ni=ilaenv(1,"zhetrd","U",(pINT)n,(pINT)n,(pINT)n,(pINT)n);
          if (ni>n) wblog(FL,"WRN %s() got lwork=(%d+2)*%d",FCT,ni,n);
          ni=(ni+1)*n;
       }
       if (size_t(ni)>wz.len) { wz.init(ni); }

       nr=3*n-2; if (nr>wr.len) { wr.init(nr); }

       zheev(
          'V',      
          'L',      
           n,       
           V.data,  
           n,       
           E.data,  
           wz.data, ni,    
           wr.data, 
           q        
       );
       wblog(FL,"ERR ZHEEV returned e=%d", q);
    }

    V.Conj();
};

template<class T>
inline void GESVD_M(  
   wbMatrix<T> &A,
   wbMatrix<T> &U,
   wbvector<double> &S,
   wbMatrix<T> &Vt
){
   wblog(FL,"ERR GESVD_M not defined for datatype %s",TSTR(T));
};

template<>
inline void GESVD_M(
   wbMatrix<double> &A, 
   wbMatrix<double> &U, 
   wbvector<double> &S,
   wbMatrix<double> &Vt
){
   unsigned M=A.dim1, N=A.dim2, K=S.len;
   wbvector<double> wd;
   double nd;
   pINT q, ni=8*MIN(M,N); 
   wbvector<pINT> wi(ni);

#ifdef WB_CLOCK
    wbc_dgesvd.resume();
#endif

   dgesdd('S',N,M,A.data,N,S.data,Vt.data,N,U.data,K,&nd,-1,wi.data,q);
   if (nd<1 || q) wblog(FL,
      "ERR DGESDD returned lwork=%g (e=%d) !?",nd,q);
   wd.init(nd);

   dgesdd(
     'S', N, M, A.data, N, S.data,
      Vt.data, N,   
      U.data, K,    
      wd.data, wd.len, wi.data,
      q
   );

#ifdef WB_CLOCK
    wbc_dgesvd.stop();
#endif

   if (q) wblog(FL,"ERR DGESDD returned e=%d !?", q);
};

template<>
inline void GESVD_M(
   wbMatrix<wbcomplex> &A, 
   wbMatrix<wbcomplex> &U, 
   wbvector<double> &S,
   wbMatrix<wbcomplex> &Vt
){
   unsigned M=A.dim1, N=A.dim2, K=S.len;
   wbvector<wbcomplex> wz;

   pINT q, n=MIN(M,N);
   wbcomplex nz;

   wbvector<double> wd(n*(5*n+7)); 
   wbvector<pINT> wi(12*n); 

   if (!M || !N) wblog(FL,"ERR %s() got %dx%d matrix!?",FCT,M,N);

#ifdef WB_CLOCK
    wbc_zgesvd.resume();
#endif

   zgesdd(
     'S',N,M,A.data,N,S.data,Vt.data,N,U.data,K,
      &nz, -1, wd.data, wi.data, q);
   if (nz.r<1 || q) wblog(FL,
      "ERR ZGESDD() returned lwork=%g (e=%d) !?",nz.r,q);
   wz.init(nz.r);

   zgesdd(
     'S', N, M, A.data, N, S.data, Vt.data, N, U.data, K,
      wz.data, wz.len, wd.data, wi.data,
      q
   );

#ifdef WB_CLOCK
    wbc_zgesvd.stop();
#endif

   if (q) wblog(FL,"ERR ZGESDD returned e=%d !?", q);
};

template<class T>
inline void wbSVD(
   const wbMatrix<T> &A0,
   wbMatrix<T> &U,
   wbvector<double> &S,
   wbMatrix<T> &Vt 
){
   wbMatrix<T> A(A0);

   unsigned M=A.dim1, N=A.dim2, K=MIN(M,N);

   U.init(M,K); S.init(K); Vt.init(K,N);
   if (!M || !N) return;

   GESVD_M(A,U,S,Vt);

#ifdef SVD_BUG_SAFEGUARD 
   { double a=SQRT(A0.norm2()), s=SQRT(S.norm2()), e=fabs(a-s)/a;
     if (e>1E-10) {
        wblog(FL,"WRN %s() need to fix SVD bug: %dx%d @ e=%.3g !?",FCT,M,N,e);
     }
   }
#endif
};

#endif

