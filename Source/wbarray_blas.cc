/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbarray (array class, col-major)
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

#ifndef __WB_ARRAY_BLAS_CC__
#define __WB_ARRAY_BLAS_CC__

/* ================================================================= */
// NB! call DZGEMM after all checks are made (such as size, rank, ...)
// C = (alpha=afac) * op(A)*op(B) + (beta=cfac)*C
// e.g., also C must already be fully allocated // DZGEMM_C_ALLOCATED
/* ================================================================= */

template<class TA, class TB, class TC > 
void DZGEMM(
   const wbarray<TA> &A,
   const wbarray<TB> &B, wbarray<TC> &C,
   const unsigned k,  
   char aflag='N', char bflag='N', 
   TA afac=1, TC cfac=0
);

#ifdef QS_USING_MPFR

template<> inline
void DZGEMM(
   const wbarray<Wb::quad> &A,
   const wbarray<Wb::quad> &B, wbarray<Wb::quad> &C,
   const unsigned k, char aflag, char bflag,
   Wb::quad afac, Wb::quad cfac
){
   wblog(FL,"ERR %s() not yet defined for data type Wb::quad");
   wblog(FL,"ERR %s() %p %p %p, %d %d %d, %p %p",
   FCT,&A, &B, &C, k, aflag, bflag, &afac, &cfac);
};

#endif

template<> inline
void DZGEMM(
    const wbarray<double> &A,
    const wbarray<double> &B, wbarray<double> &C,
    const unsigned k,  
    char aflag, char bflag,
    double afac, double cfac
){
   if (A.SIZE.len!=2 || B.SIZE.len!=2 || C.SIZE.len!=2) wblog(FL,
      "ERR %s() requires matrices (got: %s * %s = %s)",
      FCT, SSTR(A), SSTR(B), SSTR(C));

   unsigned
      a1=A.SIZE[0],
      b1=B.SIZE[0],
      c1=C.SIZE[0], c2=C.SIZE[1]; 

   if (!a1 || !b1 || !c1 || !c2) {
      if (A.SIZE[1]==b1 && B.SIZE[1]==c2 && a1==c1) { 
         if (cfac!=double(1)) C*=cfac;
         return;
      }
      wblog(FL,"ERR %s() got empty matrices (A: %s, B: %s; %c%c)",
      FCT,SSTR(A),SSTR(B),aflag,bflag); return;
   }

#ifdef WB_CLOCK
   unsigned a2=A.SIZE[1], b2=B.SIZE[1];

   stat_dgemm.account(a1*a2*(a2==b1 ? b2 : b1));

   Wb::Clock clk("dgemm",1); 
#endif

   dgemm(
      aflag, bflag, (pINT)c1, (pINT)c2, (pINT)k,
      afac, A.data, (pINT)a1, B.data, (pINT)b1,
      cfac, C.data, (pINT)c1
   );
};

template<> inline
void DZGEMM(
    const wbarray<wbcomplex> &A,
    const wbarray<wbcomplex> &B, wbarray<wbcomplex> &C,
    const unsigned k,  
    char aflag, char bflag,
    wbcomplex afac, wbcomplex cfac
){
   if (A.SIZE.len!=2 || B.SIZE.len!=2) wblog(FL,
      "ERR %s() matrices required (%s * %s)",FCT,
   A.sizeStr().data, B.sizeStr().data);

   unsigned
      a1=A.SIZE[0],
      b1=B.SIZE[0],
      c1=C.SIZE[0], c2=C.SIZE[1]; 

   if (!a1 || !b1 || !c1 || !c2) {
      if (A.SIZE[1]==b1 && B.SIZE[1]==c2 && a1==c1) { 
         if (cfac!=double(1)) C*=cfac;
         return;
      }
      wblog(FL,"WRN %s() got empty matrices (A: %s, B: %s)",
      FCT, A.sizeStr().data, B.sizeStr().data); return;
   }

#ifdef WB_CLOCK
   unsigned a2=A.SIZE[1], b2=B.SIZE[1];

   stat_dgemm.account(a1*a2*(a2==b1 ? b2 : b1));
   Wb::Clock clk("zgemm",1); 
#endif

   zgemm (
      aflag, bflag, (pINT)c1, (pINT)c2, (pINT)k,
      afac, A.data, (pINT)a1, B.data, (pINT)b1,
      cfac, C.data, (pINT)c1
   );

#ifdef WB_CLOCK
 { Wb::Clock *clk;
   if (clk.stop())
        { Wb::Clock("zgeNX",1,0,&Wb::Clocks,&clk,0); } 
   else { Wb::Clock("zgeNN",1,0,&Wb::Clocks,&clk,0); } 

   clk.ncall+=2; 
   clk.tcpu+= (long unsigned)
   (double(a1*a2*b2)*(4E-9*CLOCKS_PER_SEC)+1.5);
 }
#endif
};

template<> inline
void DZGEMM(
    const wbarray<wbcomplex> &A,
    const wbarray<double> &B0, wbarray<wbcomplex> &C,
    const unsigned k,  
    char aflag, char bflag,
    wbcomplex afac, wbcomplex cfac
){
    wbarray<wbcomplex> B; B.set(B0, wbarray<double>());
    DZGEMM(A,B,C,k,aflag,bflag,afac,cfac);
};

template<> inline
void DZGEMM(
    const wbarray<double> &A0,
    const wbarray<wbcomplex> &B, wbarray<wbcomplex> &C,
    const unsigned k,  
    char aflag, char bflag,
    double afac, wbcomplex cfac
){
    wbarray<wbcomplex> A; A.set(A0, wbarray<double>());
    DZGEMM(A,B,C,k,aflag,bflag,wbcomplex(afac,0),cfac);
};

template<> inline
void DZGEMM(
    const wbarray<double> &A,
    const wbarray<double> &B, wbarray<wbcomplex> &C,
    const unsigned k,  
    char aflag, char bflag, double afac, wbcomplex cfac
){
    wbarray<double> C0(C.SIZE); 
    DZGEMM(A,B,C0,k,aflag,bflag,afac);

    C.Plus(C0,double(1),0,cfac);
};

template<class T>  
void Wb::VMatProd(
    const wbarray<T> &A0, const wbarray<T> &B0, wbarray<T> &C,
    char aflag, char bflag, const T afac, const T cfac,
    const char cforce
){
    char isvec=0; wbarray<T> A,B;

    if (A0.SIZE.len==1) { 
       isvec++; A0.adjustVec(FL,aflag,A,1);
    }
    else A.init2ref(A0); 

    if (B0.SIZE.len==1) { 
       isvec++; B0.adjustVec(FL,bflag,B,2);
    }
    else B.init2ref(B0); 

    Wb::MatProd(A,B,C,aflag,bflag,afac,cfac,cforce);

    if (isvec) { C.skipSingletons(); }
};

template<class T>
void Wb::MatVProd(
    const wbarray<T> &A, const wbvector<T> &v, wbarray<T> &x,
    char aflag, char bflag, const T afac, const T cfac,
    const char cforce
){
    wbarray<T> B;
    opFlags<T> Aflag(aflag), Bflag(bflag);

    if (Bflag.conj()) { B.init(v.len,1,v.data).Conj(); } 
    else {
       if (Bflag.trans()) wblog(FL,
          "WRN %s() ignoring bflag=%c",FCT,bflag);
       B.init2ref(v.len,1,v.data);
    }

    Wb::MatProd(A,B,x,aflag,'N',afac,cfac,cforce);
    x.skipSingletons();
};

template<class T>
T Wb::VMatVprod(
   const wbarray<T> &v1, const wbarray<T> &M, const wbarray<T> &v2
){
   unsigned i,dim1,dim2;
   const T *p1=v1.data, *p2=v2.data, *dd=M.data;
   T x=0;

   M.getMatSize(FL,dim1,dim2);

   if (v1.numel()!=dim1 || v2.numel()!=dim2) wblog(FL,
      "ERR %s() severe size mismatch: %s x %s x %s",FCT,
       v1.sizeStr().data, M.sizeStr().data, v2.sizeStr().data);

   for (i=0; i<dim2; i++, dd+=dim1) { 
      x+=Wb::overlap(p1,dd,dim1)*p2[i];
   }

   return x;
};

template<class T>
T Wb::VMatVprod(
   const wbvector<T> &v1, const wbarray<T> &M, const wbvector<T> &v2
){
   unsigned i,dim1,dim2;
   const T *p1=v1.data, *p2=v2.data, *dd=M.data;
   T x=0;

   M.getMatSize(FL,dim1,dim2);

   if (v1.len!=dim1 || v2.len!=dim2) wblog(FL,
      "ERR %s() severe size mismatch: <%d| %s |%d>",
       FCT, v1.len, M.sizeStr().data, v2.len);

   for (i=0; i<dim2; i++, dd+=dim1) { 
      x+=Wb::overlap(p1,dd,dim1)*p2[i];
   }

   return x;
};

template<class TA, class TB, class TC>  
wbarray<TC>& Wb::MatProd(
   const wbarray<TA> &A, const wbarray<TB> &B, wbarray<TC> &C,
   char aflag, char bflag, TA afac, TC cfac,
   char cforce
){
   if ((void*)&A==(void*)&C || (void*)&B==(void*)&C) {
      wbarray<TC> XC; if (cfac!=TC(0)) XC.init(C);
      Wb::MatProd(A,B,XC,aflag,bflag,afac,cfac,cforce);
      XC.save2(C); return C;
   }

   unsigned a1,a2,b1,b2;
   char iflag=0; 

#ifdef WB_CLOCK
   Wb::Clock clk("MatProd",1); 
#endif

   if (A.rank()!=2 || B.rank()!=2) wblog(FL,
      "ERR %s() rank-2 tensors required (%s; %s)",FCT,SSTR(A),SSTR(B));

   a1=A.SIZE[0]; a2=A.SIZE[1]; if (aflag!='N') SWAP(a1,a2);
   b1=B.SIZE[0]; b2=B.SIZE[1]; if (bflag!='N') SWAP(b1,b2);

   if (a2!=b1) wblog(FL,
      "ERR %s() size mismatch %dx%d * %dx%d",FCT,a1,a2,b1,b2);

   if (cfac!=TC(0)) {
       if (C.data==NULL) { iflag=1;
          if (cforce) wblog(FL,          
          "WRN C = A*B + c*[] with c=%s",toStr(double(cfac)).data);
       }
       else if (!C.isMatrix() || C.SIZE[0]!=a1 || C.SIZE[1]!=b2) {
          wblog(FL,"ERR %s() dimension mismatch: C=(%s) =? (%d,%d).",
          FCT, SSTR(C), a1, b2); iflag=1;
       }
   } else iflag=1;

   if (iflag) {
      C.init(a1,b2);
   }

   DZGEMM(A,B,C,a2,aflag,bflag,afac,cfac);

   return C;
};

template<class TA, class TB, class TC>  
void Wb::DMatProd(
   const wbarray<TA> &A0, const wbarray<TB> &B0, wbarray<TC> &C,
   char aflag, char bflag, const TA afac, const TC cfac,
   const char cforce
){
   if (A0.SIZE.len==2 && B0.SIZE.len==2) {
      Wb::MatProd(A0,B0,C, aflag,bflag, afac,cfac, cforce);
      return;
   }

#ifdef WB_CLOCK
   Wb::Clock clk("MatProd",1); 
#endif

   if (A0.SIZE.len<1 || A0.SIZE.len>2 || B0.SIZE.len<1 || B0.SIZE.len>2)
   wblog(FL,"ERR %s() rank <= 2 required (got %d,%d)", FCT,
   A0.SIZE.len,B0.SIZE.len);

   unsigned i,j,k,sa,a1,a2,sb,b1,b2;
   wbarray<TA> A;
   wbarray<TB> B;

   if (A0.SIZE.len==1 && B0.SIZE.len==2) {
      A0.adjustDMat(FL,aflag, afac, A);
      B0.adjustMMat(FL,bflag, 1.  , B);

      sa=A.SIZE[0]; b1=B.SIZE[0]; b2=B.SIZE[1]; sb=b1*b2;
      if (sa!=b1) wblog(FL,
         "ERR %s() sever size mismatch (%d * %dx%d)",
          FCT,sa,b1,b2);

      C.adjustCMat (FL,cfac,b1,b2,cforce);

      if (C.isEmpty()) { C=B;
         for (k=j=0; j<b2; j++)
         for (  i=0; i<b1; i++, k++) C.data[k]*=A.data[i];
      }
      else {
         for (k=j=0; j<b2; j++)
         for (  i=0; i<b1; i++, k++) C.data[k]+=A.data[i]*B.data[k];
      }
   }
   else
   if (A0.SIZE.len==2 && B0.SIZE.len==1) {

      TB bfac; Wb::safeConvert(FL,afac,bfac);

      A0.adjustMMat(FL,aflag, 1.  , A);
      B0.adjustDMat(FL,bflag, bfac, B);

      a1=A.SIZE[0]; a2=A.SIZE[1]; sa=a1*a2; sb=B.SIZE[0];
      if (a2!=sb) wblog(FL,
         "ERR %s() sever size mismatch (%dx%d * %d)",
          FCT,a1,a2,sb);

      C.adjustCMat (FL,cfac,a1,a2,cforce);
      if (C.isEmpty()) { C=A;
         for (k=j=0; j<a2; j++)
         for (  i=0; i<a1; i++, k++) C.data[k]*=B.data[j];
      }
      else {
         for (k=j=0; j<a2; j++)
         for (  i=0; i<a1; i++, k++) C.data[k]+=A.data[k]*B.data[j];
      }
   }
   else wblog(FL,"ERR !?");
};

template<>
wbarray<double>& wbInverse(const char *F, int L, wbarray<double> &M) {

   if (M.rank()!=2 || M.SIZE[1]>M.SIZE[0]) wblog(F_L,
      "ERR %s() invalid rank-2 tensor (%s)",FCT,M.sizeStr().data);

   pINT e=0, l=-1, m=M.SIZE[0], n=M.SIZE[1], N=m*n;
   wbvector<pINT> ipiv(MIN(m,n));
   wbvector<double> work(1);

   if (n<m) wblog(FL,
      "ERR DGETRI() requires dim1>=dim2 (%dx%d)",m,n);

   dgetrf(m,n,M.data,m,ipiv.data,e);
   if (e) wblog(FL,"ERR DGETRF() returned e=%d !?",FCT,e);

   dgetri(n,M.data,m,ipiv.data,work.data,l,e);
   if (e) wblog(FL,"ERR DGETRI() returned e=%d !?",FCT,e);

   l=MIN(pINT(work[0]),N); work.init(l);
   dgetri(n,M.data,m,ipiv.data,work.data,l,e);
   if (e) {
      if (e<0)
           wblog(FL,"ERR DGETRI() got invalid argument #%d",-e);
      else wblog(FL,"ERR DGETRI() got singular matrix (i=%d)",e);
   }

   return M;
};

template<>
wbarray<double>& wbInverse(
   const char *F, int L,
   const wbarray<double>&M, wbarray<double> &Minv
){ Minv=M; return wbInverse(F,L,Minv); };

inline void wbEigenS(
   const wbarray<double> &M, wbarray<double> &V, wbvector<double> &E
){
   unsigned n, r=M.SIZE.len, r2=r/2;
   pINT ni=0, q=0; 
   double nd=0;

   if (r==0) { V.init(); E.init(); return; }

   if (r%2) wblog(FL,"ERR %s() even-rank required (%d)",FCT,r);
   if (!M.isHConj(-1E-12)) { 
      MXPut(FL,"ans").add(M,"M").add(V,"V").add(E,"E").add(r,"r");
      wblog(FL,"ERR %s() got non-symmetric %s matrix",
      FCT, M.sizeStr().data);
   }

   n=Wb::prodRange(M.SIZE.data,r2);

   V=M; E.init(n); if (n==0) return;

   if (r2>1) {
      V.SIZE[r2]=n; V.SIZE.len=r2+1;
   }

   dsyevd('V','L',n,V.data,n,E.data,&nd,-1,&ni,-1,q); 
   if (nd<1 || ni<1 || q) wblog(FL,
      "ERR DSYEVD returned lwork=%g/%d (e=%d) !?",nd,ni,q);

   wbvector<double> wd((pINT)nd);
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

   if (q>0) { pINT q_=q; V.Reset(M.data);
      dsyevd('V','U',n,V.data,n,E.data,wd.data,wd.len,wi.data, wi.len,q);

      if (!q) wblog(FL,
        "WRN DSYEVD() ok only for uplo=L->U (i=%d/%d) !?",q_,n);
      else wblog(FL,
        "ERR DSYEVD() failed to converge (info=%d/%d) !?\n"
        "(D=%d @ l=%ld/%ld)",q_,q,n,E.len,wd.len,wi.len);
   }

   if (q) {
      MXPut(FL,0,"tmpfile").add(M,"M").add(q,"q");
      wblog(FL,"ERR DSYEVD() returned "
      "info=%d (D=%d; l=%ld/%ld) !?",q,n,wd.len,wi.len);
   }
};

inline void wbEigenS(
   const wbarray<wbcomplex> &M, wbarray<wbcomplex> &V, wbvector<double> &E
){
   unsigned n, r=M.SIZE.len, r2=r/2;
   pINT ni=0, q=0; 
   wbcomplex nz=0; double nr=0;

   if (r==0) { V.init(); E.init(); return; }

   if (r%2) wblog(FL,"ERR %s() even-rank required (%d)",FCT,r);
   if (!M.isHConj(-1E-12)) { 
      MXPut(FL,"ans").add(M,"M").add(V,"V").add(E,"E").add(r,"r");
      wblog(FL,"ERR %s() got non-hermitian %s matrix",FCT, SSTR(M));
   }

   n=Wb::prodRange(M.SIZE.data,r2);

   V=M; E.init(n);
   if (n<=1) {
      if (n) {
         double x=V.SkipTiny_imag(1E-12); 
         if (x<0) wblog(FL,"ERR %s() got complex data @ %.3g",FCT,x);
         E.data[0]=V.data[0]; V.data[0]=1;
      }
      return;
   }

   if (r2>1) {
      V.SIZE[r2]=n; V.SIZE.len=r2+1;
   }

   zheevd('V','U',n,V.data,n,E.data,&nz,-1,&nr,-1,&ni,-1,q); 
   if (nz.r<1 || nr<1 || ni<1 || q) {
      MXPut(FL,0,"tmpfile").add(M,"M").add(q,"q");
      wblog(FL,"ERR ZHEEVD returned lwork=%g/%g/%d (e=%d)",nz.r,nr,ni,q);
   }

   wbvector<wbcomplex> wz(nz.r);
   wbvector<double> wr(nr);
   wbvector<pINT> wi(ni);

   zheevd(
      'V',      
      'U',      
       n,       
       V.data,  
       n,       
       E.data,  
       wz.data, wz.len, 
       wr.data, wr.len, 
       wi.data, wi.len, 
       q        
   );

   if (!q) { return; }
   if (q<0) wblog(FL,"ERR %s() ZHEEVD returned e=%d",FCT,q);

   V=M; E.init(n);

   wblog(FL,"WRN ZHEEVD returned i=%d / falling back to ZHEEV",q);
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
      'U',      
       n,       
       V.data,  
       n,       
       E.data,  
       wz.data, ni,    
       wr.data, 
       q        
   );

   if (q) wblog(FL,"ERR %s() ZHEEV returned e=%d",FCT,q);
};

void wbEigen_CS (
   wbarray<wbcomplex> &M,  
   wbarray<wbcomplex> &V,
   wbvector<wbcomplex> &E,
   char wjob,
   char issym, 
   char tnorm,  
   char qflag
){
   unsigned n=0, r=M.SIZE.len, r2=r/2;
   char jobvl='N', jobvr='V', vflag=1; pINT q=0;

   if (wjob) { wjob=toupper(wjob); jobvr=jobvl='N';
      switch (wjob) {
         case 'L': jobvl='V'; break;
         case 'R': jobvr='V'; break;
         case 'N': vflag=0; break;
         default : wblog(FL,"ERR invalid wjob=%c<%d>",wjob,wjob);
      }
   }

   if (r==0) { V.init(); E.init(); return; }

   if (!M.isOpS(n)) wblog(FL,
      "ERR %s() (generalized) square matrix required (%s)",
       FCT,M.sizeStr().data);

   E.init(n); if (vflag) V.init(n,n);

   if (n<=0) {
      if (n) { E[0]=M[0]; V[0]=1; }
      return;
   }

   if (r2>1) {
      V.SIZE[r2]=n; V.SIZE.len=r2+1;
   }

   if (issym) { tnorm=1;
      double x=1E-12;
      if (issym=='T') M.Symmetrize(FL,&x,0); 
      else if (issym=='H') M.Symmetrize(FL,&x,1); 
      else wblog(FL,"ERR invalid issym=%c<%d>",issym,issym);

      if (x>1E-8) {
#ifdef MATLAB_MEX_FILE
         M.put("M_");
#endif
         wblog(FL,"ERR operator not symmetric (%.3g)",x);
      }
      else if (x>1E-12) wblog(FL,
      "WRN operator not quite symmetric (%.3g)",x);
   }

   wbvector<wbcomplex> aux;
   wbvector<double> aux2(2*n);

   unsigned lwork=MAX(2U,n/4)*n; 

   if (n>127) { 
      zgeev (
         jobvl, jobvr, n, M.data, n, E.data, V.data, n, V.data, n,
         aux.data, -1, aux2.data, q);
      lwork=(unsigned)aux[0].r; {
         unsigned b=lwork/n; if (b<2 || b>n)
         wblog(FL,"WRN %s() got lwork %g*%d",FCT,double(lwork)/n,n);
      }
   }

   aux.init(lwork);

#ifdef WB_CLOCK
   Wb::Clock clk("zgeev",1); 
#endif

 { wbarray<wbcomplex> X(M); 
   zgeev (
       jobvl,     
       jobvr,     
       n,         
       X.data,    
       n,         
       E.data,    
       V.data,    
       jobvl=='V' ? n : 1,         
       V.data,    
       jobvr=='V' ? n : 1,         
       aux.data,  
       aux.len,   
       aux2.data, 
       q          
   ); }

   if (q) {
      MXPut(FL,0,"tmpfile").add(M,"M").add(q,"q");
      wblog(FL,"ERR %s() CGEEV returned %d (%d/%d)",FCT,q,aux.len,n);
   }

   if (tnorm) {
      wbEigen_CS_regen(M,V,E,wjob);

      wbarray<wbcomplex> X; double eps=10E-8;
      Wb::MatProd(V,V,X,'T');

      if (!X.isIdentityMatrix(&eps)) {
#ifdef MATLAB_MEX_FILE
         printf("\n"); M.put(FL,"M_"); V.put(FL,"V_"); E.put(FL,"E_",'t');
         wbvector<wbcomplex> Es(E); Es.Sort(); Es.put("Es_");  X.put("X_");
#endif
         if (qflag)
         wblog(FL,"WRN eig() U matrix not quite orthogonal (%.3g)",eps); else
         wblog(FL,"ERR eig() U matrix not quite orthogonal (%.3g)",eps);
      }
   }
}

template<class T>
void wbEigen_CS_regen(
   const wbarray<T> &M0,
   wbarray<T> &V,
   wbvector<T> &E,
   char wjob
){
   if (E.len==0) return;
   if (!V.isRank(2) || E.len!=V.SIZE[1] || V.SIZE[0]!=V.SIZE[1])
   wblog(FL,"ERR invalid usage (%s; %d)",V.sizeStr().data,E.len);

   unsigned i,i0=0, n=E.len, count=0;
   WBINDEX S(2); S[0]=n;
   double a, eps=1E-8, nmin=0, nmax=0;

   wbperm P; E.Sort(P); V.Select0(P,1); 

   V.NormalizeCols(FL,&nmin,&nmax,'T');

   if (nmin>1+1E-14) wblog(FL,
      "ERR zggeev() orthogonal matrix with t-norm = %.3g",nmin); else
   if (nmin<1E-4) wblog(FL,
      "ERR zggeev() orthogonal matrix with t-norm = %.3g",nmin); else
   if (nmin<0.01) wblog(FL, 
      "WRN zggeev() orthogonal matrix with t-norm = %.3g",nmin);

   for (i=1; i<n; i++) {
      a=ABSDIFF_F(E[i],E[i-1]);

      if (i+1<n) { if (a<eps) continue; }
      else { if (a<eps) i++; }

      if (i>i0+1) { T eref; double escale, x=eps;

wblog(FL,"TST fixing degenerate block: %d:%d (%d)  \r\\",i0+1,i,n);

         wbarray<T> U,X,M2,u2; wbvector<T> e2;

         S[1]=i-i0; U.init2ref(V.ref(i0),S);
         Wb::MatProd(M0,U,X); Wb::MatProd(U,X,M2,'T');

         M2.Symmetrize(FL,&x,0); 
         if (x>eps) wblog(FL,
         "WRN matrix not quite symmetric (x=%.4g) !?",x);
         M2.balanceOp(FL,eref,escale);

         wbEigen_CS (M2,u2,e2,wjob,0,0); 

         e2*=escale; e2+=eref;
         Wb::cpyRange(E.data+i0,e2.data,e2.len);

         Wb::MatProd(U,u2,X);
         memcpy(U.data,X.data,sizeof(T)*X.numel());

         U.OrthoNormalizeCols(FL,'T',T(1E-14),2);
         count++;
      }
      i0=i;
   }

};

template<class T>
void wbEigen_CS_regen_trial(
   const wbarray<T> &M0,
   wbarray<T> &V,
   wbvector<T> &E,
   char wjob
){
   if (E.len==0) return;
   if (!V.isRank(2) || E.len!=V.SIZE[1] || V.SIZE[0]!=V.SIZE[1])
   wblog(FL,"ERR invalid usage (%s; %d)",V.sizeStr().data,E.len);

   unsigned i,j,k,i0=0, n=E.len; double a, eps=1E-8;

   wbperm P; E.Sort(P); V.Select0(P,1); 
   wbcomplex *vr,*v,z,z2;
   unsigned count=0;

         wbarray<T> U,X,M2,u2; wbvector<T> e2;
         double eps0, eps2;

         Wb::MatProd(M0,V,X); Wb::MatProd(V,X,M2,'T');
         eps0=eps; M2.isDiagMatrix(&eps0);
         wbarray<T> V0(V);

   for (i=1; i<n; i++) {
      a=ABS(E[i]-E[i-1]);
      if (a<eps) {

         for (j=i0; j<i; j++) { vr=V.ref(j); v=V.ref(i);
           z=z2=0; for (k=0; k<n; k++) { z+=vr[k]*v[k]; z2+=(vr[k]*vr[k]); }
           z=z/z2; for (k=0; k<n; k++) v[k]-=z*vr[k];

           z=z2=0; for (k=0; k<n; k++) { z+=vr[k]*v[k]; z2+=(vr[k]*vr[k]); }
           z=z/z2; for (k=0; k<n; k++) v[k]-=z*vr[k];
         }

wblog(FL,"TST fixing degenerate block: %d-%d (%.3g)",i0,j,std::sqrt(z2.abs())); count++;

      }
      else i0=i;
   }

   if (count) {
        Wb::MatProd(M0,V,X); Wb::MatProd(V,X,M2,'T');
        eps2=eps; M2.isDiagMatrix(&eps2);
wblog(FL,"TST fixing degenerate block: %.6g => %.6g",eps0,eps2);
V0.put("V0_"); V.put("V_"); M0.put("M_"); E.put("E_");
wblog(FL,"ERR");
   }
};

template<class T> inline
void GESVD_old(  
   wbarray<T> &A,
   wbarray<T> &U,
   wbvector<double> &S,
   wbarray<T> &Vd
){
   wblog(FL,"ERR %s() not defined for datatype %s",
   FCT,getName(typeid(T)).c_str());
};

template<> inline
void GESVD_old(
   wbarray<double> &A, 
   wbarray<double> &U, 
   wbvector<double> &S,
   wbarray<double> &Vd
){
   unsigned dim1=A.SIZE[0], dim2=A.SIZE[1], k=S.len;
   pINT i;

   pINT lwork=MAX3(1U, 3*MIN(dim1,dim2)+MAX(dim1,dim2), 5*MIN(dim1,dim2));

   if (MIN(dim1,dim2)>=8) {
      int e=0; double w; 
      dgesvd(
        'S','S', dim1, dim2, A.data, dim1, S.data,
         U.data, dim1, Vd.data, k,
         &w, -1, i); 
      i=(pINT)w;

      if (i<lwork) e=1; else
      if (i>dim1*dim2 && i>(1<<20)) {
         unsigned m=max(dim1,dim2); e=(i<(m*m) ? 2 : 3); 
      }
      if (e) { wblog(FL,
         "WRN dgesvd() lwork=%d/%d (%.4g; %dx%d; e=%d)",
         i,lwork,double(i)/double(dim1*dim2),dim1,dim2,e);
      }
      lwork=i;
   }

   wbvector<double> W(lwork);

   dgesvd(
     'S','S', dim1, dim2, A.data, dim1, S.data,
      U.data, dim1,  
      Vd.data, k,    
      W.data, W.len, i
   );

   if (i) { if (i<0) { i=-i; wblog(FL,
      "ERR %s() dgesvd: invalid argument %d\nhint: %s",FCT,i,
      i==5 || i==8? "check input data for NAN or INF values":"??");
   }
   else { wblog(FL,
      "ERR %s() dgesvd returned e=%d\n"
      "hint: check input data for NAN or INF values",FCT,i);
   }}
};

template<> inline
void GESVD_old(
   wbarray<wbcomplex> &A, 
   wbarray<wbcomplex> &U, 
   wbvector<double> &S,
   wbarray<wbcomplex> &Vd
){
   unsigned dim1=A.SIZE[0], dim2=A.SIZE[1], k=S.len;
   wbvector<wbcomplex> W, RW(5*MIN(dim1,dim2)); 
   pINT i, lwork=MAX(1U, 2*MIN(dim1,dim2)+MAX(dim1,dim2));

   if (MIN(dim1,dim2)>=8) { wbcomplex w; 
      zgesvd(
        'S','S', dim1, dim2, A.data, dim1, S.data,
         U.data, dim1, Vd.data, k,
         &w, -1, RW.data, i  
      );

      i=unsigned(w.r);
      if (i<lwork || (i>10*dim1*dim2 && lwork>1E6)) {
         unsigned m=MIN(dim1,dim2); wblog(FL,"WRN zgesvd(%dx%d) "
        "lwork = %d -> %d (~%d*%d)",dim1,dim2,lwork,i,i/m,m);
      }
      lwork=i;
   }

   W.init(lwork);

   zgesvd(
     'S','S', dim1, dim2, A.data, dim1, S.data,
      U.data, dim1,
      Vd.data, k,
      W.data, W.len, RW.data, i
   );

   if (i) {
       if (i<0) { i=-i;
         wblog(FL,"ERR %s() zgesvd: invalid argument %d\nhint: %s",FCT,i,
         i==5 || i==8 ? "check input data for NAN or INF values":"?");
      }
      else wblog(FL,"ERR %s() zgesvd returned i=%d",FCT,i);
   }
};

template<class T> inline
int GESVD(  
   wbarray<T> &A,
   wbarray<T> &U, wbvector<double> &S, wbarray<T> &Vd,
   char lflag=0 
);

template<> inline
int GESVD(
   wbarray<double> &A, 
   wbarray<double> &U, 
   wbvector<double> &S,
   wbarray<double> &Vd, char lflag
){
   unsigned dim1=A.SIZE[0], dim2=A.SIZE[1], k=S.len;
   pINT q, ni=8*MIN(dim1,dim2); 
   double nd;

   wbvector<double> wd;
   wbvector<pINT> wi(ni);

#ifdef WB_CLOCK
   Wb::Clock clk("dgesvd",1); 
#endif

   dgesdd('S',dim1,dim2,A.data,dim1,S.data,U.data,dim1,Vd.data,k,
      &nd,-1,wi.data,q);
   if (nd<1 || q) wblog(FL,
      "ERR DGESDD() returned lwork=%g (e=%d) !?",nd,q);
   wd.init(nd);

   dgesdd(
     'S', dim1, dim2, A.data, dim1, S.data,
      U.data, dim1,  
      Vd.data, k,    
      wd.data, wd.len, wi.data,
      q
   );

   if (q) {
      if (q<0) { q=-q; wblog(FL,
         "ERR DGESDD() invalid argument %d\n%s",q, q==5 || q==8 ?
         "hint: check input data for NAN or INF values" : ""); }
      else if (lflag) {
         wblog(FL,"WRN DGESDD() did not converge (info=%d) !?",q);
      }
      else {
         MXPut(FL,0,"tmpfile").add(A,"A").add(q,"q");
         wblog(FL,"ERR DGESDD() did not converge (info=%d) !?",q);
      }
   }
   return q;
};

template<> inline
int GESVD(
   wbarray<wbcomplex> &A, 
   wbarray<wbcomplex> &U, 
   wbvector<double> &S,
   wbarray<wbcomplex> &Vd, char lflag
){
   unsigned dim1=A.SIZE[0], dim2=A.SIZE[1], k=S.len;
   pINT q, n1=MIN(dim1,dim2), n2=MAX(dim1,dim2);

   pINT nd = n1*(3*n1+6 >= 2*n2 ? 5*n1+7 : 2*(n1+n2)+1);

#ifdef WB_CLOCK
   Wb::Clock clk("zgesvd",1); 
#endif

   wbcomplex nz;
   wbvector<wbcomplex> wz;
   wbvector<double> wd(nd);
   wbvector<pINT> wi(8*n1); 

   if (!dim1 || !dim2) wblog(FL,
      "ERR %s() got %dx%d matrix!?",FCT,dim1,dim2);

   zgesdd('S',dim1,dim2,A.data,dim1,S.data,U.data,dim1,Vd.data,k,
      &nz, -1, wd.data, wi.data, q);
   if (nz.r<1 || q) wblog(FL,
      "ERR ZGESDD() returned lwork=%g (e=%d) !?",nz.r,q);
   wz.init(nz.r);

   zgesdd(
     'S',dim1,dim2,A.data,dim1,S.data, U.data,dim1,Vd.data,k,
      wz.data, wz.len, wd.data, wi.data,
      q
   );

   if (q) {
      if (q<0) { q=-q; wblog(FL,
         "ERR ZGESDD() invalid argument %d\n%s",q, q==5 || q==8 ?
         "hint: check input data for NAN or INF values" : ""); }
      else if (lflag) {
         wblog(FL,"WRN ZGESDD() failed to converge (info=%d) !?",q);
      }
      else {
         wblog(FL,"ERR ZGESDD() failed to converge (info=%d) !?",q);
      }
   }
   return q;
};

template<class T>
void wbSVD(
   const wbarray<T> &A,
   wbarray<T> &U, wbvector<double> &S,
   wbarray<T> &Vc,   
   const WBINDEX &I, 
   char bare 
){
   unsigned dim1,dim2,k, r=A.SIZE.len, l=I.len;
   int itry=0, ntry=2, q=0; 
   wbarray<T> X,Vd;
   wbperm P;

   static unsigned ncall=0; 
   static int nout=0; ++ncall;

   for (; itry<ntry; ++itry) {
      if (l==0) { 
         if (r!=2) wblog(FL,
            "ERR %s() rank-2 required, use index otherwise (%d)",FCT,r);
         X.init2ref(A);
      }
      else {
         if (l>=r) wblog(FL, 
         "ERR %s() invalid dimensions (%s; %d/%d)",FCT,I.toStr().data,l,r);
         A.toMatrixRef(X,I,1,P);
      }

      if (X.isRef()) X.Instantiate();

      dim1=X.SIZE[0]; dim2=X.SIZE[1]; k=MIN(dim1,dim2);

      U.init(dim1,k); S.init(k); Vd.init(k,dim2);
      if (!dim1 || !dim2) return;

      if (itry==0) {
         q=GESVD(X,U,S,Vd,'l');
         if (!q) break; 
      }
      else {

         char vn[8]; snprintf(vn,8,"Isvd%d",(++nout)%10);
         MXPut Ix(FL,vn,"tmpfile");
         Ix.add(ncall,"ncall").add(q,"q").add(X,"A");

         wblog(FL,"WRN GESVD returned i=%d / falling back to ZGESVD",q);
         GESVD_old(X,U,S,Vd); 
         Ix.add(U,"U").add(S,"S").add(Vd,"Vd");
      }
   }

   #ifdef SVD_BUG_SAFEGUARD 
   { double a=SQRT(A.norm2()), s=SQRT(S.norm2());
     double e=fabs(a-s)/(a*sqrt(double(dim1)*double(dim2)));

     if (e>1E-10) { char s_[128];
        if (l==0) { X.init2ref(A); } else { A.toMatrixRef(X,I,1,P); }
        if (X.isRef()) X.Instantiate();

        X.Transpose();
        dim2=X.SIZE[0]; dim1=X.SIZE[1]; k=MIN(dim1,dim2);

        Vd.init(dim2,k); S.init(k); U.init(k,dim1);

        GESVD(X,Vd,S,U); Vd.Transpose(); U.Transpose();

        s=SQRT(S.norm2());
        a=fabs(a-s)/(a*sqrt(double(dim1)*double(dim2)));
        if (ISCOMPLX_(T))
             strcpy(s_,"ZGESVD");
        else strcpy(s_,"DGESVD");
        snprintf(s_+10,118,"%dx%d @%9.3g -> %.2g",dim1,dim2,e,a);

        if (a<1E-10)
             wblog(FL,"WRN fixed %s bug: %s (ok)",   s_,s_+10);
        else wblog(FL,"ERR failed to fix %s bug\n%s",s_,s_+10);
     }
   }
   #endif

   if (l) {
      WBINDEX q,Q; A.SIZE.select(P,Q);
      q.init(l+1,  Q.data    ); q[l]=S.len; U.Reshape(q);
      q.init(r-l+1,Q.data+l-1); q[0]=S.len; Vd.Reshape(q);
   }

   if (bare) { Vd.save2(Vc); } 
   else {
      P.initFirstTo(Vd.SIZE.len-1,Vd.SIZE.len);
      Vd.permute(Vc,P); 
   }
};

#endif

