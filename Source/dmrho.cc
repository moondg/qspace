/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace NRG routines
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

#ifndef __WB_DMRHO_CC__
#define __WB_DMRHO_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

template <class TQ, class TD>
double initRHO(
   wbMatrix<double> &rhoNorm,
   const double T,
   NRGData<TQ,TD> &H,
   wbvector<double> E0,
   const unsigned dloc,
   const unsigned NRho=0,
   char vflag=1
);

template <class TQ, class TD>
void getEData(
   const QSpace<TQ,TD> &HK,
   wbvector<double> &E, wbvector<unsigned> &D,
   wbvector<unsigned> *dz=NULL 
);

template <class TQ, class TD>
TD getErange(const QSpace<TQ,TD> &HK, double *emin=NULL, double *emax=NULL);

template <class TQ, class TD>
unsigned getGSDegeneracy(const QSpace<TQ,TD> &H);

template <class TQ, class TD>
bool getBoltzman(
   const QSpace<TQ,TD> &H,
   wbvector<double> &R, wbvector<unsigned> &D, unsigned iter,
   wbvector<unsigned> *DZ=NULL, char useE0=1
);

template <class TD>
bool getBoltzman_base(
   const wbvector<TD>& E, const wbvector<unsigned> &dz,
   wbvector<TD>& R, unsigned iter, double T=0., char useE0=1
);

template <class TQ, class TD>
double initRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &HK,
   const int iter,
   double efac=1.
);

template <class TQ, class TD>
void nrgUpdateRHO(
   const wbMatrix<double> &R, 
   NRGData<TQ,TD> &A,
   NRGData<TQ,TD> &H,
   NRGData<TQ,TD> &RHO,
   unsigned NRho,
   unsigned dloc, 
   Wb::SigHandler *sig=NULL
);

template <class TQ, class TD>
void updateRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &AK
);

template <class TQ, class TD, class TX>
void updateOp_SR(QSpace<TQ,TX> &X, const QSpace<TQ,TD> &AK);

template <class TQ, class TD, class TZ>
void initG_SR(
   QSpace<TQ,TZ> &G2, const QSpace<TQ,TD> &H2,
   wbcomplex z, unsigned iter
);

template <class TQ, class TD>
void checkRho (const char *F, int L,
    const QSpace<TQ,TD> &Rho,
    const double traceRho=1.
){
    unsigned i,n=Rho.QIDX.dim1;
    double rsum=0, rmin=0;

    wbvector< wbvector<double> > E(n);
    wbvector<double> rr;
    wbarray<TD> U;

    char cgflag=Rho.gotCGS(FL);

    if (cgflag>0 && Rho.qtype.permitsOM()>1) {
       cgflag=2;
    }

    if (!Rho.isConsistent_r(2) || !Rho.isBlockDiagMatrix() || !Rho.isHConj())
    wblog(F,L,"ERR %s() inconsistency in Rho QSpace (%d;%d;%d)",
    FCT, Rho.isBlockDiagMatrix(), Rho.isHConj());

    for (i=0; i<n; i++) {
       wbEigenS(*(Rho.DATA[i]),U,E[i]);
       if (cgflag>0) {
          unsigned d=Rho.cgsDimScalar(i);
          if (cgflag<=1) { 
             unsigned d0=Rho.qtype.QDim(Rho.QIDX.rec(i)); if (d0!=d)
             wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,d,d0);
          }
          E[i]*=d;
       }
    }
    rr.Cat(E);

    rmin=rr.min();
    rsum=rr.sum();

    if (rmin<-EPS) { wblog(F,L,
       "ERR eig(rho) is in (%.4g .. %.4g; %.4g)",rmin,rr.max(),rsum);
    }
    else if (fabs(rsum-traceRho)>EPS) {
       char e=(fabs(rsum-traceRho)>1E-8 && rsum>1E-8 && traceRho>1E-8);
       if (e) {
          printf("\n\n");
          MXPut(FL).add(Rho,"Rho").add(rr,"rr").add(traceRho,"t");
       }
       sprintf(str,"%s tr(rho) inconsistency !?\n"
         "%.4g / %.4g @ %.3g (%.2g; %.2g; %.2g)", e ? "ERR":"WRN",
          rsum,traceRho,rsum-traceRho,DEPS,DEPS2,EPS); printf("\n");
       wblog(F,L,str);
    }
};

template <class TQ, class TD>
double initRHO(
   wbMatrix<double> &rhoNorm, 
   const double T,
   NRGData<TQ,TD> &H,
   wbvector<double> E0, 
   const unsigned dloc,
   const unsigned NRho, 
   char vflag
){
   unsigned i=0,i1,i2, iter, isneg=0, N=NRG_N; char i0flag=0;
   wbvector<double> R;
   wbvector<unsigned> D;
   double dbl, dfac=1; 
   double DFAC=DBL_MAX*1E-2;

   if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)");

   i1 = (dloc && NRho ? NRho-1 : 0);
   i2 = (NRho ? NRho-1 : N-1);

   if (int(i1)>int(i2) || int(i1)<0) wblog(FL,
      "ERR rho_i = [%d..%d] !?",i1,i2);

   rhoNorm.init(N,2); 

   dbl=E0.min(i);
   if (i==0) { i0flag=1;
      E0[0]+=E0[1]; dbl=E0.min(i);
      E0[0]-=E0[1];
   }

   if (i<i2) {
      wblog(FL,"WRN Emin set by site %d/%d !?",i+1,N);
      MXPut(FL).add(E0,"E0");
   }

   getBoltzman_base(E0,D,E0,-1,T); 

   for (iter=i1; iter<=i2; ++iter) {
      H.init(FL,"D",iter);
      if (getBoltzman(H.D,R,D,iter)) ++isneg;
      rhoNorm(iter,1)=R.sum(); 

      if (dloc && iter<i2) continue;

      H.init(FL,"K",iter);
      if (getBoltzman(H.K,R,D,iter)) ++isneg;
      rhoNorm(iter,0)=R.sum(); 
   }

   if (isneg) wblog(FL,
      "WRN %s() negative (eigen)energies (%d) !?",FCT,isneg);
   if (i0flag && (rhoNorm(0,0)>1 || rhoNorm(0,1)>1)) wblog(FL,
      "WRN %s() density matrix from first site may be unstable (%.3g,%.3g)",
       FCT, rhoNorm(0,0), rhoNorm(0,1));

   if (!H.K.isEmpty())
        i=getGSDegeneracy(H.K); 
   else i=getGSDegeneracy(H.D);

   if (vflag) {
      if (i==1)
           wblog(FL," *  non-degenerate ground state");
      else wblog(FL,"<i> ground state degeneracy k=%d (%d)",iter-1,i);
   }

   if (!dloc) {
      if (!NRho) wblog(FL,"ERR NRho=%d !?", NRho);
      wblog(FL,"<i> generate weights for plain NRG (T=%.3g)", T);

      for (iter=0; iter<NRho; iter++) {
         dbl=rhoNorm.recSum(iter);
         if (dbl<1) wblog(FL,"WRN sum_rhoNorm(%d,:)=%.3g",iter+1,dbl);
         dbl=1/dbl;
         rhoNorm(iter,0)*=dbl;
         rhoNorm(iter,1)*=dbl;
      }
      return 0;
   }

   if (!NRho)
   for (dfac=dloc, iter=i2-1; (int)iter>=0; iter--, dfac*=dloc) {
      rhoNorm(iter,1)*=dfac;

      if (dfac>DFAC) wblog(FL,"WRN dfac reaches DBL_MAX (%.3g)", dfac);
   }

   dbl=rhoNorm.sum();
   if (dbl<=0.) wblog(FL,"ERR %s() normRH0=%g ??? DIV/0",FCT,dbl);
   rhoNorm *= (1/dbl);

   if (vflag) {
      double rmax;
      rhoNorm.recSum(R); rmax=R.max(i);
      wblog(FL," *  maxRho=%.4g at iter=%d/%d (T=%.4g)",rmax,i+1,N,T);
   }

   return (-T*log(dbl));
};

template <class TQ, class TD>
void getEData(
   const QSpace<TQ,TD> &HK,   
   wbvector<double> &E,       
   wbvector<unsigned> &D,     
   wbvector<unsigned> *DZ     
){
   unsigned i,l,m=0,d, e=0, k,K=HK.QIDX.dim1, *dz=NULL;
   char cgflag=HK.gotCGS(FL);

   if (HK.DATA.len!=K) wblog(FL,
      "ERR %s() severe QSpace inconsistency (%d,%d)",
       FCT,HK.DATA.len,HK.QIDX.dim1);
   if (!HK.isBlockDiagMatrix('d')) wblog(FL, 
      "ERR %s() Hamiltonian must be (block)diagonal!",FCT);

   D.init(K);
   for (l=k=0; k<K; k++) D[k]=HK.DATA[k]->SIZE.max();

   E.init(D.sum());
   if (DZ && cgflag>0) { DZ->init(E.len); dz=DZ->data; }

   if (cgflag>0 && HK.qtype.permitsOM()>1) {
      cgflag=2;
   }

   for (l=k=0; k<K; ++k) {
       const wbarray<TD> &Hk = *HK.DATA[k];
       const wbvector<unsigned> &S = Hk.SIZE;
       const TD* const E0=Hk.data;

       d=D[k]; if (dz) {
          m=HK.cgsDimScalar(k);
          if (cgflag<=1) { 
             unsigned m0=HK.qtype.QDim(HK.QIDX.rec(k)); if (m0!=m)
             wblog(FL,"ERR %s() cgsDim inconsistency! (%d/%d)",FCT,m,m0);
          }
       }

       if (S[0]!=S[1]) { if (S[0]!=1 && S[1]!=1) e++; else {
          Wb::cpyRange(E.data+l, E0, d);
          if (dz)
               { for (i=0; i<d; ++i,++l) dz[l]=m; }
          else { l+=d; }
       }}
       else { if (S.len!=2) e++; else {
          if (d>1 && (fabs(E0[1])>1E-12 || fabs(E0[d])>1E-12)) wblog(FL,
          "WRN HK[%d] not E diagonal (%g,%g)", k+1,E0[1], E0[d]);

          for (i=0; i<d; ++i,++l) {
             E[l]=E0[i*d+i]; if (dz) dz[l]=m;
          }
       }}
       if (e) wblog(FL,
          "ERR HK[%d] must be rank-2 tensor! (%s)",
           k+1, HK.DATA[k]->sizeStr().data
       );
   }
};

template <class TQ, class TD>
unsigned getGSDegeneracy(const QSpace<TQ,TD> &H) {

   wbvector<double> E; wbvector<unsigned> D,dz;
   unsigned i=0, n=0;
   double e0, *ee;

   getEData(H,E,D,&dz); e0=E.min(); ee=E.data;
   for (; i<E.len; ++i) if (fabs(ee[i]-e0)<EPS) n++;

   return n;
};

template <class TQ, class TD>
TD getErange(const QSpace<TQ,TD> &H, double *emin, double *emax){

   if (!H.isBlockDiagMatrix('d')) {
      QSpace<TQ,TD> AK;
      QSpace<TQ,double> HK; wbMatrix<double> EE;
      H.EigenSymmetric(AK,HK,EE);
      return getErange(HK,emin,emax);
   }

   wbvector<double> E; wbvector<unsigned> D;
   TD e1,e2;

   getEData(H,E,D);
      e1=E.min(); if (emin) (*emin)=e1;
      e2=E.max(); if (emax) (*emax)=e2;
   return (e2-e1);
};

template <class TQ, class TD>
bool getBoltzman(
   const QSpace<TQ,TD> &H,
   wbvector<double> &R, wbvector<unsigned> &D, unsigned iter,
   wbvector<unsigned> *DZ, char useE0
){
   wbvector<unsigned> dz; 
   wbvector<double> E;
   bool e=0;

   getEData(H,E,D,&dz); e=(E<0);

   if (getBoltzman_base(E,dz,R,iter,0.,useE0)) { e=1; }
   if (DZ) dz.save2(*DZ);

   return e;
};

template <class TD>
bool getBoltzman_base(
   const wbvector<TD>& E,
   const wbvector<unsigned>& dz, 
   wbvector<TD>& R, 
   unsigned iter,   
   double T0,       
   char useE0
){
   static double T=0.; 
   static wbvector<double> E0; 

   unsigned i,e=0;

   if (int(iter)<0) { T=T0; E0=R; return 0; } 

   if (iter>=E0.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,iter+1,E0.len);
   if (dz.len && dz.len!=E.len) wblog(FL,
      "ERR %s() size mismatch (%d/%d)",FCT,dz.len,E.len);

   if (E.isEmpty()) { R.init(); return 0; }

   double beta,r, Escale=gEScale(iter);
   double Eref=(useE0 ? E0[iter]/Escale : 0.);

   if (E<-Eref) { ++e; wblog(FL, 
      "WRN Eref should be lower bound (%.4g < %.4g)  @ %.4g !?",
       E.min(), Eref, Escale);
   }

   if (&R!=&E) R.init(E.len);

   if (T>0.) {
      beta=Escale/T;
      for (i=0; i<E.len; ++i) {
         r=exp(-beta*(E[i]+Eref)); if (dz.len) r*=dz[i];
         R[i]=r;
      }
   }
   else {
      for (i=0; i<E.len; i++)
      R[i]=( ABS(E[i]+Eref)<EPS ? (dz.len ? dz[i] : 1) : 0. );

      i=unsigned(R.sum());
      if (i) wblog(FL," *  T=0 ground state degeneracy %d", i);
      else   wblog(FL," *  T=0 got unique ground state");
   }

   return e;
};

template <class TQ, class TD>
double initRho(
   QSpace<TQ,TD> &Rho,
   const QSpace<TQ,TD> &H,
   const int iter,
   double rhofac
){
   unsigned i,k,l,d,n, K=H.QIDX.dim1, has0=0;
   wbvector<unsigned> D,dz;
   wbvector<double> R;
   double Z, *r;

   getBoltzman(H,R,D,iter,&dz,0);

   Z=R.sum();

   if (Z==0.) {
      MXPut(FL,"Ix").add(H,"H").add(R,"R").add(iter,"iter").add(dz,"dz");
      wblog(FL,"ERR rhoNorm[%d]=%g -> DIV/0 !?", iter+1,Z);
   }
   if (dz.len && R.len!=dz.len) wblog(FL,
      "ERR %s() length mismatch %d/%d",FCT,R.len,dz.len);

   R *= (rhofac/Z); r=R.data;

   if (dz.len) {
      for (i=0; i<R.len; ++i) { if (dz.data[i]>1) r[i]/=dz.data[i]; }
   }

   Rho=H; 

   for (l=k=0; k<K; k++, l+=d) {
       wbarray<TD> &Rk = (*Rho.DATA[k]);

       d=D[k]; Rk.init(d,d);

       for (n=i=0; i<d; i++) {
          Rk.data[i*d+i] = r[l+i];
          if (!n && r[l+i]>DEPS) n++;
       }
       if (!n) has0++;
   }

   if (has0) Rho.SkipZeroData(DEPS2,'b'); 
   return Z;
};

template <class TQ, class TD>
void updateRho(QSpace<TQ,TD> &Rho, const QSpace<TQ,TD> &AK) {

   QSpace<TQ,TD> Xk;

   if (!Rho.isBlockDiagMatrix()) wblog(FL,
      "WRN %s() Rho not blockdiag !?",FCT);

   AK.contract(2,Rho,1,Xk); 
   Xk.contract("3,2",AK,"2,3;*",Rho); 
};

template <class TQ, class TD>
void nrgUpdateRHO(
   const wbMatrix<double> &R, 
   NRGData<TQ,TD> &A,
   NRGData<TQ,TD> &H,
   NRGData<TQ,TD> &RHO,
   unsigned NRho,
   unsigned dloc, 
   Wb::SigHandler *sig
){
   unsigned iter, ilast=0, rzero, N=NRG_N; char ie;
   const size_t n=16; size_t l=0; char tag[n];
   wbvector<double> w; double w2, wtot=0;

   QSpace<TQ,TD> Rho;

   RHO.K.clearQSpace();

   if (R.dim1!=N || R.dim2!=2)
      wblog(FL, "ERR %s() size mismatch of rhoNorm (%dx%d; %d %d)",
      FCT, R.dim1, R.dim2, N, NRho);

   for (iter=N-1; (int)iter>=0; iter--) {
      if (sig) { sig->check911(); }

      R.getRec(iter,w); rzero=(w<DEPS);
      w2=w.sum();

      l=snprintf(tag,n,
         "%s %2d/%d", !NRho ? "FDM":"DMK", iter,N);
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%s; %d/%d)",FCT,tag,l,n);

      RHO.updatePara(FL,"K",iter);

      if (rzero && RHO.K.isEmpty()) {
         wblog(FL,"%s skip ...%16s\r\\", tag,"");
         continue;
      }

      A.init(FL,"K",iter);  
      A.init(FL,"D",iter);  
      H.init(FL,"D",iter);  

      A.check_dloc(FLF,dloc,iter); 

      if (!(ie=RHO.K.isEmpty())) {
      updateRho(RHO.K, A.K); } 

      if (!NRho) {
         if (H.D.isEmpty()) wblog(FL,
            "%s no discarded space (%.3g) %16s\r\\",tag,wtot,"");
         else if (w2!=0) wblog(FL,
            "%s dRho=%.3g (%.3g) %16s\r\\",tag,w2,wtot,"");
         else wblog(FL,
            "%s dRho=%.3g %16s\r\\",tag,wtot,"");
      }
      else wblog(FL,"%-30s\r\\",tag);

      if (!rzero) {
         if (w[1]) {
            initRho(Rho, H.D, iter, w[1]);
            updateRho(Rho, A.D); 
            RHO.K+=Rho; if (!ilast) { ilast=iter; }
         }

         if (w[0] && (!ilast || ilast==iter)) { H.init(FL,"K",iter);
            initRho(Rho, H.K, iter, w[0]);
            updateRho(Rho, A.K);
            RHO.K+=Rho; if (!ilast) { ilast=iter; }
         }

         RHO.K.SkipZeroData(DEPS2,'b');
      }

      wtot+=w2;

      if (!RHO.K.isEmpty()) { try { 
         checkRho(FL, RHO.K, Wb::addRange(R.ref(iter,-1),N-iter,R.dim2));
      }
      catch (...) { 
         wblog(FL,"ERR %s() iter=%d/%d%s",FCT,iter+1,N,
         iter+1==N ? "\nhint: got kept space at last iteration?":"");
      }}
      else {
         if (ie && rzero) wblog(FL,
            "%s empty RHO (will not contribute) %16s\r\\",tag,"");
         else wblog(FL,"WRN %2d/%d got empty RHO !? %16s",iter,N,"");
      }
      doflush(); 
   }
};

template <class TQ, class TD, class TZ>
void initG_SR(
   QSpace<TQ,TZ> &G2, const QSpace<TQ,TD> &H2,
   wbcomplex z, unsigned iter
){
   unsigned k=0, K=H2.QIDX.dim1;
   size_t i,dk; TZ *gk;

   double Escale=gEScale(iter);

   G2.initT(H2); 

   for (; k<K; ++k) {
      gk=G2.DATA[k]->data;
      dk=G2.DATA[k]->dim1();
      for (i=0; i<dk; ++i) { gk[i]=1./(z-Escale*gk[i]); } 
      G2.DATA[k]->ExpandDiagonal();
   }
};

template <class TQ, class TD, class TX>
void updateOp_SR(QSpace<TQ,TX> &X, const QSpace<TQ,TD> &AK) {

   QSpace<TQ,TX> Xk;

   if (!X.isBlockDiagMatrix()) wblog(FL,
      "ERR %s() input X not blockdiag !?",FCT);
   if (X.itags.len!=2) wblog(FL,
      "ERR %s() input X not a scalar op (%s) !?",FCT,IT2STR(X));

   AK.contract(2,X,1,Xk); 
   Xk.contract("32",AK,"23*",X); 
};

#endif

