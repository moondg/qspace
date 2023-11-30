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

/* COMMENTS / CHANGE LOG ============================================= *
   Wb,Dec20,10 : switched wbArray (deprecated) => wbarray (col-major)
   Wb,Feb18,19 : changed (C1,C2) -> (B,C)
 * =================================================================== */

char USAGE[]=""; // outsourced to fgrNRG.m // Wb,Feb14,19

char USAGE_2[] = 
"  Usage: fgrNRG('<NRG0>','<NRG2>','<setup>.mat','B','C')          \n\
";

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "fgrNRG"
#endif

#define  EPS 1E-12
#define DEPS DBL_EPSILON
#define NLEN 128  

#define LOAD_CGC_QSPACE
#include "wblib.h"

#include "spectral.hh" 
#include "nrgdata.hh"
#include "dmrho.cc"

template <class TQ, class TD, class TB, class TS>
void fgrNRGIter(
   const NRGData<TQ,TB> &B, const NRGData<TQ,TB> &C,
   const QSpace<TQ,TD> &Rho, char wrho, 
   const NRGData<TQ,TD> &H0, const NRGData<TQ,TD> &H2, const double dEk,
   TDSpectral<TS> &fgrData, int iter
);

template <class TQ, class TD, class TS>
void updateSpec(
   const char *F, int L,
   TDSpectral<TS> &fgrData,
   unsigned iter, unsigned io,
   QSpace<TQ,TS> &RX, 
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2,
   const double dEk
);

template <class TQ, class TD>
void updateSpec(
   const char *F, int L,
   TDSpectral<double> &fgrData, 
   unsigned iter, unsigned io,
   QSpace<TQ,wbcomplex> &RX,    
   const QSpace<TQ,TD> &H1, const QSpace<TQ,TD> &H2, const double dEk
){
   wblog(FL,"ERR %s() invalid type setting\n"
   "wbcomplex Xop into real spectral data !?",FCT);
};

template <class TQ, class TD>
void updateSpec(
   const char *F, int L,
   TDSpectral<wbcomplex> &fgrData, 
   unsigned iter, unsigned io,
   QSpace<TQ,double> &RX,    
   const QSpace<TQ,TD> &H1, const QSpace<TQ,TD> &H2, const double dEk
){
   wblog(FL,"ERR %s() invalid type setting\n"
   "real Xop into complex spectral data !?",FCT);
};

#ifdef MATLAB_MEX_FILE
static void myCleanUp(void) {
}
#endif

   double Lambda=1;

template <class TS>
void MEX_FUNCTION(
   int nargout, mxArray** &argout,
   int nargin, const mxArray** &argin, TDSpectral<TS> &fgrData);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   int i=0;
   for (; i<nargin; ++i) {
      if (mxIsChar(argin[i]) && Mx::IsEqual(argin[i],"--rixs")) break;
   }

   if (i>=nargin) {
      TDSpectral<double> fgrData;
      MEX_FUNCTION(nargout, argout, nargin, argin,fgrData);
   }
   else {
      TDSpectral<wbcomplex> fgrData;
      MEX_FUNCTION(nargout, argout, nargin, argin,fgrData);
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in fgrNRG"); }
   aclu.Check();
}; 

template <class TS>
void MEX_FUNCTION(
   int nargout, mxArray** &argout,
   int nargin, const mxArray** &argin, TDSpectral<TS> &fgrData
){
   unsigned i,k,m,nloc=0,iter, e=0, N=0, QDIM, NRho=0, dloc=0, nlog=256;
   char disp=0, store=1, partial=0, RAW=0, raw=0, c=0; 
   char rixs=0, calcRho=0, locRho=0, calcOps=0;
   double gE0, gE2, TN,dbl, emin, dEk=0, Eoffset=0, emax=10., T=-1;
   double alpha=0, sigma=0;

   char vtag[4], itag[16];
   mxArray *a, *ar=NULL;
   const mxArray *ac=NULL, *bc=NULL;

   wbMatrix<TS> a0,aa;
   wbMatrix<wbcomplex> at;

   wbstring vstr(64), istr(64), NRG0, NRG2;

   NRGData<gTQ,gTD> A0("A"), A2("A"), H0("H"), H2("H");

   NRGData<gTQ,gTD> B("B"), C("C"); 
   NRGData<gTQ,gTD> RHO("RHO");
   NRGData<gTQ,wbcomplex> B02("BH"), VX0("VX");

   wbvector< wbvector< QSpace<gTQ,gTD> > > CC,DD;
   wbvector< QSpace<gTQ,gTD> > C_,D_,Rho;
   QSpace<gTQ,gTD> Z0;

   wbMatrix<gTD> rhoNorm, dd;
   wbvector<double> om,E0,E2,w2,t0;
   wbvector<char> zz;
   OPTS opts;

   Wb::SigHandler SIG(FL); 

   MX_CHECK_HELPER_NARGS(0,-1,-1);
   NRG_N=0; NRG_ITER=0; STRICT_ITER0=0; gES.init();

   memset(str,0,STRLEN); 

#ifdef WB_CLOCK
   Wb::Clock clk(myname,0); 
#endif
#ifdef MATLAB_MEX_FILE
   mexAtExit(myCleanUp);
#endif

#ifndef MAIN

   if (nargin<5 || nargout>3) wblog(FL,
      "ERR invalid number of I/O arguments (%d/%d)",nargin,nargout);

   for (k=0; k<2; ++k) if (mxIsChar(argin[k])) {
      i=mxGetString(argin[k],str,NLEN);
      if (i || !str[0]) wblog(FL,
         "ERR invalid NRG file specs (arg #%d: %d,%d)",k+1,i,strlen(str));
      if (k==0) NRG0=str; else NRG2=str;
   }

   if (!Mx::IsDblScalar(FL,argin[k])) wblog(FL,
      "invalid Lambda (arg#%d)",k+1);
   if (mxGetNumber(argin[k++],Lambda)) wbdie(FL,str);

   ac=argin[k++]; c=(ac && mxIsCell(ac) ? 1 : 0);
   if (c)
        { Wb::initLocalOp(FL,ac,CC,k,  '!'); }
   else { Wb::initLocalOp(FL,ac,C_,k,  '!');  }

   if (nargin>int(k+2) &&
      ((mxIsStruct(argin[k+1]) || mxIsCell(argin[k+1])))
   ){
      bc=argin[k++]; if (c ^ (bc && mxIsCell(bc))) wblog(FL,
         "ERR %s() operator cell-structure mismatch",PROG);

      if (c)
           { Wb::initLocalOp(FL,bc,DD,k+1,'!'); }
      else { Wb::initLocalOp(FL,bc,D_,k+1,'!'); }
   }

   opts.init(argin+k, nargin-k);

#else
   if (nargout || nargin<4 || nargin>5) wblog(FL,
      "ERR invalid number of I/0 arguments (%d/%d)",nargin,nargout);

   { char *sk=str;
     for (k=0; k<(unsigned)nargin; ++k, sk+=NLEN) {
        if (!mxIsChar(argin[k])) wblog(FL,
           "ERR got non-string input in usage #2 "
           "(arg#%d/%d)%N%N%s%N",k+1,nargin,USAGE_2);
        i=mxGetString(argin[k],sk,NLEN-1);
        if (i || !sk[0]) wblog(FL,
           "ERR unexpected string (arg#%d/%d: %d/%d)",
           k+1,nargin,i,strlen(sk)
        );
     }
   }

   NRG0=str;
   NRG2=str+NLEN;

   opts.init(str+2*NLEN); k=3;

   opts.getOpt("Lambda", Lambda,'!');

   B.name=str+(k++)*NLEN; opts.getOpt(B.name.data,a,'!');

   if (c)
        { Wb::initLocalOp(FL,a,CC,k,'!'); }
   else { Wb::initLocalOp(FL,a,C_,k,'!'); }

   if (nargin>4) { 
      mxArray *b=NULL;
      C.name=str+(k++)*NLEN; opts.getOpt(C.name.data,b,'!');

      c=(a && mxIsCell(a)) if (c ^ (b && mxIsCell(b))) wblog(FL,
         "ERR %s() operator cell-structure mismatch (usage #2)",PROG);

      if (c)
           { Wb::initLocalOp(FL,b,DD,k+2,'!'); }
      else { Wb::initLocalOp(FL,b,D_,k+2,'!'); }
   }

#endif

   A0 .setupIO(FL, NRG0.data, argin[0], "0");
   H0 .setupIO(FL, NRG0.data, argin[0]);
   RHO.setupIO(FL, NRG0.data, argin[0]);

   A2 .setupIO(FL, NRG2.data, argin[1], "2");
   H2 .setupIO(FL, NRG2.data, argin[1]);
   B  .setupIO(FL, NRG2.data, argin[1]); 
   C  .setupIO(FL, NRG2.data, argin[1]);

   rixs=opts.getOpt("--rixs"); 
   if (rixs ^ ISCOMPLX_(TS)) { 
      wblog(FL,"ERR %s() got data type %s in %s mode !?",
      myname,TSTR(TS),rixs? "RIXS":"FGR");
   }

   if (rixs) {
      VX0.setupIO(FL, NRG0.data, argin[0]); 
   }

   if (Lambda<=1.) wblog(FL,"ERR invalid Lambda=%g",Lambda);

   H0.checkVec(FL,"K",2,"NRG data (HK0)"); 
   H2.checkVec(FL,"K",2,"NRG data (HK)" ); 

   H0.checkVec(FL,"D",2,"NRG data (HD0)");
   H2.checkVec(FL,"D",2,"NRG data (HD)" );
   A0.checkVec(FL,"D",3,"NRG data (AD0)");
   A2.checkVec(FL,"D",3,"NRG data (AD)" );

   if (NRG_N<=0 ) wblog(FL,"ERR got empty NRG data !?"); else
   if (NRG_N>999) wblog(FL,"ERR got NRG chain length N=%d !?",NRG_N);
   N=NRG_N; i=0;

   H0.getMxInfo(FL,"E0",E0);
   H2.getMxInfo(FL,"E0",E2);

   if (E0.len!=N || E2.len!=N) wblog(FL,
      "ERR length mismatch in E0 (N=%d/%d/%d)",E0.len,E2.len,N);

   H0.getMxInfo(FL,"EScale",gES);
   H2.getMxInfo(FL,"EScale",w2); 
   if (N!=gES.len || N!=w2.len) wblog(FL,
      "ERR mismatch in chain length (EScale: %d/%d/%d)",gES.len,w2.len,N);
   if (gES!=w2) { dbl=gES.normDiff(w2); wblog(FL,
      "ERR %s() got different EScale data\nfrom %s vs %s @ e=%.3g !?",
      myname,H0.wnrg(),H2.wnrg(),dbl); };
   wblog(FL," *  using EScale from %s_info.mat",H0.wnrg());

   dloc=getDLoc(A2); 

   A2.init(FL,"K",0); QDIM=A2.K.QDIM;

   TN=gEScale(N-1); 
   emin=TN*1E-3;

#ifndef MATLAB_MEX_FILE
   if (Wb::GetEnv(FL,"FDM_T",T)>=0)
#endif
   opts.getOpt("T",T);

   partial=opts.getOpt("partial");
   raw=opts.getOpt(FL,"raw");
   RAW=opts.getOpt(FL,"RAW");
   opts.getOpt(FL,"alpha",alpha); sigma=(!rixs && alpha>0 ? 0 : 0.3);
   opts.getOpt(FL,"sigma",sigma);

   if (opts.getOpt("tt",a)) {
      t0.init(FL,a);
      if (rixs && t0.len) wblog(FL,"ERR %s() no tt option for RIXS",FCT);
   }

   opts.getOpt(FL,"Eoffset",Eoffset);

   if (!rixs) { if (!RAW && !raw && alpha<=0 && sigma>0) {
      alpha=1; 
   }}

   if (opts.getOpt("nostore")) { store=0; }  
   locRho =opts.getOpt("locRho" ); 
   calcRho=opts.getOpt("calcRho"); 
   calcOps=opts.getOpt("calcOps"); 

   opts.getOpt("Z0",a); if (a) {
      if (mxIsQSpace(a)<=0) wblog(FL,"ERR invalid QSpace Z0");
      Z0.init(FL,a);
   }

   opts.getOpt("disp", disp);
   opts.getOpt("nlog", nlog);
   opts.getOpt("emin", emin);
   opts.getOpt("emax", emax);
   opts.getOpt("NRho", NRho);

   opts.getOpt("zflags",a); if (a) zz.init(FL,a,"zflags");

   if (NRho>N) {
      wblog(FL,"WRN NRho out of bounds (%d -> %d)", NRho,N);
      NRho=N;
   }

   if (rixs) {
      gRX.init(); 
      if (!opts.getOpt("Einc",gRX.Einc)) wblog(FL,
         "ERR %s() missing RIXS incoming photon enegy Einc",myname);
      opts.getOpt("mode",gRX.wrxs);
      opts.getOpt("afac",gRX.afac);
      gRX.sigma=sigma;
      gRX.dloc=dloc;

      gRX.check(FL,1);
   }
   opts.checkAnyLeft(); 

   if (NRho) { if (T<0) T=0; }
   else {
      dbl=sqrt(2/(log(Lambda)*log(double(dloc)))); 
      m=(unsigned)ceil(5*dbl+3); 
      if (m>=N) m=N-1;

      dbl=gEScale(N-m-1);
      if (T<0) { T=dbl;
         wblog(FL," *  using T=%.4g (= T[N-%d] = %.4g TN)",T,m,T/TN); }
      else if (T<dbl) {
         wblog(FL,"WRN T=%.3g < T_{N-%d}=%.3g",T,m,dbl);
      }
   }

   if (rixs) {
      if (NRho) { strcpy(vtag,"rxK"); vstr.printf(FL,
      "RIXS (using single-shell @ NRho=%d!)",NRho); }
      else { strcpy(vtag,"RXS"); vstr.cpy(FL,
      "RIXS energy loss (using Kramers-Heisenberg)"); }
   }
   else {
      if (NRho) { strcpy(vtag,"fgK"); vstr.printf(FL,
      "Fermi golden rule (using single-shell @ NRho=%d!)",NRho); }
      else { strcpy(vtag,"FGR"); vstr.cpy(FL,
      "absorption/emission (using Fermi golden rule)"); }
   }

   str[0]=0;
   if (disp) sprintf(str," disp=%d", disp);
   if (!store ) strcat(str," nostore");
   if (calcOps) strcat(str," calcOps");
   if (locRho ) strcat(str," locRho" );
   if (calcRho) strcat(str," calcRho");

   printf("\n"
      "  %-18s %s computing %s\n"
      "  %-18s N=%d, QDIM=%d, d=%d%s%s\n"
      "  %-18s %.4g (= %.4g TN)\n"
      "  %-18s %.2g .. %.2g (%d/dec)\n",
   "program:",     myname, vstr.data,
   "parameters:",  N, QDIM, dloc, str[0] ? "," : "", str[0] ? str : "",
   "temperature:", T, T/TN, emin,
   "omega range:", emax, nlog);

   if (t0.len) printf(
   "  %-18s %.2g .. %.2g (len=%d)\n","time range",t0[0], t0.last(), t0.len);

   if (rixs) {
      printf("\n  RIXS incoming photon energy: %s\n",STR(gRX));
   }; printf("\n");

   if (!CC.isEmpty()) {
      nloc=0; 

      if (zz.len || !Z0.isEmpty()) wblog(FL, 
         "WRN %Nsince C is specified as cell array, z-ops must be\n"
         "already applied in C (zflags and Z0 will be ignored!)%N"
      );

      if (!DD.len) {
         for (k=i=0; i<CC.len; ++i) {
            nloc=MAX(nloc,unsigned(CC[i].len));
            if (CC[i].len) { if (k<i) { CC[k]=CC[i]; }; ++k; }
         }
         if (k<i) { wblog(FL,
            "WRN ignoring empty operators (%d/%d)",i-k,i);
            CC.len=k;
         }
      }
      else {
         if (DD.len!=CC.len) wblog(FL,
            "ERR length mismatch of C with B (%d/%d)",DD.len,CC.len);

         for (k=i=0; i<CC.len; ++i) {
            nloc=MAX(nloc, unsigned(MAX(CC[i].len, DD[i].len) ));
            if (CC[i].len && DD[i].len) {
               if (k<i) { CC[k]=CC[i]; DD[k]=DD[i]; }
               ++k;
            }
         }
         if (k<i) { wblog(FL,
            "WRN ignoring empty operator pairs (%d/%d)",i-k,i);
            CC.len=k; DD.len=k;
         }
      }

      B.initOp(FL,CC, NULL,"1st");
      C.initOp(FL,DD, &B,  "2nd");
   }
   else { 

      if (D_.len && D_.len!=C_.len) wblog(FL,
         "ERR length mismatch of C with B (%d/%d)",D_.len,C_.len);

      B.initOp(FL,C_); B.initSym("KK");
      C.initOp(FL,D_); C.initSym("KK");

      if (calcOps) { B.forceCalc(); C.forceCalc(); }

      if (!zz.data) { zz.init(C_.len); zz.set(1); }
      else {
         if (zz.len!=C_.len) {
            if (zz.len==1) { zz.Resize(C_.len); zz.set(zz[0]); }
            else wblog(FL,
               "ERR length mismatch of zflags with C (%d/%d)",
               zz.len,C_.len
            );
         }
      }

      if (zz.anyUnequal(0)) { e=0;
         wblog(FL," *  applying Z0 [zflags = %s]", zz.toStr().data);

         if (Z0.isEmpty()) wblog(FL,
            "ERR Z0 required (having zflags)");
         if ((B.KK.len && B.KK[0].rank()!=Z0.rank()) ||
             (C.KK.len && C.KK[0].rank()!=Z0.rank()) ) wblog(FL,
            "ERR rank mismatch of Z0 with ops (%d/%d) !?",
             Z0.rank(),B.KK[0].rank());

         if (B.calc>0) {
            for (i=0; i<B.KK.len; ++i) {
               if (zz[i]) { e+=applyZ0(Z0,B.KK[i]); }
            }
         }

         if (C.calc>0) {
            for (i=0; i<C.KK.len; i++) {
               if (zz[i]) { e+=applyZ0(Z0,C.KK[i]); }
            }
         }

         if (e) wblog(FL,
            "ERR QIDX of Z0 does not match QIDX of %s%s%s\n"
            "(hint: unset Z0 to empty if not relevant.", B.name.data,
            C.KK.len ? " or " : "", C.KK.len ? C.name.data : ""
         );
      }
   }

   if (B.KK.isEmpty()) wblog(FL,"ERR got empty B.KK !?");

   if (calcOps) {
      B.forceCalc();  
      C.forceCalc();
   }
   if (!store) {
      B.store=0; if (C.nrgIdx.len) {
      C.store=0; }
   }

   if (C.KK.isEmpty())
        { B.name="C";     B.info(FL,"C" ); }
   else { B.info(FL,"B"); C.info(FL,"C"); }

   if (1 || B.store) {
      if (CC.isEmpty())
           B.updateInfop(FL,"i",C_.toMx());
      else B.updateInfop(FL,"i",B.CI.toMx());
   }
   if (1 || (C.store && C.CI.len)) {
      if (DD.isEmpty())
           C.updateInfop(FL,"i",D_.toMx());
      else C.updateInfop(FL,"i",C.CI.toMx());
   }
   if (!B.calc && !C.calc) wblog(FL,"<i> using stored operator sets");

   fgrData.init(FL,vtag, t0, 
      N, B.KK.len, emin, emax, nlog
   );

   fgrData.checkAvgOm();

   if (RAW) fgrData.raw=1;

   gE0=nrgGetE0(E0,gES, E0); 
   gE2=nrgGetE0(E2,gES, E2); 

   dbl=E0[NRho>0 ? NRho-1 : N-1];
   E0-=dbl; initRHO(rhoNorm, T, H0, E0, dloc, NRho);
   E2-=dbl;

   dbl = E2.last() - E0.last();

   if (dbl || Eoffset) {
      wblog(FL,"<i> Eoffset = %.4g %+.4g = %.8g",Eoffset,dbl,Eoffset+dbl);
      Eoffset += dbl;
      if (rixs) gRX.Eoffset=Eoffset;
   }

   dd.init(RHO.getMxInfo("rhoNorm"));

   if (calcRho || dd!=rhoNorm || RHO.checkVec(FL,"",2) ||
      (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T)
   ){
      str[0]=0;
      if (locRho && !RHO.MX) {
         sprintf(str," internally (locRho=%d)",locRho);
         RHO.MX=mxCreateStructMatrix(1,NRG_N,0,NULL);
         RHO.locMX=1; locRho=99;
      }
      wblog(FL,"<i> update density matrices%s%N",str);

      nrgUpdateRHO(rhoNorm, A0,H0, RHO, NRho, dloc, &SIG);
      printf("\n");
   }
   else wblog(FL,"<i> using stored DM set%N");

   if (locRho!=99) locRho=0;
   Rho.initDef(2);

   if (rixs) {
      if (C.KK.len && !B.KK.isEqual(C.KK)) wblog(FL,
         "WRN %s() RIXS is expected to have C==B !?",FCT);

      if (!B.store) {
         if (B.MX) wblog(FL,"ERR %s() already got MX=%p !?",FCT,B.MX);
         B  .setupMX(FL,NRG_N,2,1); 
      }; B02.setupMX(FL,NRG_N,2,1); 

      for (iter=0; iter<N; ++iter) { NRG_ITER=iter; SIG.check911();
         A0.initX(FL); A2.initX(FL);

         wblog(FL,"--> %s %2d/%d matrix elements B^XX' \r\\",vtag,iter+1,N);

         if (iter+2<=nloc && !A2.D.isEmpty()) wblog(FL,
            "WRN %Nlocal operator extends to trancated sites (%d/%d)",
             iter+1, nloc-1
         );

         B.updateOp(FL,A0,A2,4,iter); 
         B.saveOp(FL,iter,4);
      }
      printf("\n");

      VX0.setupMX(FL,NRG_N,2,1); 

      B02.KK.init(B.KK.len);  

      for (--iter; iter<N; --iter) { NRG_ITER=iter; SIG.check911();
         H0.initX(FL); H2.initX(FL);
         if (iter+1<N) { A0.initX(FL); A2.initX(FL); B.initXX(FL); }

         wblog(FL,"<-- %s %2d/%d (BHeisenberg)^K \r\\",vtag,iter+1,N);

         gRX.setE0(iter, E0[iter], E2[iter]); 

         B02.initBG(FL,B,H0,H2,3); 
         B02.saveOp(FL,iter,4); 

         if (iter) {
         B02.backPropagateOp(FL,A0,A2); } 

      }
      printf("\n");

      if (1 && locRho) { 
         if (RHO.MX) { RHO.PUT(FL,"RHO","caller"); }
         if (B02.MX) { B02.PUT(FL,"B02","caller"); }
         if (B  .MX) { B  .PUT(FL,"B",  "caller"); }
      }
   }

   for (iter=0; iter<N; ++iter) { NRG_ITER=iter; SIG.check911();

      snprintf(itag,16,"%s %2d/%d",vtag,iter+1,N);

      A0.initX(FL); H0.initX(FL); 
      rhoNorm.getRec(iter,w2); 

      if (!rixs) {
         A2.initX(FL); H2.initX(FL);
         if (iter+2<=nloc && !A2.D.isEmpty()) wblog(FL,
            "WRN %Nlocal operator extends to trancated sites (%d/%d)",
            iter+1, nloc-1);

         B.updateOp(FL,A0,A2,4,iter); B.saveOp(FL,iter,4); 
         C.updateOp(FL,A0,A2,4,iter); C.saveOp(FL,iter,4);

         dEk = Eoffset + E0[iter] - E2[iter]; 
      }
      else {
         B02.initXX(FL); B.initXX(FL);
         VX0.updateOpProd(FL,iter,B02,B); 
         if (VX0.MX) { VX0.saveOp(FL,iter,4); } 
      }

      if (w2.anyGT(DEPS)) { 
         if (w2[0]>DEPS) {  
            if (w2[1]>DEPS) wblog(FL,"WRN %s "
               "kept (%.3g) + truncated (%.3g) ",itag,w2[0],w2[1]);
            else wblog(FL, "WRN %s kept (%.3g) ",itag,w2[0]);
         }
         else if (w2[1]>DEPS) { wblog(FL,
            "%s truncated (%.3g)%10s\r\\",itag,w2[1],"");
         }

         for (i=0; i<2; ++i) { if (w2[i]) {
            const QSpace<gTQ,gTD> &Hi = (i ? H0.D : H0.K);
            char wi=(i? 'D':'K');

            initRho(Rho[i], Hi, iter, w2[i]);

            if (rixs)
                 { fgrNRGIter(VX0,VX0,Rho[i],wi,H0,H0,0., fgrData,iter); }
            else { fgrNRGIter(B,  C,  Rho[i],wi,H0,H2,dEk,fgrData,iter); }
         }}
      }
      else {
         if (disp>2) wblog(FL,"%s (RHO)",itag);
         else wblog(FL,"%s (RHO)%30s\r\\",itag,"");
      }

      RHO.init(FL,"K",iter);
      if (rixs)
           { fgrNRGIter(VX0,VX0,RHO.K,'K',H0,H0,0., fgrData,iter); }
      else { fgrNRGIter(B,  C,  RHO.K,'K',H0,H2,dEk,fgrData,iter); }
   }

   if (disp<=2) printf("\n\n");

   if (getSpecDimOps(B.DD.len ? B.DD :
      (B.KK.len ? B.KK : (C.DD.len ? C.DD : C.KK)), w2))
   fgrData.applyIROPfac(FL,w2,rixs);

   if (RAW) {
      fgrData.addPartialFourier();
      ar=fgrData.At.toMx('r');
   }

   fgrData.getSpecData(om,a0);

   if (!rixs && t0.len && !raw) { 
      if (sigma>0) {
         wblog(FL,"<i> frequency broadened TDM data (sigma=%g)",sigma);

         fgrData.getSmoothSpec_t(sigma,0.); 
         fgrData.Fourier(om,aa,t0,at,alpha);
      }
      else if (alpha>0) {
        if (Lambda<=1) wblog(FL,"ERR alpha requires Lambda (%g)",Lambda);

         fgrData.Fourier(om,at,alpha,Lambda); 
      }
      else {
         wblog(FL,"<i> plain Fourier transformed data (no broadening)");

         Wb::Fourier(om,a0,t0,at);
      }
   }

   if (VX0.MX) { VX0.PUT(FL,"VX0","caller"); } 

   MXPut Iout(0,0);

   RHO.init(FL,"K",0);
   A0.init(FL,"K",0);

   Iout.add(vstr,"info").add(Wb::TimeStamp(),"finished");
   if (rixs)
        { Iout.add(gRX,"Irxs"); } 
   else { Iout.add(Eoffset,"Eoffset"); }

   Iout.add(rhoNorm,"rho").add(T,"T")
       .add(A0.K,"A0").add(RHO.K,"RHO");

   { MXPut Sx(0,0); wbvector<double> rr(2); rr[0]=raw; rr[1]=RAW;
     Sx.add((char*)vtag,"vtag");
     Sx.add(emin,"emin").add(emax,"emax").add(nlog,"nlog")
       .add(E0,"E0").add(E2,"E2")
       .add(gES,"EScale")
       .add(gE0,"gE0").add(gE2,"gE2")  
       .add(fgrData.getAvgOm(),"avgOm")
       .add(rr,"raw");
     if (partial) Sx.add(fgrData,"RAW");
     if (ar)      Sx.addP(ar, "cc_raw");

     mxArray *F=mxCreateCellMatrix(1,2);
     mxSetCell(F,0,H0.NAME.toMx());
     mxSetCell(F,1,B.NAME.toMx());
     Sx.addP(F,"NRG");

     Iout.addP(Sx.toMx(),"Inrg");
   }

   if (!t0.isEmpty() || !at.isEmpty()) {
      Iout.addP(MXPut(0,0)
         .add(alpha,"alpha")
         .add(sigma,"sigma")
         .addP(t0.toMx('t'),"tt")
         .addP(at.toMx('r'),"at")
         .addP(aa.toMx('r'),"aa")
      .toMx(),"Itdm");
   }

   if (nargout<=1) { Iout.addP(a0.toMx('r'),"a0"); }
   if (nargout<=2) { Iout.addP(om.toMx('t'),"om"); }

   if (!B.NAME.isEmpty()) {
      B.updateInfo(FL,"om", om.toMx('t'));
      B.updateInfo(FL,"a0", a0.toMx('r'));
      B.updateInfo(FL,"Ifgr",Iout.S,'k'); 
   }

   if (!RHO.NAME.isEmpty()) {
      RHO.updateInfo(FL,"rhoT",    numtoMx(T));
      RHO.updateInfo(FL,"rhoNorm", rhoNorm.toMx());
   }

   argout[MAX(0,nargout-1)]=Iout.toMx(); 
   if (nargout>=2) { argout[nargout-2]=a0.toMx('r'); }
   if (nargout>=3) { argout[0]=om.toMx('t'); }

#ifdef MATLAB_MEX_FILE
   myCleanUp();
#endif

   Wb::Clock *clio=fdmClocks.get("data:I/O",1);
   wblog(FL,"I/O %s time usage: %s%N", myname,
      clio ?  SEC2STR(clio->gettime()) : "0"); 
   if (clio) { clio->reset(); }

}; 

template <class TQ, class TD, class TB, class TS>
void fgrNRGIter(
   const NRGData<TQ,TB> &B, const NRGData<TQ,TB> &C,
   const QSpace<TQ,TD> &Rho, char w, 
   const NRGData<TQ,TD> &H0, const NRGData<TQ,TD> &H, const double dEk,
   TDSpectral<TS> &fgrData,
   int iter
){
   unsigned i,k,n,m=0;
   QSpace<TQ,TB> R2; 
   char tags[][3] = { "*K","*D" }; 

   if (iter<0) iter=NRG_ITER;

   n=B.KK.len;
   if (n!=B.KD.len || n!=B.DK.len || n!=B.DD.len) wblog(FL,
     "ERR severe operator dim mismatch [%d %d %d %d]", 
      B.KK.len, B.KD.len, B.DK.len, B.DD.len
   );
   if (&C!=&B && C.KK.len)
   if (n!=C.KD.len || n!=C.DK.len || n!=C.DD.len) wblog(FL,
     "ERR severe operator dim mismatch [%d %d %d %d]", 
      C.KK.len, C.KD.len, C.DK.len, C.DD.len
   );

   if (w!='K' && w!='D') wblog(FL,"ERR invalid w=%c<%d> !?",w,w);
   for (i=0; i<2; ++i) tags[i][0]=w;

   if (w=='K') m=1; 

   for (i=m; i<2; ++i) {
      const QSpace<TQ,TD>
         &H1 = H0.getQSpace(tags[i][0]),
         &H2 = H .getQSpace(tags[i][1]);

      for (k=0; k<n; ++k) {
         const QSpace<TQ,TB>
            &C1k = B.getQSpace(tags[i],k),
            &C2k = (C.KK.len ? C.getQSpace(tags[i],k) : C1k);

         Rho.contract(2,C2k,1,R2);

         R2.TimesEl(C1k,'*',2); 

         if (R2.isEmpty()) {
            continue;
         }

         updateSpec(FL,fgrData,iter,k, R2,H1,H2, dEk);
      }
   }
};

template <class TQ, class TD, class TS>
void updateSpec(
   const char *F, int L,
   TDSpectral<TS> &fgrData,
   unsigned iter, unsigned io,
   QSpace<TQ,TS> &RX, 
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2,
   const double dEk
){
   unsigned i,j,k,m,n;

   wbvector<unsigned> I1,I2;
   wbvector<TD> E1,E2;
   wbindex i1h,i2h;

   wbarray<TD> Omega;
   TD *e1,*e2;

   matchH12(FL,RX,H1,H2,iter, i1h,E1,I1, i2h,E2,I2,'r');

   E2 -= dEk; 

   for (k=0; k<i1h.len; ++k) {
      i=i1h[k]; m=I1[i+1]-I1[i];
      j=i2h[k]; n=I2[j+1]-I2[j]; Omega.init(m,n);

      if (!RX.DATA[k]->sameSize(Omega)) { 
         MXPut(FL,"Ix").add(RX,"RX").add(k+1,"k")
           .add(H1,"H1").add(H2,"H2").add(Omega,"Omega");
         RX.DATA[k]->info(FL,"RX",k,3); Omega.info(FL,"Omega");
         wblog(FL,"ERR %s() dimension mismatch (%s <> %s)",
         FCT,SSTR_(RX.DATA[k]),SSTR(Omega));
      }

      e1=E1.data+I1[i]; e2=E2.data+I2[j];

      for (i=0; i<m; ++i)
      for (j=0; j<n; ++j) Omega(i,j)=e2[j]-e1[i];

      fgrData.Add(Omega, *(RX.DATA[k]), iter,io);
   }
};

