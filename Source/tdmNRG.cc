
/* CHANGE LOG ======================================================== *

   NB! C2^DAGGER also important when applying Z0
   from the left for both C1 and C2 (eg. C1 AND C2 being destruction
   or creation operators!  Wb,Jun27,08 #TFLAG

   Wb,Feb25,19: changed T -> D (for discarded state space)
   outsource updateSpec() section in to matchH12()

 * =================================================================== */

char USAGE[]=""; // outsourced to tdmNRG.m // Wb,Jan12,19

char USAGE_2[] = 
"  tdmNRG('nrgdata0', 'nrgdata2', 'setup.mat', 'B', 'C')                \n\
";

#define PROG_TAG "TDM"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "tdmNRG"
#endif

#define  EPS 1E-12
#define DEPS DBL_EPSILON
#define NLEN 128  

#define LOAD_CGC_QSPACE
#include "wblib.h"

#include "spectral.hh" 
#include "nrgdata.hh"
#include "dmrho.cc"

template <class TQ, class TD>
void initOverlap(
   NRGData<TQ,TD> &P,
   const NRGData<TQ,TD> &A0,
   const NRGData<TQ,TD> &A
);

template <class TQ, class TD, class TS>
void tdmNRGIter(
   NRGData<TQ,TD> &C1, 
   NRGData<TQ,TD> &C2,
   NRGData<TQ,TD> &P,
   QSpace<TQ,TD> &Rho, char wrho, 
   NRGData<TQ,TD> &H,
   TDSpectral<TS> &tdmData,
   int iter,
   const char* KKflag=""
);

template <class TQ, class TD, class TS>
void updateSpec(
   const char *F, int L,
   TDSpectral<TS> &tdmData,
   unsigned iter, unsigned io,
   QSpace<TQ,TD> &RX, 
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2
);

#ifdef MATLAB_MEX_FILE
static void myCleanUp(void) {
}
#endif

double Lambda=1;
char fgrFlag=0;

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::gpara.init(); try {

   unsigned i,k,m,nloc=0,iter, e=0, N=0, QDIM, NRho=0, dloc=0, nlog=256;
   char disp=0, store=1, partial=0, RAW=0, lflag, raw=0, gotC2=0;
   char calcRho=0, locRho=0, calcOps=0;
   double TN,dbl, emin, emax=10., wx=0., T=-1;
   double alpha=0, sigma=0;

   char vstr[64], vtag[4], itag[16];
   mxArray *S, *a, *ar=NULL;
   const mxArray *ac, *bc=NULL;

   wbMatrix<double> a0,aa;
   wbMatrix<wbcomplex> at;
   TDSpectral<double> tdmData;

   wbstring istr(128), NRG0, NRG2;

   NRGData<gTQ,gTD> A0("A"), H0("H"), A("A"), H("H"), RHO("RHO");
   NRGData<gTQ,gTD> C1("C1"), C2("C2"), P("P"); 
   wbvector< wbvector< QSpace<gTQ,gTD> > > CC,DD;
   wbvector< QSpace<gTQ,gTD> > C0,D0,Rho;
   QSpace<gTQ,gTD> Z0;

   wbMatrix<gTD> rhoNorm, dd;
   wbvector<double> om,DE,E0,w,t0;
   wbvector<char> zz;
   OPTS opts;

   Wb::SigHandler SIG(FL); 

   MX_CHECK_HELPER_NARGS(1,-1,-1); 

   NRG_N=0; NRG_ITER=0; gES.init();

#ifdef WB_CLOCK
   Wb::Clock clk(myname,0); 
#endif
#ifdef MATLAB_MEX_FILE
   mexAtExit(myCleanUp);
#endif

#ifndef MAIN

   if (nargin<5 || nargout>2 || mxIsChar(argin[2])) wblog(FL,
   "ERR Invalid number of I/O arguments (%d,%d).",nargin,nargout);

   for (k=0; k<2; k++) if (mxIsChar(argin[k])) {
      i=mxGetString(argin[k],str,NLEN);
      if (i || !str[0]) wblog(FL,
         "ERR Error reading file tag from arg #%d (%d,%d)",
          k+1, i, strlen(str));
      if (k==0) NRG0=str; else NRG2=str;
   }

   if (!Mx::IsDblScalar(FL,argin[k])) wblog(FL,
   "Invalid argument %d (Lambda).",k+1);
   if (mxGetNumber(argin[k++], Lambda)) wbdie(FL,str);

   gotC2 = nargin>int(k+2) &&
     (mxIsStruct(argin[k+1]) || mxIsCell(argin[k+1]));

   ac=argin[k++]; if (gotC2) {
   bc=argin[k  ]; }

   if (ac && mxIsCell(ac)) {
      Wb::initLocalOp(FL,ac,CC,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,bc,DD,k+1,'!'); }
   }
   else {
      Wb::initLocalOp(FL,ac,C0,k,  '!'); if (gotC2) {
      Wb::initLocalOp(FL,bc,D0,k+1,'!'); }
   }

   if (gotC2 || (int(k+1)<nargin && mxIsEmptyQSpace(argin[k]))) k++;

   t0.init(FL,argin[k++]);

   opts.init(argin+k, nargin-k);

#else

   if (nargout || nargin<4 || nargin>5) wblog(FL,
   "ERR Invalid number of I/0 arguments (%d,%d).", nargin, nargout);

   for (i=0; i<(unsigned)nargin; i++) {
      if (!mxIsChar(argin[i])) wblog(FL,
      "ERR Invalid arg #%d (not a string)%N%N%s%N", i+1, USAGE_2);
   }

   for (k=0; k<2; k++) {
      i=mxGetString(argin[k],str,NLEN);
      if (i || !str[0]) wblog(FL,
         "ERR Error reading file tag from arg #%d (%d,%d)",
          k+1, i, strlen(str));
      if (k==0) NRG0=str; else NRG2=str;
   }

   if ((i=mxGetString(argin[k++],str,NLEN))) wblog(FL,
   "ERR Error reading string of arg #%d (%d,%d)", k+1,i,strlen(str));

   opts.init(str); 

   opts.getOpt("tt",a,'!'); t0.init(FL,a);  
   opts.getOpt("Lambda", Lambda,'!');

   for (i=0; i<2 && k<nargin; i++) {
      if (mxGetString(argin[k++],str,STRLEN))
      wblog(FL,"ERR Error reading string of arg #%d", k);

      opts.getOpt(str,a,'!');
      if (i==0) C1.init(FL,a,C0,k,'!',str);
      else      C2.init(FL,a,D0,k,'!',str);
   }

#endif

   if ((NRG0.isEmpty() || NRG2.isEmpty()) && nargout<2) wblog(FL,
   "ERR %s requires at least two output arguments", myname);

   A0.setupIO(FL, NRG0.data, argin[0], "0");
   H0.setupIO(FL, NRG0.data, argin[0]);
   RHO.setupIO(FL,NRG0.data, argin[0]);

   A .setupIO(FL, NRG2.data, argin[1], "2");
   H .setupIO(FL, NRG2.data, argin[1]);
   C1.setupIO(FL, NRG2.data, argin[1]); 
   C2.setupIO(FL, NRG2.data, argin[1]);

   if (Lambda<=1.) wblog(FL,"ERR Invalid Lambda=%g.",Lambda);

   H0.checkVec(FL,"K",2,"NRG data (HK0)"); 
   H .checkVec(FL,"K",2,"NRG data (HK)" ); 

   H0.checkVec(FL,"D",2,"NRG data (HT0)");
   H .checkVec(FL,"D",2,"NRG data (HD)");
   A0.checkVec(FL,"D",3,"NRG data (AT0)");
   A .checkVec(FL,"D",3,"NRG data (AD)");

   if (NRG_N==0) wblog(FL, "ERR Input NRG data is empty.");
   N=NRG_N;

   ac=H0.getMxInfo("E0"); if (ac) DE.init(FL,ac);
   else wblog(FL,"ERR failed to retrieve E0 from NRG data");

   if (N!=DE.len) wblog(FL,
      "ERR Mismatch in chain length (N=%d, E0:%d)", N, DE.len);

   ac=RHO.getMxInfo("EScale");
   if (ac==NULL) wblog(FL,"ERR failed to read EScale (info.mat)");
   else {
      gES.init(FL,ac);
      if (N!=gES.len) wblog(FL,
         "ERR mismatch in chain length (EScale: %d/%d)",gES.len,N);
      else wblog(FL," *  using EScale from info.mat");
   }

   dloc=getDLoc(A); 

   A.init(FL,"K",0); QDIM=A.K.QDIM;

   TN=gEScale(N-1); 
   emin=TN/1000.;

#ifndef MATLAB_MEX_FILE
   if (Wb::GetEnv(FL,"FDM_T",T)>=0)
#endif
   opts.getOpt("T", T);

   fgrFlag=opts.getOpt("-fgr");
   partial=opts.getOpt("partial");
   raw=opts.getOpt(FL,"raw");
   RAW=opts.getOpt(FL,"RAW");
   opts.getOpt(FL, "alpha",alpha); sigma=(alpha>0 ? 0 : 0.3);
   opts.getOpt(FL, "sigma",sigma);

   if (!RAW && !raw && alpha<=0 && sigma>0) {
      alpha=1; 
   }

   store =!opts.getOpt("nostore"); 
   locRho =opts.getOpt("locRho" );
   calcRho=opts.getOpt("calcRho");
   calcOps=opts.getOpt("calcOps");

   opts.getOpt("Z0",a); if (a) {
      if (mxIsQSpace(a)<0) wblog(FL,"ERR Invalid Z0 (%s)", str);
      Z0.init(FL,a);
   }

   opts.getOpt("disp", disp);
   opts.getOpt("nlog", nlog);
   opts.getOpt("emin", emin);
   opts.getOpt("emax", emax);
   opts.getOpt("NRho", NRho);

   opts.getOpt("zflags",a); if (a) zz.init(FL,a,"zflags");

   if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)", NRho, N);

   opts.checkAnyLeft(); 

   if (NRho==0) {
      dbl=sqrt(2/(log(Lambda)*log(double(dloc))));
      m=(unsigned)ceil(5*dbl+3); 

      dbl=gEScale(N-m-1);
      if (T<0) {
         T=dbl; wblog(FL," *  taking T=%.4g (=T_{N-%d})", T, m);
      }
      else if (T>0 && T<dbl)
      wblog(FL,"WRN T=%.3g < T_{N-%d}=%.3g", T, m, dbl);
   }
   else if (T<0) T=0;

   if (NRho) {
       sprintf(vstr,"tDM-NRG (single shell NRho=%d)", NRho);
       strcpy(vtag,"tDK");
   }
   else {
       strcpy(vstr,"tDM-NRG (generalized DM)");
       strcpy(vtag,"tDM");
   }

   str[0]=0;
   if (disp) sprintf(str," disp=%d", disp);
   if (!store ) strcat(str, " nostore");
   if (locRho ) strcat(str, " locRho" );
   if (calcOps) strcat(str, " calcOps");
   if (calcRho) strcat(str, " calcRho");

   wblog(FL,"<i> %47R\n%-12s: %s\n"
       "parameters  : N=%d, QDIM=%d, d=%d%s%s\n"
       "temperature : %.4g (= %.4g TN)\n"
       "tdmData     : %.2g .. %.2g (%d/dec)",
   "=", myname, vstr, N, QDIM, dloc, str[0] ? "  ":"", str[0] ? str:"",
   T, T/TN, emin, emax, nlog);

   if (t0.len) wblog(FL,
      "<i> time range  : %.2g .. %.2g (%d vals)",
      t0[0], t0.last(), t0.len
   );

   wblog(FL,"<i> %47R","=");

   if (!CC.isEmpty()) { nloc=0;

      if (zz.len || !Z0.isEmpty()) wblog(FL, 
         "WRN %Nsince C is specified as cell array,\n"
         "zflags and Z0 will be ignored!%N"
      );

      if (!DD.len) {
         for (k=i=0; i<CC.len; i++) {
            nloc=MAX(nloc,(unsigned)CC[i].len);
            if (CC[i].len) { if (k<i) CC[k]=CC[i]; k++; }
         }
         if (k<i) { wblog(FL,
            "WRN empty operator sets will be ignored (%d/%d)",i-k,i);
            CC.len=k;
         }
      }
      else {
         if (DD.len!=CC.len) wblog(FL,"ERR C2 since specified "
         "must have the same length as C1 (%d/%d)", DD.len, CC.len);

         for (k=i=0; i<CC.len; i++) {
            nloc=MAX(nloc, (unsigned)MAX(CC[i].len, DD[i].len) );
            if (CC[i].len && DD[i].len) {
               if (k<i) { CC[k]=CC[i]; DD[k]=DD[i]; }
               k++;
            }
         }
         if (k<i) { wblog(FL,
            "WRN empty operator sets will be ignored (%d/%d)",i-k,i);
            CC.len=k; DD.len=k;
         }
      }

      C1.initOp(FL,CC, NULL,"1st");
      C2.initOp(FL,DD, &C1, "2nd");
   }
   else { 
      if (D0.len && D0.len!=C0.len) wblog(FL,
         "ERR C2 since specified must have same length as C1 (%d/%d)",
          D0.len, C0.len
      );

      C1.initOp(FL,C0  ); C1.initSym("KK");
      C2.initOp(FL,D0  ); C2.initSym("KK");

      if (calcOps) { C1.forceCalc(); C2.forceCalc(); }

      if (!zz.data) { zz.init(C0.len);  }
      else {
         if (zz.len!=C0.len) {
            if (zz.len==1) { zz.Resize(C0.len); zz.set(zz[0]); }
            else wblog(FL,
               "ERR Length mismatch of zflags with C1 (%d/%d)",
                zz.len, C0.len
            );
         }
      }

      if (zz.anyUnequal(0)) {
         if (Z0.isEmpty()) wblog(FL,
         "ERR Z0 not defined (required since zflags is specified!)");

         if ((C1.calc>0 && C1.KK.len) || (C2.calc>0 && C2.KK.len))
         wblog(FL," *  applying Z0 [zflags = %s]", zz.toStr().data);

         for (k=0; k<2; k++) {
            NRGData<gTQ,gTD>  &Ck = (k==0 ? C1 : C2); e=0;
            if (Ck.calc>0) {
               for (i=0; i<Ck.KK.len; i++) if (zz[i])
               e+=applyZ0(Z0,Ck.KK[i]);
            }
            if (e) wblog(FL,
               "ERR QIDX of Z0 does not match QIDX of %s\n"
               "Hint: unset Z0 to empty if not relevant.",
                Ck.name.data
            );
         }
      }
   }

   if (calcOps) { C1.forceCalc(); C2.forceCalc(); } 
   if (!store ) {
      C1.store=0;
      if (C2.nrgIdx.len) C2.store=0;
   }

   if (C2.KK.len==0) C2.store=0;

   C1.info(FL,"C1"); if (C2.KK.len) {
   C2.info(FL,"C2"); }

   if (1 || C1.store) {
      if (CC.isEmpty())
           C1.updateInfop(FL,"i",C0.toMx());
      else C1.updateInfop(FL,"i",C1.CI.toMx());
   }
   if (1 || (C2.store && C2.CI.len)) {
      if (DD.isEmpty())
           C2.updateInfop(FL,"i",D0.toMx());
      else C2.updateInfop(FL,"i",C2.CI.toMx());
   }

   if (!C1.calc && !C2.calc) wblog(FL,
   "<i> using stored operator sets (%s,%s).", C1.name.data, C2.name.data);

   if (fgrFlag) {
      if (C2.KK.len)
         wblog(FL,"NB! calculating t-dependence in FGR modus");
      else {
         wblog(FL,"WRN FGR flag will be ignored (C2 is empty)");
         fgrFlag=0;
      }
   }

   tdmData.init(FL,"TDM", t0,
      N, C1.KK.len, emin, emax, nlog
   );

   if (RAW) tdmData.raw=1;

   nrgGetE0(DE,gES,E0); 

   E0-=E0[NRho>0 ? NRho-1 : N-1];
   initRHO(rhoNorm, T, H0, E0, dloc, NRho);

   dd.init(RHO.getMxInfo("rhoNorm"));

   if (calcRho || dd!=rhoNorm || RHO.checkVec(FL,"",2) ||
      (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T)
   ){
      wblog(FL,"<i> update density matrices\\");

      if (locRho && !RHO.MX) {
         printf(" internally (locRho=%d)",locRho);
         RHO.MX=mxCreateStructMatrix(1,NRG_N,0,NULL);
         RHO.locMX=1; locRho=99;
      }

      printf("\n\n");
      nrgUpdateRHO(rhoNorm, A0,H0, RHO, NRho, dloc, &SIG);
      printf("\n");
   }
   else wblog(FL,"<i> using stored DM set.%N");

   if (locRho!=99) locRho=0;

   Rho.initDef(2); P.initOp(1);

   for (iter=0; iter<N; ++iter) {
      NRG_ITER=iter; SIG.check911();

      rhoNorm.getRec(iter,w); 
      wx+=(fabs(w[0])+fabs(w[1]));

      A0.init(FL,"K",iter);  H0.init(FL,"K",iter);
      A0.init(FL,"D",iter);  H0.init(FL,"D",iter); 
      A .init(FL,"K",iter);  H .init(FL,"K",iter);
      A .init(FL,"D",iter);  H .init(FL,"D",iter);

      snprintf(itag,16,"%s %2d/%d", vtag, iter+1, N); lflag=0;

      C1.updateOp(FL,A,4,iter);    C1.saveOp(FL,iter,4);
      C2.updateOp(FL,A0,A,4,iter); C2.saveOp(FL,iter,4);

      if (iter) {
         P.updateOp(FL,A0,A,4,iter); 
         P.skipOpZeros(DEPS,'b');
      } else { initOverlap(P,A0,A); }

      if (wx>DEPS) { 
         if (w[0]>DEPS) {
            if (w[1]>DEPS) wblog(FL,
            "%s kept (%.3g) + truncated (%.3g)",itag,w[0],w[1]);
            else wblog(FL, "%s kept (%.3g)",itag,w[0]);
         }
         else if (w[1]>DEPS)
         wblog(FL,"%s truncated (%.3g)%20s\r\\",itag,w[1],"");
         lflag++;

         for (i=0; i<2; ++i) if (w[i]) {
            const QSpace<gTQ,gTD> &HX = (i==0 ? H0.K : H0.D);

            initRho(Rho[i], HX, iter, w[i]);
            tdmNRGIter(
               C1,C2,P, Rho[i],(i==0 ? 'K':'D'), H, tdmData, iter,
               "includeKK"
            );
         }
      }

      RHO.init(FL,"K",iter);
      tdmNRGIter(C1,C2,P, RHO.K,'K', H, tdmData, iter);

      if (!lflag) {
         if (disp>2) wblog(FL,"%s (RHO)",itag);
         else wblog(FL,"%s (RHO)%30s\r\\",itag,"");
      }
   }

   if (disp<=2) printf("\n\n");

   if (RAW) {
      tdmData.addPartialFourier();
      ar=tdmData.At.toMx('r');
   }

   tdmData.getSpecData(om,a0);

   if (raw) {} 
   else if (sigma>0 && t0.len) {
      istr.printf(FL,"frequency broadened TDM data (sigma=%g)", sigma);

      tdmData.getSmoothSpec_t(sigma, 0.); 
      tdmData.Fourier(om,aa,t0,at,alpha);
   }
   else if (alpha>0) {
      if (Lambda<=1) wblog(FL,"ERR alpha requires Lambda (%g)",Lambda);
      istr.printf(FL,"time-domain broadened TDM data (alpha=%g)",alpha);

      tdmData.Fourier(om,at,alpha,Lambda);
   }
   else if (t0.len) {
      istr.printf(FL,"plain Fourier transformed data (no broadening)");
      wblog(FL,"<i> %s",istr.data);

      Wb::Fourier(om,a0,t0,at);
   }

   S=mxCreateStructMatrix(1,1,0,NULL);
   mxAddField2Scalar(FL,S,"ver",    wbstring(vstr).toMx());

   mxAddField2Scalar(FL,S,"T",      numtoMx(T));
   mxAddField2Scalar(FL,S,"alpha",  numtoMx(alpha));
   mxAddField2Scalar(FL,S,"sigma",  numtoMx(sigma));
   mxAddField2Scalar(FL,S,"rho",    rhoNorm.toMx());
   mxAddField2Scalar(FL,S,"E0",     E0.toMx());
   mxAddField2Scalar(FL,S,"EScale", gES.toMx());

   mxAddField2Scalar(FL,S,"om",     om.toMx('t'));
   mxAddField2Scalar(FL,S,"a0",     a0.toMx('r'));
   mxAddField2Scalar(FL,S,"aa",     aa.toMx('r'));

   if (ar)
   mxAddField2Scalar(FL,S,"cc_raw", ar);

   if (partial)
   mxAddField2Scalar(FL,S,"RAW",    tdmData.toMx());
   mxAddField2Scalar(FL,S,"emin",   numtoMx(emin));
   mxAddField2Scalar(FL,S,"emax",   numtoMx(emax));
   mxAddField2Scalar(FL,S,"nlog",   numtoMx(nlog));

   { mxArray *F=mxCreateCellMatrix(1,2);
     mxSetCell(F,0,RHO.NAME.toMx());
     mxSetCell(F,1,C1 .NAME.toMx());
     mxAddField2Scalar(FL,S,"NRG",F);
   }

   mxAddField2Scalar(FL,S,"stamp",  Wb::TimeStamp().toMx());

   if (!C1.NAME.isEmpty()) {
      C1.updateInfo(FL, "tt",      t0.toMx('t'));
      C1.updateInfo(FL, "cc",      at.toMx('t'));
      C1.updateInfo(FL, "Itdm",    S, 'k'); 
   }

   if (!RHO.NAME.isEmpty()) {
      RHO.updateInfo(FL, "rhoT",    numtoMx(T));
      RHO.updateInfo(FL, "rhoNorm", rhoNorm.toMx());
   }

   argout[0]=at.toMx('t'); 

   if (nargout>1)
        argout[1]=S;
   else mxDestroyArray(S);

#ifdef MATLAB_MEX_FILE
   myCleanUp();
#endif

   Wb::Clock *clio=fdmClocks.get("data:I/O",1);
   wblog(FL,"I/O %s time usage: %s%N", myname,
      clio ? SEC2STR(clio->gettime()) : "0"); 
   if (clio) { clio->reset(); }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in tdmNRG"); }
};

template <class TQ, class TD>
void initOverlap(
   NRGData<TQ,TD> &P,
   const NRGData<TQ,TD> &A0,
   const NRGData<TQ,TD> &A
){
   unsigned i,j,k=0;
   const char* tags[4] = { "KK","KD","DK","DD" };

   const QSpace<TQ,TD> *a0[] = { &A0.K, &A0.D };
   const QSpace<TQ,TD> *aa[] = { &A .K, &A .D };

   for (i=0; i<2; ++i)
   for (j=0; j<2; ++j) {
      wbvector< QSpace<TQ,TD> > &Pk=P.getOpsXX(tags[k++]);
      if (Pk.len!=1) Pk.initDef(1);
      a0[i]->contract("13*",*(aa[j]),"13",Pk[0]); 
   }
}

template <class TQ, class TD, class TS>
void tdmNRGIter(
   NRGData<TQ,TD> &C1, 
   NRGData<TQ,TD> &C2,
   NRGData<TQ,TD> &P,
   QSpace<TQ,TD> &Rho, char w, 
   NRGData<TQ,TD> &H,
   TDSpectral<TS> &tdmData,
   int iter,
   const char* KKflag 
){
   unsigned i,k,n,m=1;
   const char* tags[4] = { "KK","KD","DK","DD" };
   QSpace<TQ,TD> R1,R2;

   if (iter<0) iter=NRG_ITER;

   if (KKflag && KKflag[0]!=0) {
      if (strcmp(KKflag,"includeKK"))
      wblog(FL,"ERR invalid KKflag=%s !?",KKflag);
      m=0;
   }

   n=C1.KK.len;
   if (n!=C1.KD.len || n!=C1.DK.len || n!=C1.DD.len) wblog(FL,
   "ERR severe operator dim mismatch [%d %d %d %d]",
   C1.KK.len, C1.KD.len, C1.DK.len, C1.DD.len); 

   for (i=m; i<4; i++) {
      const QSpace<TQ,TD>
         &H1 = H.getQSpace(tags[i][0]),
         &H2 = H.getQSpace(tags[i][1]);

      for (k=0; k<n; k++) {
         const QSpace<TQ,TD>
         &Ck = C1.getQSpace(tags[i],k),
         &P1 =(C2.KK.len ? C2.getQSpace(w,tags[i][1],k) 
                          : P.getQSpace(w,tags[i][1])),
         &P2 = (fgrFlag  ? C2.getQSpace(w,tags[i][0],k) 
                         :  P.getQSpace(w,tags[i][0]));

         P2.contract(FL,"1*",Rho,"1",R1);
         R1.contract(FL,2,P1,1,R2);

         R2.TimesEl(Ck,0,2); if (R2.isEmpty()) continue;

         updateSpec(FL, tdmData, iter, k, R2, H1, H2);
      }
   }
}

template <class TQ, class TD, class TS>
void updateSpec(
   const char *F, int L,
   TDSpectral<TS> &tdmData,
   unsigned iter, unsigned io,
   QSpace<TQ,TD> &RX, 
   const QSpace<TQ,TD> &H1,
   const QSpace<TQ,TD> &H2
){
   unsigned i,j,k,m,n;
   wbvector<unsigned> I1,I2;
   wbindex i1h,i2h;

   wbvector<TD> E1,E2;
   wbarray<TD> EE;
   TD *e1,*e2;

   matchH12(FL,RX,H1,H2,iter, i1h,E1,I1, i2h,E2,I2,'r');

   for (k=0; k<i1h.len; k++) {
       i=i1h[k]; m=I1[i+1]-I1[i];
       j=i2h[k]; n=I2[j+1]-I2[j]; EE.init(m,n);

       if (!RX.DATA[k]->sameSize(EE)) { 
          RX.DATA[k]->info("RX[k]"); EE.info("EE");
          wblog(FL,"WRN dimension mismatch!? (%d)",k);
       }

       e1=E1.data+I1[i]; e2=E2.data+I2[j];

       for (i=0; i<n; ++i)
       for (j=0; j<m; ++j) EE(j,i) = e1[j]-e2[i];

       tdmData.Add(EE, *(RX.DATA[k]), iter, io);
   }
};

