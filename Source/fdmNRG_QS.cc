/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace NRG routines
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

/* CHANGE LOG ======================================================== *

   changed (C1,C2) or (F,C) to (B,C) consistently
   Wb,Nov24,14

 * =================================================================== */

char USAGE[]=""; // outsourced to fdmNRG_QS.m // Wb,Jan12,19

char USAGE_2[] =
"Alternative usage: fdmNRG('nrgdata', 'setup.mat', 'B', 'C')";

#define PROG_TAG "FDM"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "fdmNRG"
#endif

#include <float.h>    

#define TST_TCHI_0

#define  EPS 1E-12
#define DEPS DBL_EPSILON

#define LD_CLEBSCH_QS
#include "wblib.h"

#include "nrgdata.hh"  
#include "dmrho.cc"    
#include "spectral.hh"

template <class TQ, class TD>
void dmNRGIter(
    const char *F_, int L,
    const unsigned iter,
    Spectral<TD> &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &B,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf=0
);

template <class TQ, class TD>
inline void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral<TD> &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf=0
){
    dmNRGIter(F,L,iter, ASpec, Rho, H, C, C, cc, XX, isbuf);
}

template <class TQ, class TD>
void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral<TD> &ASpec,
    const QSpace<TQ,TD> &Rho,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,TD> &HD,
    const wbvector< QSpace<TQ,TD> > &BKT,
    const wbvector< QSpace<TQ,TD> > &BTK,
    const wbvector< QSpace<TQ,TD> > &CKT,
    const wbvector< QSpace<TQ,TD> > &CTK,
    const wbvector<char> &cc,
    const char *isbuf=0,
    char dblock=0
);

template <class TQ, class TD>
void updateSpec(
    const char *F, int L,
    Spectral<TD> &ASpec, unsigned iter, const unsigned is,
    QSpace<TQ,TD> &RX,
    const QSpace<TQ,TD> &H1,
    const QSpace<TQ,TD> &H2,
    const char *isbuf=0,
    char dblock=0,
    char mflag=0
);

mxArray* addTOA(const char *F, int L, const mxArray *S,
   double T,
   const wbvector<double> &om,
   const wbMatrix<double> &A0,
   const char *Fname="",
   const char *Cname=""
);

   static void myCleanUp(void) { Wb::ResSummary(FL); };

   double Lambda=1;

template<class TQ, class TD>
void FDM_NRG(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[], char mat_setup
);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   char cmplx=0, mat_setup=0;  int i=0;
   const mxArray *ap=NULL;

   NRGData<gTQ,double> XD("XD");

   MX_CHECK_HELPER_NARGS(2,-1,-1); 

   if (mxIsChar(argin[0])) {
      mxGetString(argin[0],str,8); str[7]=0;
      if (!strcmp(str,"--mat")) { mat_setup=1; --nargin; ++argin; }
   }

   if (mat_setup) {
      if (nargout || (nargin!=3 && nargin!=4)) wblog(FL,
         "ERR invalid number of I/0 arguments (%d/%d)",nargout,nargin);

      for (; i<nargin; ++i) {
         if (!mxIsChar(argin[i])) wblog(FL,
         "ERR invalid arg #%d %N%N%s%N",i+1,USAGE_2);
      }
   }
   else if (mxIsChar(argin[0])) { ++i; }
   else {
      if (!mxIsStruct(argin[1])) wblog(FL,"ERR %s() invalid "
         "in-memory usage\n(expecting struct Inrg for arg #2)",FCT);
      ap=mxGetField(argin[1],0,"paras"); 
   }

   if (i) {
      wbstring fname(128);

      i=mxGetString(argin[0],fname.data,fname.len-1); 
      if (i || !fname) wblog(FL,
         "ERR reading NRGdata string from arg #1 (l=%d)",fname.len);
      XD.setupIO(FL,fname.data,NULL);
      ap=XD.getMxInfo("paras"); 
   }

   if (ap) {
      double q=0; 
      const mxArray *ac=mxGetField(ap,0,"complex");
      if (!ac) { cmplx=-2; }
      else { mxGetNumber(ac,q); if (q) { cmplx=1; } }
   }
   else { cmplx=-1; }

   if (cmplx<0) wblog(FL,
      "WRN %s() failed to determine cmplx (%d)",myname,cmplx);

   if (cmplx>0)
        { FDM_NRG<gTQ,wbcomplex>(nargout,argout,nargin,argin,mat_setup); }
   else { FDM_NRG<gTQ,double>   (nargout,argout,nargin,argin,mat_setup); }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in fdmNRG_QS"); }
   aclu.Check();
}; 

template<class TQ, class TD>
void FDM_NRG(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[], char mat_setup
) {

   unsigned i,e,k,m,iter,ilast=0, N=0, lflag, NRho=0, dloc=0, Dk;
   int vflag=Wb::envVRB; 
   char store=1, gotops=1, mspec=0, gotwk=0;
   char noRHO=0, calcRho=0, locRho=0, calcOps=0, version=4;
   char rzero, keven, kodd, Tflag=0, partial=0, PARTIAL=0;
   double TN, Om=0, dbl=0, Delta=0., T=-1;
   const unsigned r1=2, r2=4;

   wbvector<unsigned> D;
   wbvector<double> w;

   Spectral<TD> ASpec;
   wbvector<TD> a1, a2; 
   wbMatrix<TD> A0;

   mxArray *a;
   const mxArray *ac, *Inrg=NULL;

   wbstring vstr(32), fname(128);
   char vtag[4], vbuf[4], itag[16];

   NRGData<TQ,TD> H("H"), A("A"), RHO("RHO"), B("B"), C("C");
   wbvector< QSpace<TQ,TD> > B0, C0;
   QSpace<TQ,TD> Rho, Z0;

   wbMatrix<double> rhoNorm, rhoNormI, dd;
   wbvector<double> om,DE,E0;
   wbvector<char> zz,cc;
   OPTS opts;

   Wb::SigHandler SIG(FL); 

   NRG_N=0; NRG_ITER=0; gES.init();

   Wb::Clock *clt;
   Wb::Clock clt_("FDM:all",0,0,&fdmClocks,&clt); 

#ifdef MATLAB_MEX_FILE
   mexAtExit(myCleanUp);  
#endif
   Wb::ResSummary(NULL,1);

   if (mat_setup) {

      i=mxGetString(argin[0],fname.data,fname.len-1);
      if (i || fname.isEmpty()) wblog(FL,
         "ERR reading string of arg #1 (%d)", fname.len);

      if (mxGetString(argin[1],str,STRLEN)) wblog(FL,
         "ERR error reading string of arg #2 (%d)", STRLEN);
      opts.init(str); 

      if (nargin==4) {
         k=2; 
         if (mxGetString(argin[k], str, STRLEN)) wblog(FL,
            "ERR error reading string of arg #%d",k+1);
         B.initName(str); 
         opts.getOpt(str,a,'!'); if (a) mxInitQSpaceVec(FL,a,B0,0,&r1,&r2);
      }

      k=nargin-1;
      if (mxGetString(argin[k],str,STRLEN)) wblog(FL,
         "ERR error reading string of arg #%d", k+1);
      C.initName(str); 
      opts.getOpt(str,a,'!'); if (a) mxInitQSpaceVec(FL,a,C0,0,&r1,&r2);

      opts.getOpt("Z0",a); if (!a) opts.getOpt("Z",a,'!');
      if (a) { unsigned r=-1;
         if (mxIsQSpace(FL,a,r)<=0) wblog(FL,"ERR invalid Z0 (%s)", str);
         Z0.init(FL,a);
      }
   }
   else {
      if (nargout>4 || nargin<4) wblog(FL,
         "ERR invalid number of I/O arguments (%d/%d)",nargout,nargin);

      if (mxIsChar(argin[0])) {
         i=mxGetString(argin[0], fname.data, fname.len-1);
         if (i || fname.isEmpty()) wblog(FL,
            "ERR error reading file tag from arg #1");
      }

      ac=argin[1];
      if (!ac || mxGetNumberOfElements(ac)!=1
         || mxGetFieldNumber(ac,"Lambda")<0
         || mxGetFieldNumber(ac,"EScale")<0
         || mxGetFieldNumber(ac,"E0"    )<0
      ){
         k=1; if (fname.isEmpty())
         wblog(FL,"ERR valid Inrg structure required (2nd argument)");
      }
      else { Inrg=ac; k=2; }

      for (i=1; i<=2; ++i, ++k) {
         if (mxIsEmpty(argin[k])) {
            continue;
         }

         if (mxIsChar(argin[k])) { try {
            if (mxGetString(argin[k],str,STRLEN))
            wblog(FL,"ERR error reading string of arg #%d",k+1);

            ac=mexGetVariablePtr("caller",str); 

            if (!ac) wblog(FL,
              "ERR reading QSpace `%s' (arg #%d) from workspace",str,k+1);

            if (i==1)
                 { mxInitQSpaceVec(FL,ac,B0,0,&r1,&r2); B.initName(str); }
            else { mxInitQSpaceVec(FL,ac,C0,0,&r1,&r2); C.initName(str); }
            }
            catch (...) {
               wblog(FL,"ERR reading %s operator (arg #%d)",
               i==1? "1st":"2nd", k+1);
            }
         }
         else { try {
            wbvector< QSpace<TQ,TD> > &X=(i==1 ? B0 : C0);
            mxInitQSpaceVec(FL,argin[k],X,0,&r1,&r2);
         } catch (...) {
            wblog(FL,"ERR reading %d operator (arg #%d)",
            i==1? "1st":"2nd", k+1);
         }}
      }
      Z0.init(FL,argin[k++]);

      opts.init(argin+k, nargin-k);
   }

   if (opts.getOpt("-q")) vflag=0; else
   if (opts.getOpt(FL,"-v")) vflag|=2; else
   if (opts.getOpt(FL,"-V")) vflag|=4;

   printf("\n");
   if (!fname.isEmpty()) {
   }
   else {
      wblog(FL,"<i> internal mode (no file I/O)");
      if (nargout<4) wblog(FL,
         "ERR internal mode requires 4 output args (%d)",nargout);
      argin[0]=argout[3]=mxDuplicateArray(argin[0]); 
   }

   A  .setupIO(FL, fname.data, argin[0], vflag ? "" : 0); 
   H  .setupIO(FL, fname.data, argin[0]);
   RHO.setupIO(FL, fname.data, argin[0]);
   B  .setupIO(FL, fname.data, argin[0]);
   C  .setupIO(FL, fname.data, argin[0]);

   if (fname.isEmpty()) {
      if (!Inrg) { opts.getOpt("Inrg", a); Inrg=(const mxArray*)a; }
      if (!Inrg) wblog(FL,"ERR Inrg structure required (2nd argument)\n"
         "since no NRG file structure was specified");
      if (mxGetNumberOfElements(Inrg)!=1) wblog(FL,
         "ERR invalid input Inrg (scalar structure expected)");
      wblog(FL," *  reading Lambda, E0, EScale from Inrg structure");

      ac=mxGetField(Inrg,0,"Lambda"); if (ac) mxGetNumber(ac,Lambda);
      else wblog(FL,"ERR invalid input Inrg (missing field 'Lambda')");

      ac=mxGetField(Inrg,0,"E0"); if (ac) DE.init(FL,ac);
      else wblog(FL,"ERR invalid input Inrg (missing field 'E0')");

      ac=mxGetField(Inrg,0,"EScale"); if (ac) gES.init(FL,ac);
      else wblog(FL,"ERR invalid input Inrg (missing field 'EScale')");
   }
   else {
      if (vflag) {
         for (i=fname.len; i>0; --i) { if (fname[i]=='/') { ++i; break; }}
         wblog(FL," *  reading Lambda, E0, EScale from %s_info.mat",
         fname.data+i);
      }

      ac=RHO.getMxInfo("Lambda"); if (ac) mxGetNumber(ac,Lambda);
      else wblog(FL,"ERR failed to read Lambda from <NRG>_info.mat");

      ac=RHO.getMxInfo("E0"); if (ac) DE.init(FL,ac);
      else wblog(FL,"ERR failed to read E0 from <NRG>_info.mat");

      ac=RHO.getMxInfo("EScale"); if (ac) gES.init(FL,ac);
      else wblog(FL,"ERR failed to read EScale from <NRG>_info.mat");
   }

   if (WbUtil<TD>::isComplex()) wblog(FL,  
      " *  complex mode (based on Inrg)"); 

   if (Lambda<=1.) wblog(FL,"ERR invalid Lambda=%g.",Lambda);

   H.checkVec(FL,"K",2,"NRG data"); 

   N=NRG_N;

   if (!N) wblog(FL,"ERR empty NRG data set (N=%d)",N);
   if (N!=DE.len || N!=gES.len) wblog(FL,
      "ERR mismatch in chain length (E0[%d], EScale[%d], N=%d)",
       DE.len,gES.len, N);

   SIG.check911();
   dloc=getDLoc(A,&SIG,&Dk, 0);

   A.init(FL,"K",0); 

   TN=gEScale(N-1);  
   ASpec.Emin=TN*1E-3;

   store =!opts.getOpt("nostore"); 
   locRho =opts.getOpt("locRho" );
   noRHO  =opts.getOpt("noRHO"  ); 
   calcRho=opts.getOpt("calcRho"); if (locRho) calcRho=1;
   calcOps=opts.getOpt("calcOps");

   opts.getOpt("nlog", ASpec.nlog);
   opts.getOpt("emin", ASpec.Emin);
   opts.getOpt("emax", ASpec.Emax);
   opts.getOpt("binOffset",Delta);
   opts.getOpt("NRho", NRho);

   keven=opts.getOpt("keven");
   kodd =opts.getOpt("kodd" ); m=0;

   if (opts.getOpt("bul")) { ++m; version=0; } 
   if (opts.getOpt("hof")) { ++m; version=1; } 
   if (opts.getOpt("fDM")) { ++m; version=4; } 

   if (m>1) wblog(FL,
   "ERR multiple definitions of methods (%d)",m);

   opts.getOpt(FL,"mspec",mspec);

   opts.getOpt("cflags",a); if (a) { cc.init(FL,a,"cflags"); }
   opts.getOpt("zflags",a); if (a) { zz.init(FL,a,"zflags"); }

   if (NRho) {
      if (NRho>N) wblog(FL,"ERR NRho out of bounds (%d/%d)", NRho, N);
      if (version>3) version=3;
   }
   else if (version<4) NRho=N;

#ifndef MATLAB_MEX_FILE
   if (Wb::GetEnv(FL,"FDM_T",T)>=0)
   Tflag=1;
#endif

   opts.getOpt("T", T);
   partial=opts.getOpt("partial");
   PARTIAL=opts.getOpt("PARTIAL");

   if (opts.getOpt("rhoNorm",a)) { rhoNormI.init(FL,a); }

   opts.checkAnyLeft(); 

   if (ASpec.Emin<=0) {
      ASpec.Emin=1E-3*pow(Lambda,-double(N)/2);
   }

   if (1) {
      double e1=0., e2=0.;

      H.init(FL,"K",0); e1=getErange(H.K)*gEScale(0); 
      H.init(FL,"K",1); e2=getErange(H.K)*gEScale(1); dbl=MAX(e1,e2);

      if (ASpec.Emax<=0) {
         if (dbl<10) { ASpec.Emax=10; }
         else { ASpec.Emax=10*ceil(0.12*dbl); }
      }
      else { e1=ASpec.Emax;
         if (e1<10) { ASpec.Emax=10; } else
         if (e1<1.2*dbl) {
            ASpec.Emax=10*ceil(0.12*dbl);
         }

         if (e1!=ASpec.Emax) wblog(FL,
         "WRN adjusting Emax=%g => %g (having E<=%.3g)",e1,ASpec.Emax,dbl);
      }
   }

   if (vflag) {
      strcpy(str,myname);
      wblog(FL,"=== %s",Wb::strpad(str,'=',48)); str[0]=0;
      if ((vflag &0xf) !=1) sprintf_str(", vflag=%d", vflag);
      if (!store )    strcat(str,", nostore");
      if (locRho )    strcat(str,", locRho");
      if (calcOps)    strcat(str,", calcOps");
      if (version==0) strcat(str,", bul"); else
      if (version==1) strcat(str,", hof");
      if (keven  )    strcat(str,", keven");
      if (kodd   )    strcat(str,", kodd");
      if (partial)    strcat(str,", partial");
      if (PARTIAL)    strcat(str,", PARTIAL");

      wblog(FL," *  Lambda=%g, L=%d, D=%d, dloc=%d\nsym=%s",
         Lambda,N,Dk, dloc, A.K.qtype.toStr('V').data);
      if (str[0]) { wblog(FL," *  %s",str+2); }
      fflush(0);
   }

   switch (version) {
     case 0:

       vstr.cpy(FL,"bare NRG (R. Bulla)");
       strcpy(vtag,"NRG"); break;

     case 1:

       if (keven || kodd)
       wblog(FL,"WRN keven/kodd with conventional fdmNRG !?");

       if (NRho==N)
            vstr.cpy(FL,"original DM-NRG (W. Hofstetter)");
       else vstr.printf(FL,"DM-NRG (W. Hofstetter; N=%d/%d)",NRho,N);
       strcpy(vtag,"HOF"); break;

     case 3:

       vstr.printf(FL,"fdm-NRG (single shell NRho=%d)", NRho);
       strcpy(vtag,"DMK"); break;

     case 4:

       vstr.cpy(FL,"fdm-NRG (full density matrix)");
       strcpy(vtag,"FDM"); break;

     default: wblog(FL,"ERR invalid version %d.", version);
   }

   if (version<=1) strcpy(vbuf,"buf"); else vbuf[0]=0;

   if ((keven || kodd) && version>1) {
      wblog(FL,"WRN %s ignored unless using conventional NRG\n"
        "hof or bul; using %s (%d)", keven ? "keven":"kodd",
         vstr.data, version
      ); keven=kodd=0;
   }
   else if (keven && kodd) {
      wblog(FL,"WRN Both keven and kodd are set (ignoring kodd)");
      kodd=0;
   }

   if (version>=4) {
      dbl=sqrt(2/(log(Lambda)*log(double(dloc))));
      m=(unsigned)ceil(5*dbl); 

      k=( m+1<N ? N-m-1 : 0 );
      if (int(m)<0) { m=0; }; dbl=gEScale(k);
      if (T<0) { T=dbl;
         if (vflag) wblog(FL," *  using default T=%.4g (N-%d)", T, N-k);
      }
      else if (T>0 && T<dbl)
      wblog(FL,"WRN T=%.3g < T_{N-%d}=%.3g", T, N-k, dbl);
   }
   else if (T<0) T=0;

   if (vflag) {
      str[0]=0; if (Delta!=0) sprintf_str(", binOffset=%g",Delta);
      wblog(FL,"\b *  %48R\n"
      " *  method      : %s\n"
      " *  temperature : %.4g (= %.4g TN)\n"
      " *  omega range : %.4g .. %.4g (%d/dec%s)\n"
      "=== %48R", "-",
      vstr.data, T, T/TN, ASpec.Emin, ASpec.Emax, ASpec.nlog, str, "=");
      fflush(0);
   }

   SIG.check911();
   if (version>1) {
      A.checkVec(FL,"D",3,"NRG data"); SIG.check911();
      H.checkVec(FL,"D",2,"NRG data"); SIG.check911();
   }

   if (C0.isEmpty()) gotops=0;

   if (B0.len && B0.len!=C0.len) wblog(FL,
      "ERR length mismatch of 1st with 2nd operator (%d,%d).",
       B0.len, C0.len);
   if (B0.len && B0[0].QDIM!=C0[0].QDIM) wblog(FL,
      "ERR QDIM mismatch of 1st with 2nd operator (%d,%d).",
       B0[0].QDIM, C0[0].QDIM);

   C.initOp(FL, C0, version>1, NULL, "op2=C"); 
   B.initOp(FL, B0, version>1, &C,   "op1=B"); 

   if (calcOps){
      C.forceCalc(); 
      B.forceCalc();
   }
   if (!store) { C.store=0; if (B.nrgIdx.len) B.store=0; }

   if (vflag) {
      if (!gotops)
         wblog(FL," *  got empty ops: calculate FDM/Z only");
      else {
         if (B.nrgIdx.len)
         B.info(FL ); else {
            wblog(FL," *  B = C -> calculating <C(t)||C'>"); }
         C.info(FL );
      }
   }

   if (C.store) C.updatePara(FL,"KK", 0, "KK_"); 
   if (B.store) B.updatePara(FL,"KK", 0, "KK_");

   if (vflag) {
      if (!C.calc && !B.calc)
      wblog(FL,"<i> using stored operator set");
   }

   if (B.nrgIdx.len && B.KK.isEqual(C.KK)) {
      wblog(FL,"--> computing C only since B=C");
      B.nrgIdx.init();
   }

   if (!zz.len || !cc.len) {
      if (cc.len) { 
         zz=cc; for (i=0; i<zz.len; ++i) zz[i]=(cc[i] ? 0:1);
      }
      else if (zz.len) { 
         cc=zz; for (i=0; i<cc.len; ++i) cc[i]=(zz[i] ? 0:1);
      }
      else {
         zz.init(C0.len); cc.init(C0.len);
         if (Z0.isEmpty()) cc.set(1); else zz.set(1);
      }
   }

   if (zz.len!=C0.len) {
      if (zz.len==1) { zz.Resize(C0.len); zz.set(zz[0]); } else wblog(FL,
      "ERR length mismatch of zflags with [BC] (%d/%d)",zz.len,C0.len);
   }
   if (cc.len!=C0.len) {
      if (cc.len==1) { cc.Resize(C0.len); cc.set(cc[0]); } else wblog(FL,
      "ERR length mismatch of cflags with [BC] (%d/%d)",cc.len,C0.len);
   }

   if (vflag) {
      if (zz^cc) {
         if (zz.len==1) 
              wblog(FL," *  zflags=%d", zz[0]);
         else wblog(FL," *  zflags=[%s]", zz.toStrf("","").data);
      }
      else
         wblog(FL," *  zflags=[%s], cflags=[%s]",
         zz.toStrf("","").data, cc.toStrf("","").data
      );
   }

   if (!Z0.isEmpty()) { e=0;
      try {
         if (C.calc>0) {
            for (i=0; i<C.KK.len && !e; ++i)
            if (zz[i]) e+=applyZ0(Z0,C.KK[i]);
         }
         if (B.calc>0) {
            for (i=0; i<B.KK.len && !e; ++i)
            if (zz[i]) e+=applyZ0(Z0,B.KK[i]);
         }
         if (e) wblog(FL,"ERR %s() e=%d",FCT,e); 
      }
      catch(...) {
         wblog(FL,
           "ERR QIDX of Z0 does not match Q-spaces of %s\n"
           "Hint: unset Z0 to empty if not relevant.",
            e==1 ? C.name.data : B.name.data
         );
      }
   }
   else
   if (zz.anyUnequal(0)) wblog(FL,"ERR got zflags for empty Z0!");

   ASpec.init("A", 2*C.KK.len, Delta, T, PARTIAL ? N:0); if (mspec) {
   ASpec.calcSpectralMoments(FL,mspec); }

   if (version<=1) ASpec.initBuf();

   nrgGetE0(DE,gES,E0);

   if (version>=1) {
      double Eref=E0[NRho>0 ? NRho-1 : N-1];
      E0-=Eref;
      Om=initRHO(rhoNorm,T,H,E0,dloc,NRho,vflag) + Eref;
   }
   else {
      wbvector<double> Ek;
      H.init(FL,"K",0); getEData(H.K,Ek,D);
      E0.set(0); E0[0]=-Ek.min()*gEScale(0);

      initRHO(rhoNorm, T, H, E0, 0, NRho, vflag);
   }
   SIG.check911();

   if (!rhoNormI.isEmpty()) {
      if (rhoNorm.sameSize(rhoNormI)) {
         wblog(FL,"NB! using rhoNorm from input args");
         rhoNorm=rhoNormI;
      }
      else wblog(FL,
      "ERR size mismatch of input rhoNorm\n%dx%d vs. %dx%d",
       rhoNorm.dim1,rhoNorm.dim2,rhoNormI.dim1,rhoNormI.dim2);
   }
   if (rhoNorm.anyLT(0)) {
      MXPut(FL,"Ix").add(rhoNorm,"rhoNorm").add(T,"T");
      wblog(FL,"ERR %s() got negative weights !?",FCT);
   }
   gotwk=rhoNorm.col_anyGT(0,DEPS); 
   if (gotwk) wblog(FL,
      "NB! skipping RHO contributions in favor of RhoK");

   if (version && (gotops || !noRHO)) { wbtop().runningLarge(FL);

      dd.init(RHO.getMxInfo("rhoNorm"));

      if (calcRho || dd!=rhoNorm || RHO.checkVec(FL,"",2) ||
         (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T)
      ){
         int q=-1; char s[16]; s[0]=0;
            if (calcRho) q=1; else
            if (dd!=rhoNorm) q=2; else
            if (RHO.checkVec(FL,"",2)) q=3; else
            if (mxGetNumber(RHO.getMxInfo("rhoT"), dbl) && dbl!=T) q=4;

         if (locRho && !RHO.MX) {
            strcpy(s," internally");
            RHO.MX=mxCreateStructMatrix(1,NRG_N,0,NULL);
            RHO.locMX=1; locRho=99;
         }

         wblog(FL,"%s %scalculate density matrices RHO%s (q=%d)",
            vtag, !dd.isEmpty() ? "(re)":"",s,q); 

         SIG.check911();
         nrgUpdateRHO(rhoNorm, A, H, RHO, NRho, dloc, &SIG); 
         printf("\n"); SIG.check911();
      }
      else wblog(FL,"%s using stored density matrices",vtag);
   }

   wbtop().runningLarge(FL);

   if (locRho!=99) locRho=0;

   if (vflag & 248) { 
      wblog(FL,"TST CG_EPS     =[%g, %g] \r\\",CG_EPS1, CG_EPS2);
      wblog(FL,"TST CG_SKIP_EPS=[%g, %g] \r\\",CG_SKIP_EPS1, CG_SKIP_EPS2);
   }

   if (gotops)
   for (iter=0; iter<N; ++iter) { SIG.check911();

      NRG_ITER=iter;
      rhoNorm.getRec(iter,w); rzero=(w<DEPS);

      H.init(FL,"K",iter);
      H.init(FL,"D",iter); 

      if (C.calc>0 || B.calc>0) {
         A.init(FL,"K",iter);
         A.init(FL,"D",iter);
      }

      m=4; 

      C.updateOp(FL,A,m,iter); 
      B.updateOp(FL,A,m,iter);

      C.saveOp(FL,iter,m);
      B.saveOp(FL,iter,m);

      snprintf(itag,16,"%s %02d",vtag,iter+1); lflag=0;

      if (version<=1 && H.D.isEmpty()) {
         wblog(FL,"%s: all kept -> skip to next iteration%s",
         itag, vflag &12 ? "":"\r\\"); 
         continue;
      }

      if (version==0) {
         wblog(FL,"%s: (plain NRG)%-30s\r\\", itag,
         C.calc>0 || B.calc>0 ? " + operator update" : ""); ++lflag;
      }

      if (w[1]>=DEPS) { if (!ilast) { ilast=iter; }
         if (version) {
            if (vflag &12)
                 wblog(FL,"%s: (RhoT)", itag);
            else wblog(FL,"%s: (RhoT)%30s\r\\", itag,"");
            ++lflag;
         }

         initRho(Rho, H.D, iter, w[1]);
         if (B.nrgIdx.len==0) {
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "DD", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "DK", vbuf);
         }
         else {
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "DD", vbuf);
            dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "DK", vbuf);
         }
      }
      else if (H.D.isEmpty()) {
         if (vflag &12)
              wblog(FL,"%s: kept all states ", itag);
         else wblog(FL,"%s: kept all states %8s\r\\",itag,"");
         ++lflag;
      }

      if (w[0]>=DEPS) { 
         if (version) { ++lflag; wblog(FL,
            "%s: include RhoK @ w=%.4g", itag, w[0]); }

         initRho(Rho, H.K, iter, w[0]);

         if (B.nrgIdx.len==0)
              dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "KD", vbuf);
         else dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "KD", vbuf);

         if (!ilast || ilast==iter) {
         if (B.nrgIdx.len==0)
              dmNRGIter(FL,iter, ASpec, Rho, H,    C, cc, "KK", vbuf);
         else dmNRGIter(FL,iter, ASpec, Rho, H, B, C, cc, "KK", vbuf); }

         if (!ilast) { ilast=iter; }
      }

      if (rzero && !NRho && version!=1 && !H.D.isEmpty()) {
         if (vflag &12) 
              wblog(FL,"%s: (dRho=%.3g)",itag,w.sum());
         else wblog(FL,"%s: (dRho=%.3g)%8s\r\\",itag,w.sum(),"");
         ++lflag;
      }

      if (version>=1 && !gotwk) {
         RHO.init(FL,"K",iter);

         if (!lflag) {
            if (RHO.K.isEmpty()) {
               sprintf_str("%s: (EMPTY RHO)",itag);
            }
            else {
               sprintf_str("%s: (RHO)",itag);
               ++lflag;
            }

            if (vflag &12) wblog(FL,"%s",str);
            else wblog(FL,"%-50s\r\\",str);
         }

         if (!RHO.K.isEmpty()) {
         if (!H.D.isEmpty()) {
            if (B.nrgIdx.len==0)
                 dmNRGIter(FL,iter, ASpec, RHO.K, H,    C, cc, "KD", vbuf);
            else dmNRGIter(FL,iter, ASpec, RHO.K, H, B, C, cc, "KD", vbuf);
         }
         if (version==1) { 
            if (B.nrgIdx.len==0)
                 dmNRGIter(FL,iter, ASpec, RHO.K, H,    C, cc, "KK", vbuf);
            else dmNRGIter(FL,iter, ASpec, RHO.K, H, B, C, cc, "KK", vbuf);
         }}
      }

      if (vbuf[0]) {

         if ((!keven || (iter%2)==1) && (!kodd || (iter%2)==0)) {
            ASpec.crossAddBuf(2, iter<NRho && version);
         }
         else {
            sprintf_str(
            "%s: skipped (%s)", itag, keven ? "keep even" : "keep odd");
            if (vflag &12)
                 wblog(FL,"%s",str);
            else wblog(FL,"%s \r\\",str);
         }
         ASpec.initBuf(); 
      }

      if (PARTIAL) ASpec.saveIter(FL,iter);

      if (!lflag) {
         if (!B.calc && !C.calc) break;
         if (vflag &12)
              wblog(FL,"%s: (update operators)",itag);
         else wblog(FL,"%s: (update operators)%20s\r\\",itag,"");
         ++lflag;
      }
      doflush(); 
   }

   if (vflag &14) printf("\n"); 

   if (getSpecDimOps(B.DD.len ? B.DD :
      (B.KK.len ? B.KK : (C.DD.len ? C.DD : C.KK)), w))
   ASpec.applyIROPfac(FL,w);

   MXPut Iout(0,0);

   RHO.init(FL,"K",0);
   A.init(FL,"K",  0);

   Iout.add(vstr,"info").add(wbstring(vtag),"vtag")
      .add(Wb::TimeStamp(),"finished")
      .addP(MXPut(0,0)
         .add(B0,   "B"        ) 
         .add(C0,   "C"        ) 
         .add(zz,   "zflags"   )
         .add(cc,   "cflags"   )
         .add(E0,   "E0"       )
         .add(gES,  "EScale"   )
         .add(ASpec.Emin,"emin")
         .add(ASpec.Emax,"emax")
         .add(ASpec.nlog,"nlog")
         .add(Delta,"binOffset")
         .add(EPS,"eps").add(DEPS,"deps").add(CG_VERBOSE,"cg_verbose")
         .add(CG_EPS1,"cg_eps1").add(CG_EPS2,"cg_eps2")
         .add(CG_SKIP_EPS1,"cg_skip_eps1").add(CG_EPS2,"cg_skip_eps2")
        .toMx(),"paras")
      .add(A.K,"A0").add(RHO.K,"RHO") 
      .add(rhoNorm,"rho")
      .add(Om,"Om") 
      .add(T,"T")
      .add(w,"symfac");

   ASpec.getRawData(om,A0); 
   if (partial) Iout.addP(A0.toMx('r'),"A4");
   if (PARTIAL) {
      wbarray<TD> A1,A2;
      ASpec.getPARTIAL(A1,A2); Iout.add(A1,"A1").add(A2,"A2");
   }
   a1=A0.recSum();

   ASpec.pairRawData();
   ASpec.detailedBalance(T,om,A0);

   if (vflag)
   ASpec.dispSumRule();

   if (partial) { Iout.addP(A0.toMx('r'),"A4f"); }
   a2=A0.recSum();

   if (a1.len==a2.len) {
      A0.init(2,a1.len);
      A0.recSetP(0,a1.data);
      A0.recSetP(1,a2.data); Iout.add(A0,"a4");
   }
   else wblog(FL,"ERR size inconsistency (%d/%d)",a1.len,a2.len);

   Iout.add(ASpec.Ar,"reA0").add(ASpec.mspec,"mspec");

   a=A0.toMx();
   ASpec.getRawData(om,A0); 

   if (vflag) {
      doflush();
   }

   clt_.done();

   if (!fname.isEmpty()) {
      C.updateInfo(FL, "om",      om.toMx('t'));
      C.updateInfo(FL, "a0",      A0.toMx('r'));
      C.updateInfo(FL, "a4",      a);
      C.updateInfo(FL, "ISpec",   Iout.S,'k'); 

      if (!locRho) {
      C.updateInfo(FL, "rhoT",    numtoMx(T));
      C.updateInfo(FL, "rhoNorm", rhoNorm.toMx());
      }

      if (Tflag) {
         a=addTOA(FL,
           C.getMxInfo("TOA"), T, om, A0, B.name.data, C.name.data
         );
         C.updateInfo(FL,"TOA",a);
      }
   }

   if (nargout<=1) {
      argout[0]=A0.toMx('r'); 
   }
   else {
      if (!A0) { om.init(); } 
      argout[0]=om.toMx('t');
      argout[1]=A0.toMx('r');
   }

   if (nargout>2) { argout[2]=Iout.toMx(); }

   if (vflag) {
      doflush();
   }

   Wb::Clock *clio=fdmClocks.get("data:I/O",1);

   if (vflag) wblog(FL,"FIN NRGData I/O time usage: %s",
      clio ? SEC2STR(clio->gettime()) : "0");

   if (vflag &12 || CG_VERBOSE>5 || (clt && clt->gettime()>3600)) {
       if (clt ) { clt ->info(); }
       if (clio) { clio->info(); }
   }
   if (clt ) { clt ->reset(); }
   if (clio) { clio->reset(); }

   if (vflag) { myCleanUp(); printf("\n"); }

   if (vflag) {
      doflush();
   }

#ifdef __WB_MEM_CHECK__
   A0.init(); w.init(); a1.init(); a2.init();
   D.init(); fname.init();
   H.init(); A.init(); RHO.init(); Rho.init(); Z0.init();
   B.init(); C.init(); B0.init(); C0.init();
   ASpec.init(); opts.init(); gES.init();

   rhoNorm.init(); rhoNormI.init(); dd.init();
   om.init(); DE.init(); E0.init();
   zz.init(); cc.init();

   getBoltzman_base(E0,D,E0,-1,-1); 

   Wb::MemCheck(FL,"LIST");
#endif

   if (vflag) {
      doflush();
   }
};

mxArray* addTOA(const char *F, int L, const mxArray *S,
   double T,
   const wbvector<double> &om,
   const wbMatrix<double> &A0,
   const char *Fname,
   const char *Cname
){

   mxArray *a, *S2=mxCreateStructMatrix(1,1,0,NULL);
   const char* fn[] = {"T","Om","A0","B","C"};
   unsigned nf=5;

   if (!S) {
      mxAddField2Scalar(F,L, S2, fn[0], numtoMx( T ));
      mxAddField2Scalar(F,L, S2, fn[1], om.toMx('t'));
      mxAddField2Scalar(F,L, S2, fn[2], A0.toMx('r'));
      mxAddField2Scalar(F,L, S2, fn[3], mxCreateString(Fname ? Fname :""));
      mxAddField2Scalar(F,L, S2, fn[4], mxCreateString(Fname ? Cname :""));
   }
   else {
      unsigned i,j, m=mxGetM(S), n=mxGetN(S);
      int id[nf], e=mxGetNumberOfDimensions(S)>2 || (m!=1 && n!=1);

      m*=n;

      for (i=0; i<nf; ++i) {
         id[i]=mxGetFieldNumber(S,fn[i]);
         if (id[i]<0) ++e;
      }

      if (e) wblog(F,L,"ERR invalid structure TOA !?");

      n=mxGetNumberOfFields(S);
      S2=mxCreateStructMatrix(m+1,1,0,NULL);
      for (i=0; i<n; ++i) mxAddField(S2,mxGetFieldNameByNumber(S,i));

      for (i=0; i<m; ++i)
      for (j=0; j<n; ++j) {
          a=mxGetFieldByNumber(S,i,j); if (!a) continue;
          mxSetFieldByNumber(S2,i,j,mxDuplicateArray(a));
      }

      mxSetFieldByNumber(S2, i, id[0], numtoMx( T ));
      mxSetFieldByNumber(S2, i, id[1], om.toMx('t'));
      mxSetFieldByNumber(S2, i, id[2], A0.toMx('r'));
      mxSetFieldByNumber(S2, i, id[3], mxCreateString(Fname ? Fname :""));
      mxSetFieldByNumber(S2, i, id[4], mxCreateString(Fname ? Cname :""));
   }

   return S2;
}

template <class TQ, class TD>
inline void dmNRGIter(
    const char *F_, int L,
    const unsigned iter,
    Spectral<TD> &ASpec,
    const QSpace<TQ,TD> &Rho,
    const NRGData<TQ,TD> &H,
    const NRGData<TQ,TD> &B,
    const NRGData<TQ,TD> &C,
    const wbvector<char> &cc,
    const char *XX,
    const char *isbuf
){
    if (!strcmp(XX,"KK"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.K, H.K,
       B.KK, B.KK,
       C.KK, C.KK, cc, isbuf, 'K'
    ); else
    if (!strcmp(XX,"KD"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.K, H.D,
       B.KD, B.DK,
       C.KD, C.DK, cc, isbuf
    ); else
    if (!strcmp(XX,"DK"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.D, H.K,
       B.DK, B.KD,
       C.DK, C.KD, cc, isbuf
    ); else
    if (!strcmp(XX,"DD"))
    dmNRGIter(F_,L,iter, ASpec, Rho, H.D, H.D,
       B.DD, B.DD,
       C.DD, C.DD, cc, isbuf, 'D'
    );
    else wblog(FL,"ERR invalid tag `%s'", XX);
};

template <class TQ, class TD>
void dmNRGIter(
    const char *F, int L,
    const unsigned iter,
    Spectral<TD> &ASpec,
    const QSpace<TQ,TD> &Rho,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,TD> &HD,
    const wbvector< QSpace<TQ,TD> > &BKT,
    const wbvector< QSpace<TQ,TD> > &BTK,
    const wbvector< QSpace<TQ,TD> > &CKT,
    const wbvector< QSpace<TQ,TD> > &CTK,
    const wbvector<char> &cc,
    const char *isbuf,
    char dblock
){
   unsigned i, d=BKT.len, e=0;
   QSpace<TQ,TD> RX;

   if (HK.isEmpty() || HD.isEmpty()) {
      for (i=0; i<BTK.len; ++i) if (!BTK[i].isEmpty()) ++e;
      for (i=0; i<BKT.len; ++i) if (!BKT[i].isEmpty()) ++e;
      for (i=0; i<CTK.len; ++i) if (!CTK[i].isEmpty()) ++e;
      for (i=0; i<CKT.len; ++i) if (!CKT[i].isEmpty()) ++e;
      if (e) { dbstop(FL); wblog(FL,
      "ERR Inconsistency in empty data space (%d) !!", e); }

      return;
   }

   if (BTK.len!=d || CTK.len!=d || CKT.len!=d) wblog(FL,
      "ERR severe length mismatch (B-ops(%d): %d/%d; %d/%d)",
       iter+1, BKT.len, BTK.len, CKT.len, CTK.len);
   if (cc.len!=d) wblog(FL,
      "ERR cc length mismatch (%d/%d)", cc.len, d);

   if (Rho.isEmpty()) return;

   char conjA=1; 

   for (i=0; i<d; ++i) {
       Rho.contract(2,CKT[i],1,RX);
       RX.TimesEl(BKT[i],conjA,2);

       updateSpec(F_L, ASpec, iter, 
       2*i, RX, HK, HD, isbuf, dblock, cc[i]);
   }

   for (i=0; i<d; ++i) {
       wbperm P; if (CTK[i].isR3Op()) { P.initStr(FL,"132"); }

       CTK[i].contract(2,Rho,1,RX,P);
       RX.TimesEl(BTK[i],conjA,2);

       if (cc[i]) RX*=-1;

       updateSpec(F_L, ASpec, iter, 
       2*i+1, RX, HD, HK, isbuf, dblock, cc[i]);
   }
};

template <class TQ, class TD>
void updateSpec(
    const char *F, int L,
    Spectral<TD> &ASpec, unsigned iter, const unsigned is,
    QSpace<TQ,TD> &RX, 
    const QSpace<TQ,TD> &H1,
    const QSpace<TQ,TD> &H2,
    const char *isbuf,
    char dblock, 
    char mflag   
){
    unsigned i,j,k,m,n, tobuf=0;
    wbvector<unsigned> I1,I2;
    wbindex i1h, i2h;

    wbvector<double> E1, E2;
    wbarray<double> Omega, CC;
    double *e1, *e2;

    if (isbuf && isbuf[0]) {
       if (!strcmp(isbuf,"buf")) tobuf=1;
       else wblog(F,L,"ERR %s() invalid isbuf=>%s<",FCT,isbuf);
    }
    if (is>=ASpec.Ap.dim1) wblog(F,L,
       "ERR %s() index out of bounds (%d/%d)",FCT,is,ASpec.Ap.dim1);

    matchH12(FL,RX,H1,H2,iter, i1h,E1,I1, i2h,E2,I2,'r');

    for (k=0; k<i1h.len; ++k) {
        i=i1h[k]; m=I1[i+1]-I1[i];
        j=i2h[k]; n=I2[j+1]-I2[j]; Omega.init(m,n);

        if (!RX.DATA[k]->sameSize(Omega)) { 
           MXPut(FL,"i").add(*RX.DATA[k],"RX").add(Omega,"Omega").add(k+1,"k");
           RX.DATA[k]->info("RX[k]"); Omega.info("Omega");
           wblog(FL,"ERR %s() dimension mismatch (%d) !?",FCT,k);
        }

        e1=E1.data+I1[i]; e2=E2.data+I2[j];

        for (i=0; i<m; ++i)
        for (j=0; j<n; ++j) Omega(i,j) =  e2[j]-e1[i];

        if (tobuf) {
           ASpec.Add2buf(Omega, *(RX.DATA[k]), is); }
        else {
           ASpec.Add(
              Omega, *(RX.DATA[k]), is,
              dblock && RX.hasSameQ(k,0,k,1) ? dblock : 0,
              mflag
           );
        }
    }
};

