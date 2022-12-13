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

/* CHANGE LOG ======================================================== *

=> Wb,Apr22,10: added/updated deps and db to markSet

=> Wb,Apr13,13: keep A0 and H0 for iter=0 exactly as provided
   by input to avoid confusion when reading NRG_00.mat data.
   tags: KEEP_AH0

 * =================================================================== */

char USAGE[]=""; // outsourced to NRGWilsonQS.m // Wb,Jan12,19

#define PROG_TAG "NRG"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "NRGWilsonQS"
#endif

   #define WB_SPARSE_CLOCK

#define LOAD_CGC_QSPACE
#include "wblib.h"

   Wb::Clock nrgTime_0("NRG all"); 
   Wb::Clock nrgTime_5("NRG data I/O");

unsigned NRG_N, NRG_ITER;

template <class TQ, class TD>
void nrgBuildH4_abelian(
   QSpace<TQ,TD> &H4,
   const QSpace<TQ,TD> &HK,
   const QSpace<TQ,double> Es,
   const wbvector<TD> &f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FG
);

template <class TQ, class TD>
void nrgBuildH4_cg(
   QSpace<TQ,TD> &H4,
   QSpace<TQ,TD> &A4,
   const QSpace<TQ,TD> &HK,
   const wbvector<TD> &f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   char addHC,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FG
);

template <class TQ, class TD>
void updateFOps(
   wbvector< QSpace<TQ,TD> > &F12,
   const QSpace<TQ,TD> &A1,
   const QSpace<TQ,TD> &A2,
   const wbvector< QSpace<TQ,TD> > &FC,
   char lflag=0
);

template <class TQ, class TD>
void nrgDispIter(
   const unsigned iter,
   const unsigned itermax,
   const wbvector<TD> &E4, const TD &EK,
   const QSpace<TQ,TD> &HK,
   const char vflag
);

void checkGSDeg(const wbvector<double> &E4);

double getPhysE0(const wbvector<double> &DE);

   static void myCleanUp(void) { Wb::ResSummary(FL); };

   double Lambda=1;

   double FNfac=-1;

inline double nrgScale(int iter) { 

   if (FNfac<=0) wblog(FL,"ERR %s() invalid FNfac=%g",FCT,FNfac);

   double x = pow(Lambda, -double(iter)/2.) * 0.5*(Lambda+1) * FNfac;

   return (iter ? x : 1.0); 
};

template<class TQ, class TD>
void NRG_Wilson(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::gpara.init(); try {

   MX_CHECK_HELPER_NARGS(1,-1,-1);

   char cmplx=0;
   unsigned r; int t;

   if (nargin<=2) {
      const char *var[] = { "H0","A0","FC","Z","FL" };
      OPTS opts; mxArray *a;

      if (!mxIsChar(argin[0])) wblog(FL,
         "ERR usage #2 (1st argument expected mat-file name for setup)");
      if (mxGetString(argin[0],str,128)) wblog(FL,
         "ERR failed to read arg #1 (string, usage #2)");

      opts.init(str); 
      for (int i=0; i<5; ++i) { opts.getOpt(var[i],a);
         if (a) { r=-1;
            t=mxIsQSpace(0,0,a,r,'c',-2);
            if (t>0 && t&4) { cmplx|=1; break; }
         }
      }
      if (!cmplx) { opts.getOpt("ff",a);
         if (a) { r=-1;
            t=Mx::IsNumArray(0,0,argin[3],r,'c');
            if (t>0 && t&4) { cmplx|=1; }
         }
      }
      opts.checkAnyLeft(); 
   }
   else {
      wbindex iq("1 2 5 6 8", 1); 

      for (unsigned i=0; i<iq.len && iq[i]<(unsigned)nargin; ++i) { r=-1;
         t=mxIsQSpace(0,0,argin[iq[i]],r,'c',-2);
         if (t>0 && t&4) { cmplx|=1; break; }
      }
      if (!cmplx) { r=-1;
         t=Mx::IsNumArray(0,0,argin[3],r,'c'); 
         if (t>0 && t&4) { cmplx|=1; }
      }
   }

   if (cmplx)
        { NRG_Wilson<gTQ,wbcomplex>(nargout,argout,nargin,argin); }
   else { NRG_Wilson<gTQ,double>   (nargout,argout,nargin,argin); }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in NRGWilsonQS"); }
}; 

template<class TQ, class TD>
void NRG_Wilson(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
) {
    unsigned i,j,k,l,r,m=0,n,N, iter=0, Nkeep=256;
    double dbl, Etrunc1, Etrunc=0, Estop=0, deps=1E-12, db=-1;
    int vflag=Wb::envVRB; 
    char nostore=0, toFile=0, wf=-1, checkIdA=1, cmplx=ISCOMPLX_(TD);
    char akflag=0, zflag=1, cgflag=0, fflag=0;
    unsigned gotNK=0, gotTR=0;
    int dmax=-1, NEE=-1;

    wbvector<double> E4,E0,ED,dd;
    wbMatrix<unsigned> D4;
    wbMatrix<double> gg,EE,EK, E4_;
    wbMatrix<TD> ff;

    wbvector<TD> f;
    wbMatrix<TQ> QS;
    wbMatrix<int> NK;

    int idAK,idAD,idHK,idHD, idE0, idES;

    wbvector< QSpace<TQ,TD> > FC, FG, F1, F2, F1K, HKALL; 
    QSpace<TQ,TD> AK, AD, H4,A4, FCC, FX, Z;
    QSpace<TQ,TD> H0, A0, HK_;
    QSpace<TQ,double> HK, HD, ID;

    const wbperm PA; 

    OPTS opts;
    wbstring fout;
    mxArray *S, *a;

    const char* vpass[] = { "param", "Gamma" };
    unsigned npass=sizeof(vpass)/sizeof(const char*);
    mxArray* apass[npass];

    time_t tstart=time(NULL);

    Wb::SigHandler SIG(FL);
    Wb::Clock nrgTime_1("build H4"); 
    Wb::Clock nrgTime_2("eig(H4)");
    Wb::Clock nrgTime_3("compute AK,AD");
    Wb::Clock nrgTime_4("update ops");

    gCPUTime.init();

    FNfac=-1; 

    Wb::Clock_resume sw0(&nrgTime_0); 

    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"start");
    #endif

#ifdef MATLAB_MEX_FILE
    mexAtExit(myCleanUp);  
#endif
    Wb::ResSummary(NULL,1);

    if (nargin<=2) {
       if (mxGetString(argin[0],str,128)) wblog(FL,
          "ERR Error reading string of arg #1");

       opts.init(str); 

       opts.getOpt("H0",a,'!'); H0.init(FL,a); 
       opts.getOpt("A0",a,'!'); A0.init(FL,a);
       opts.getOpt("Lambda", Lambda,'!');

       opts.getOpt("ff",a,'!'); ff.init(FL,a);
       opts.getOpt("FC",a,'!'); mxInitQSpaceVec(FL,a,FC);
       opts.getOpt("Z", a,'!'); if (!mxIsEmpty(a)) Z.init(FL,a);

       opts.getOpt("gg",a);  if (a) gg.init(a);
       opts.getOpt("FL",a);  if (a) mxInitQSpaceVec(FL,a,FG);

       for (i=0; i<npass; ++i) {
          opts.getOpt(vpass[i],a);
          apass[i] = a ? mxDuplicateArray(a) : NULL;
       }

       if (nargin>1) {
          if (!mxIsChar(argin[1]) || mxGetString(argin[1],str,128))
          wblog(FL,"ERR reading string of arg #2 (fout)");
          fout=str;
       }
    }
    else {
       m=6; if (nargin<int(m) || nargout>2)
          usage(FL,"invalid number of I/O arguments");

       wbindex iq("1 5 6", 1); 
       if (mxIsEmpty(argin[5])) iq.len=2; 

       for (i=0; i<iq.len; ++i) {
          try { r=-1; mxIsQSpace(FL,argin[iq[i]],r,cmplx,-2); }
          catch (...) {
             wblog(FL,"ERR invalid QSpace at arg #%d",iq[i]+1);
          }
       }

       H0.init(FL,argin[0]);

       i=1; if (!mxIsEmpty(argin[i])) { r=-1;
          if (mxIsQSpace(FL,argin[i],r)<=0) wblog(FL,
             "ERR invalid QSpace at arg #%d", i+1);
          A0.init(FL,argin[i]); 
       }

       if (mxGetNumber(argin[2], Lambda)) wberror(FL,str);
       ff.init(FL,argin[3]);

       mxInitQSpaceVec(FL,argin[4], FC);
       if (!mxIsEmpty(argin[5])) Z.init(FL,argin[5]);

       if (nargin>=int(m+2)) {
          if (mxIsEmpty(argin[m]) && mxIsEmpty(argin[m+1])) { m+=2; } else
          if (mxIsQSpaceVec(0,0,argin[m+1],-1,NULL,NULL,'c')) {
             if (!mxIsNumeric(argin[m]) || mxIsComplex(argin[m]))
                wblog(FL,"ERR %s() real gg data expected",myname);
             gg.init(argin[m]);
             mxInitQSpaceVec(FL,argin[m+1], FG); m+=2;
          }
       }

#ifdef MATLAB_MEX_FILE
       for (i=0; i<npass; ++i)
       apass[i]=mexGetVariable("caller",vpass[i]); 
#else
       npass=0;
#endif

       opts.init(argin+m, nargin-m);
    }

    if (vflag &15) {
    wblog(FL,"=== %s ====================================",myname); }

    if (opts.getOpt(   "-q")) vflag= 0; else
    if (opts.getOpt(FL,"-v")) vflag|=2; else
    if (opts.getOpt(FL,"-V")) vflag|=6; 
    fflag=opts.getOpt("-f"); 

    opts.getOpt(FL,"deps",deps);
    opts.getOpt(FL,"db",  db  );
    opts.getOpt(FL,"dmax",dmax);

    opts.getOpt(FL,"Nkeep",Nkeep);
    opts.getOpt(FL,"Etrunc", Etrunc ); Etrunc1=1.2*Etrunc;
    opts.getOpt(FL,"ET1",Etrunc1);     if (Etrunc1<0) Etrunc1=Etrunc;

    opts.getOpt(FL,"zflag",zflag);
    opts.getOpt(FL,"NEE", NEE);

    if (opts.getOpt(FL,"-Estop","using Estop=1E-4")) Estop=1E-4;
    else opts.getOpt(FL,"Estop",Estop);

    if (fout.isEmpty()) {
       opts.getOpt(FL,"fout", fout);
    }

    opts.getOpt(FL,"NKEEP",NK);  
    opts.getOpt(FL,"ETRUNC",ED); 

    FNfac=-1;
    opts.getOpt(FL,"FNfac", FNfac); 

    nostore=opts.getOpt(FL,"ionly");
    if (opts.getOpt(0,0,"-IdF")) checkIdA=0;

    opts.checkAnyLeft(); 

    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"info");
    #endif

#if defined(DBG_QSX_BUF)
    BFF.flush();
    BFF.blogf(FL,"--- %s() having BUF[%ld], XBUF[%ld]",
      FCT, gCS.BUF.size(), gXS.XBUF.size());
#endif

    if (!H0.isHConj() || H0.rank(FL)!=2) { if (str[0])
       wblog(FL,"ERR H0 must be Hermitian rank-2 object\n%s",str); else
       wblog(FL,"ERR H0 must be Hermitian rank-2 object");
    }

    if (H0.itags.isEmpty()) {
       H0.init_itags(FL,"nrg:HK",0);
    }

    i=H0.skipZeroOffDiag(1E-14,'b'); if (i) {
       wblog(FL," *  skipped %d off-diagonal blocks from H0",i);
      #ifdef MATLAB_MEX_FILE
       H0.put("H0_");
      #endif
    }

    AK=A0;
    if (AK.isEmpty()) { akflag=1;
       wblog(FL,"WRN got empty A0 (building from H0)");
       AK.initIdentity(H0);
       AK.PrependSingletons(3); 
       AK.Permute("1 3 2");     
    }
    else {
       wbvector<widx_t> D,D_;
       if (!AK.isConsistent_r(3)) wblog(FL,
          "ERR invalid %s (%s)", nargin>1 ? "arg #2" : "QSpace A0", str);

       AK.getDim(D,&D_); 
       if (D.len!=3 || D_.len!=3) wblog(FL,
          "ERR %s() got r=%g/%g/3 !?",myname,D.len,D_.len);
       n=D_[0]*D_[2]; if (D_[1]!=n) {
          wbstring ss(D_.toStrf("","x"));
          if (D_[1]<n) wblog(FL,
             "WRN got truncated A0 @ %s (LRs) !?",ss.data);
          else wblog(FL,"ERR %s() got dim(A0) = %s !?",myname,ss.data);
       }

       AK.contract(FL,"13*",AK,"13",FX);
       if (!FX.isIdentityMatrix(1E-12)) {
          MXPut(FL,"a").add(AK,"AK").add(FX,"FX");
          wblog(FL,"ERR A0 does not describe orthonormal basis");
       }

       if (H0.hasQOverlap(1,AK,2,'<')<=0) wblog(FL, 
          "ERR H0 must be in R-basis of A0 (assuming LRs)");

       FX.init2DiffOp(AK,1,H0,0); 
       if (!FX.isEmpty()) { m=FX.DATA.len;
          FX.Append2AndDestroy(FL,H0); sprintf(str," *  "
            "increased H0 by %d diagonal zero-block%s (%ld->%ld)",
             m, m!=1 ? "s":"", H0.DATA.len-m, H0.DATA.len);
          wblogBuf.push(FL,str);
       }
    }

    if (A0.itags.isEmpty()) {
       A0.init_itags(FL,"nrg:AK",0);
    }

    if (Lambda<1) wblog(FL,"ERR %s() invalid Lambda=%g !?",myname,Lambda);

    if (zflag<=1 && Lambda>1 && H0.qtype.len>1 &&
       H0.qtype.allNonAbelian()) wblog(FL,
       "WRN %s() having zflag=%d with sym=%s\n"
       "(hint: particle-hole symmetry requires zflag>=2)",
       myname, zflag, H0.qtype.toStr().data
    );

    if (!(gotNK=NK.numel())) { NK.init(1,1); NK[0]=Nkeep; }

    {  int nk0=-1; 
       H0.EigenSymmetric(AK,AD,HK,HD,E4_,D4,nk0);
       if (nk0>NK[0]) { wblog(FL,
          "WRN keeping full input H0 at D=%d / Nkeep=%d",nk0, NK[0]);
       }; NK[0]=nk0; 

       E4_.getCol(0,E4);

       HK_=H0;  
       AK=A0;
    }

    for (k=0; k<FC.len; ++k) {
       cgflag=FC[k].gotCGS(FL); r=FC[k].rank(FL);

       if (!cgflag && r==3) {
          r=FC[k].Reduce2AbelianOp(FL); 
       }

       if ((cgflag<=0 && r>2) || (cgflag>0 && r>3)) {
          l=0; str[0]=0;
          if (FC[k].qtype.len) l+=snprintf(
             str+l,99,"; q=%s",FC[k].qtype.toStr().data);
          if (FC[k].otype!=QS_NONE) l+=snprintf(
             str+l,99,"; t=%s",FC[k].otype2Str().data);
          wblog(FL,
         "ERR %s() got rank-%d operator (%d%s)",myname,r,cgflag,str);
       }

       if ((cgflag>0 || FC[k].itags) && !FC[k].itags.isOp()) wblog(FL,
          "WRN FC[%d] with operator itags=%s",k+1,STR(FC[k].itags));
       FC[k].init_tags(); 
    }

    if (zflag) {
       if (Z.isEmpty()) wblog(FL,
          "ERR %s() got empty Z-operator for fermionic signs",myname);
       if (!Z.isDiagMatrix()) wblog(FL,
          "ERR %s() got not diagonal operator Z",myname);
    }

    if (zflag<0 || zflag>3) wblog(FL,
       "ERR %s() invalid zflag=%d",myname,zflag);
    m=FC.len; F1.init(m); F2.init(m);

    for (i=0; i<m; ++i) {
       if (!FC[i].isOperator()) wblog(FL, 
          "ERR invalid operator FC(%d)\n"
          "(hint: set otype='operator' for rank-3 operators)",i+1);
       if (zflag<=1)
            FC[i].hconj(F1[i]);              
       else F1[i]=FC[i]; 

       if (zflag)
          { Z.contract(FL,2,FC[i],1, F2[i]); 
            F2[i].otype=FC[i].otype; }
       else F2[i]=FC[i];
    }

    if (1) { QSpace<TQ,TD> FY;

       if (FC.len && FC[0].rank(FL)==2) {
          for (i=0; i<F1.len; ++i) {
             F1[i].TensorProd(F2[i],FX).hconj(FY);
             FX.Append2AndDestroy(FL,FCC);
             FY.Append2AndDestroy(FL,FCC);
          }
       }
       else if (FC.len && FC[0].otype==QS_OPERATOR) {
          wbperm P("1,3,2,4");

          for (i=0; i<F1.len; ++i) {
             if (zflag<=1) {
                F1[i].contract(FL,"3", F2[i],"3",FX,P);
                FX.hconj(FY);
             }
             else {
                F1[i].contract(FL,"3*",F1[i],"3",FX,P);
                F2[i].contract(FL,"3*",F2[i],"3",FY,P);
             }
             FX.Append2AndDestroy(FL,FCC);
             FY.Append2AndDestroy(FL,FCC);
          }
       }
       else wblog(FL,"ERR invalid FC operator set");

       FCC.MakeUnique(); FCC.SkipZeroData(1E-14);

       if (!FCC.isHConj()) {
          MXPut(FL,"I4").add(F1,"F1").add(F2,"F2").add(FCC,"FCC")
             .add(zflag,"zflag");
          wblog(FL,"ERR %s() coupling "
             "F'*F does not yield hermitian H\n'%s'",myname,str
          );
       }
    }

    if (ff.isEmpty() || FC.isEmpty()) wblog(FL,
       "ERR %s() got empty couplings (%dx%d,%d)",
        myname,ff.dim1,ff.dim2,FC.len);

    if (F1.len>1 && ff.isVector()) {
       unsigned n=ff.numel(), m=F1.len;
       wbMatrix<TD> fx; fx.swap(ff); ff.init(n,m);
       for (i=0; i<n; ++i)
       for (j=0; j<m; ++j) ff(i,j)=fx.data[i];
    }
    else if (ff.isVector() && ff.dim1==1) {
       SWAP(ff.dim1,ff.dim2); 
    }

    if (!gg.isEmpty() && gg.isVector()) {
       if (FG.len>1) {
          unsigned n=gg.numel(), m=FG.len;
          wbMatrix<double> gx; gx.swap(gg); gg.init(n,m);
          for (i=0; i<n; ++i)
          for (j=0; j<m; ++j) gg(i,j)=gx.data[i];
       }
       else if (gg.isVector() && gg.dim1==1) {
          SWAP(gg.dim1,gg.dim2); 
       }
    }

    N=1+ff.dim1; 

    if (!ff.isNormal()) wblog(FL,
       "ERR %s() got NaN in nearest-neighbor couplings ff !?",myname);
    if (!gg.isNormal()) wblog(FL,
       "ERR %s() got NaN in local energies gg !?",myname);

    if (ff.dim2!=F1.len && ff.dim2!=F1.len*F1.len) wblog(FL,
       "ERR size mismatch ff (%d x %d / %d FC ops)",ff.dim1, ff.dim2,F1.len);
    if (!gg.isEmpty() && gg.dim1<ff.dim1) wblog(FL,
       "ERR size mismatch gg (%d / %d ff.len) x %d",gg.dim1,ff.dim1,gg.dim2);
    if (gg.dim2!=FG.len) wblog(FL, 
       "ERR size mismatch gg (%d x %d / %d FL ops)",gg.dim1,gg.dim2, FG.len);

    if (F1.len!=F2.len) wblog(FL,"ERR %d/%d",F1.len,F2.len);
    if (!Z.isEmpty() && (Z.QDIM!=H0.QDIM || Z.rank(FL)!=2)) wblog(FL,
       "ERR severe rank mismatch of Z (%dx%d/%d)",
        Z.QIDX.dim1, Z.QIDX.dim2, Z.QDIM
    );

    for (k=0; k<F1.len; ++k) {
       if (F1[k].SkipZeroData()) wblog(FL,
          "WRN %s() FC[%d] operator contains zero-blocks",myname,k+1);
       if (F1[k].isEmpty()) wblog(FL,"ERR got empty FC[%d]",k+1);

       if (F1[k].QDIM!=H0.QDIM || !F1[k].isOperator()) wblog(FL,
          "ERR %s() symmetry or rank mismatch of FC[%d] "
          "(qdim=%d/%d; r=%d/%s/%d)",myname,k+1,F1[k].QDIM, H0.QDIM,
          F1[k].rank(FL), QS_STR[F1[k].otype], F1[k].isOperator());

       if (F2[k].QDIM!=H0.QDIM || !F2[k].isOperator()) wblog(FL,
          "ERR symmetry or rank mismatch of FC[%g] "
          "(q=%d,%d; r=%d/%s/%d)",k+1,F2[k].QDIM, H0.QDIM,
          F2[k].rank(FL), QS_STR[F2[k].otype], F2[k].isOperator());

       F1[k].checkQ(FL,F2[k]);
       if (F1[k].QIDX.dim2!=F2[k].QIDX.dim2 ||
          F1[k].otype!=F2[k].otype) wblog(FL,
          "ERR FC[%g] operator rank mismatch (qdim=%d/%d; %d/%d; %s/%s)",
          k+1, F1[k].QIDX.dim2, F2[k].QIDX.dim2, F1[k].QDIM, F2[k].QDIM,
          QS_STR[F1[k].otype], QS_STR[F2[k].otype]);

       if (k) {
          char ic=F1[k].gotCGS(FL); if (cgflag!=ic) wblog(FL,
          "ERR %s() CG inconsistency (%d/%d)",myname,ic,cgflag);
       }
       else { cgflag=F1[k].gotCGS(FL); }
    }

    for (k=0; k<FG.len; ++k) {
       if (FG[k].SkipZeroData()) { wblog(FL,
          "WRN local operator FL(%d) contains zero data - skip.", k+1);
           continue;
       }
       if (FG[k].isEmpty()) { wblog(FL,
          "WRN no contribution by local operator FL(%d) !?", k+1);
           for (i=0; i<gg.dim1; ++i) gg(i,k)=0;
           continue;
       }
       if (FG[k].QDIM!=H0.QDIM || FG[k].rank(FL)!=2) { wblog(FL,
          "ERR severe local operator inconsistency [FL(%d): %d,%d]",
           k+1, FG[k].QDIM, FG[k].rank(FL));
       }
       if (!FG[k].isHConj()) wblog(FL,
       "ERR FL(%d) not a Hermitian rank-2 object", k+1);

       if (k) {
          char ic=FG[k].gotCGS(FL); if (cgflag!=ic)
          wblog(FL,"ERR %s() CG inconsistency (%d/%d)",myname,ic,cgflag);
       }
       else cgflag=FG[k].gotCGS();

       if (FG[k].rank(FL)!=2) wblog(FL,
          "ERR %s() got non-scalar rank-%d local operator",myname,
          FG[k].rank(FL)
       );

       if (!FG[k].itags.isOp()) wblog(FL,
          "WRN FG[%d] with operator itags=%s !?",k+1,STR(FG[k].itags));
       FG[k].init_tags(); 
    }

    if (Nkeep>9999 && !fout && !fflag) {
       fout="./NRG/NRG"; wblog(FL,
       "WRN saving data to file for Nkeep=%d (/9999)", Nkeep);
    }
    toFile = !fout.isEmpty();

    if (!Nkeep) {
       wblog(FL,"<i> %NNkeep=%d life is simple - exit B)", Nkeep);
       return;
    }

    if (Estop<0) wblog(FL,"ERR invalid Estop=%g",Estop);
    if (Estop>0.1) wblog(FL,"ERR invalid Estop=%g (expected <0.1)",Estop);

    if (FNfac<0) {
       FNfac=1; 
       unsigned k=(ff.dim1>=4 ? ff.dim1-4 : 0); double q, qmin=-1;
       for (; k<ff.dim1; ++k) {
          q=ff(k,0)/nrgScale(k+1); 
          qmin = (qmin>0 ? MIN(qmin,q) : q);
       }
       FNfac=qmin;
       if (FNfac<0.1 || FNfac>10) {
       wblog(FL,"WRN %s() got FNfac=%.3g!?",myname,FNfac); }
    }

    if (vflag) {
       wblog(FL,"--- %48R","-");
       l=sprintf(str,"Lambda=%g, L=%d",Lambda,N);
       if (Etrunc>0)
            { l+=sprintf(str+l,", Etrunc=%.3g (@%d)",Etrunc,Nkeep); }
       else { l+=sprintf(str+l,", Nkeep=%d",Nkeep); }
       if (!toFile)
            { l+=sprintf(str+l,", internal"); }
       wblog(FL," *  %s\nsym=%s",str,STR2(A0.qtype,'V'));

       str[0]=0; l=0;
       if (vflag &14) l+=sprintf(str+l," vflag=%d",vflag); 
       if (nostore) l+=sprintf(str+l," %s",toFile? "NOSTORE":"noStore");
       if (l) wblog(FL," *  flags:%s",str);

       if (WbUtil<TD>::isComplex()) {
           wblog(FL," *  running in complex mode");
       }
       wblog(FL,"=== ================================================");

       if (!gg.isEmpty())
          wblog(FL," *  got %d local operator%s",FG.len,FG.len!=1 ? "s":"");
       wblogBuf.flush();
    }

    if (toFile) {
       wbstring cwd(128); if (getcwd(cwd.data, 127)==NULL) wblog(FL,
         "ERR cwd length exceeds maximum length %d\n%s",cwd.len,cwd.data);

       if (nostore) wblog(FL,"NB! saving info only");
       else {
          for (iter=0; iter<1024; ++iter) {
              sprintf(str, "%s_%02d.mat", fout.data, iter);
              if (remove(str)) break;
          }   if (iter>1000) wblog(FL,"WRN iter=%d !?", iter);
       }

       char *s1=str, *s2=str+32; s2[0]=0;

       if (nostore)
            strcpy (s1,"info");
       else sprintf(s1,"{0-%d}",N);

       if (iter) sprintf(s2," (%d files removed)",iter);

       if (vflag) printf("\n");
       printf("   pwd: %s\n"
              "   out: %s_%s.mat%s\n",
       Wb::repHome(cwd).data, Wb::repHome(fout).data, s1, s2);
       if (vflag) printf("\n");
    }
    else if (fout.data && fout.data[0]) wblog(FL,
       "WRN keeping data internally\nyet got fout='%s' !?",fout.data);
    else wblog(FL,"<i> internal mode (no file I/O)");

    { wbvector< wbMatrix<TQ> > qq(F1.len+FG.len);
      wbvector<widx_t> sz;

      for (k=i=0; i<F1.len; ++i) F1[i].getQDim(qq[k++],sz);
      for (  i=0; i<FG.len; ++i) FG[i].getQDim(qq[k++],sz);

      QS.CAT(1,qq).makeUnique();

      if (checkIdA) { wbMatrix<TQ> Qs;
         AK.getQsub(2,Qs); 
         Qs.makeUnique(); if (Qs!=QS) {
#ifdef MATLAB_MEX_FILE
         Qs.Print("A0->QS"); QS.Print("FX->QS");
#endif
         if (akflag) wblog(FL,
            "ERR local QIDX inconsistency%N%N    Hint: "
            "A0 not specified and failed to construct default from H0");
         else wblog(FL,
            "ERR local QIDX inconsistency in specified A0%N%N    Hint: "
            "ensure that H0 contains complete basis even if data is zero"
         );
      }}
    }

    ID.QIDX.Cat(2,QS,QS); ID.QDIM=QS.dim2; ID.setupDATA();

    for (k=0; k<QS.dim1; ++k) {
       TQ *qk=ID.QIDX.rec(k);

       for (i=0; i<F1.len; ++i) { if (F1[i].findDimQ(qk,n)) break; }
       if (i==F1.len) {
          for (i=0; i<FG.len; ++i) { if (FG[i].findDimQ(qk,n)) break; }
          if (i==FG.len) {
             MXPut(FL).add(ID,"ID").add(F1,"F1").add(FG,"FG");
             wblog(FL,"ERR failed to find Id.Q(%d,:) in F[1L]",k+1);
          }
       }
       ID.DATA[k]->initIdentity(n);
    }

    ID.qtype=AK.qtype;

    if (cgflag>0 || !H0.CGR.isEmpty())
         ID.initIdentityCGS(); 
    else ID.itags.init_qdir("+-");

    if (checkIdA) {
       QSpace<TQ,TD> E2;

       AK.contract(FL,"12*",AK,"12",E2); 
       if (!E2.DATA.len) wblog(FL,"ERR AK/AK contracted to empty!?");
       E2*=(1./E2.DATA[0]->data[0]); 

       try { dbl=E2.normDiff2(ID); }
       catch (...) {
          MXPut(FL).add(E2,"E2").add(ID,"ID");
          throw;
       }
       if (dbl>1E-12) { 
          wbvector< QSpace<TQ,TD> > AA(1);
          double db2=0; wbindex Ia(1); Ia[0]=2;

          AA[0].init2ref(AK);
          E2.initIdentityCG(AA,Ia);

          try { db2=E2.normDiff2(ID); }
          catch (...) {
             MXPut(FL).add(E2,"E2").add(ID,"ID");
             throw;
          }

          if (db2<1E-12)
             wblog(FL,"WRN got truncated A0 (e=%.4g,%.4g) !?",dbl,db2);
          else {
             MXPut(FL,"i").add(ID,"ID").add(E2,"id");
             wblog(FL,"ERR local state space (e=%.4g/%.4g) !?",dbl,db2);
          }
       }
    }

    dbl=1/nrgScale(0); 

    if (dbl!=1) wblog(FL,
       "ERR %s() got nrgScale(0)=%g (expecting 1.)",myname,dbl);

    dbl=pow(Lambda,-N/2.);
    if (ff.dim1 && (fabs(ff.recMax(ff.dim1-1) / ff.recMax(0)) > 10*dbl))
    wblog(FL,
       "WRN ff must fall off exponentially\n[%.3g .. %.3g; %.3g] !?",
        ff.recMax(N-1), ff.recMax(0), dbl);

    for (i=0; i<ff.dim1; ++i) { dbl=1/nrgScale(i+1);
       for (j=0; j<ff.dim2; ++j) ff(i,j)*=dbl;
       for (j=0; j<gg.dim2; ++j) gg(i,j)*=dbl;
    }

    { const char* fld[] = { "AK","AD","HK","HD","E0","ES" };
      S=mxCreateStructMatrix(1, toFile || nostore? 1:N, 6, fld);
      idAK=0; idAD=1; idHK=2; idHD=3; idE0=4; idES=5;

      if (!S) wblog(FL,"ERR %s() failed to intialize mex !?",myname);
    }

    if (!NK.isVector()) wblog(FL,"ERR invalid NKEEP");
    else {
       wbvector<int> nk(NK.numel(),NK.data);
       NK.init(N, cgflag>0 ? 4:2);
       for (n=MIN(unsigned(nk.len),N), i=0; i<n; ++i) {
          NK.setRec(i,nk[i]);
       }
       for (; i<N; ++i) NK.setRec(i,Nkeep);
       NK(N-1,0)=0; 
    }

    if (NEE<0) NEE=2*Nkeep;

    EE.init(N, (unsigned)NEE); EE.set(NAN);
    E0.init(N); EK.init(N,3); EK.set(NAN);
    HKALL.init(N); 

    #ifdef __WB_MEM_CHECK__
       Wb::MemCheck(FL,"info");
    #endif

    for (iter=0, NRG_N=N; iter<N; ++iter) { NRG_ITER=iter; SIG.check911();

        Wb::Clock_resume sw(&nrgTime_1); 

       if (iter) {
          HK_ *= (nrgScale(iter-1) / nrgScale(iter));

          if (cgflag<=0) { try {
             f.init2ref(ff.dim2, ff.rec(iter-1));
             nrgBuildH4_abelian(H4, HK_, ID,
                f, F1K, F2, 
                gg.dim1 ? gg.rec(iter-1) : NULL, FG);
             }
             catch (...) {
                wblog(FL,"ERR %s() iter=%d",myname,iter);
             }

          }
          else {
             f.init2ref(ff.dim2, ff.rec(iter-1));
             nrgBuildH4_cg(H4,A4, HK_, f, F1K,
                (zflag<=1 || !wf) ? F2 : F1, 
                 zflag<=1 ? 1 : 0, 
                gg.dim1 ? gg.rec(iter-1) : NULL, FG);

             if (iter==1 && !H4.isHConj(0,0,1E-12,'v')) { 
                MXPut(FL,"a").add(HK,"HK").add(ff.getRec(iter-1),"ff")
                .add(gg.dim1? gg.getRec(iter-1) : wbvector<double>(),"gg")
                .add(F1K,"F1").add((zflag<=1 || !wf) ? F2 : F1,"F2")
                .add(FG,"FG");
                wblog(FL,"ERR H4 got non-hermitian operator setting !?\n"
                "hint: got correct zflag=%d ?\n%s",zflag,str);
             }
          }

          sw.Switch(&nrgTime_2); 

          #ifdef __WB_MEM_CHECK__
             Wb::MemCheck(FL,"info");
          #endif

          H4.EigenSymmetric(
             AK, AD, HK, HD, E4_, D4, NK(iter,0),
             iter<ED.len && ED[iter]>0 ?  ED[iter]
             : ((Etrunc1>0 && !gotTR) ? Etrunc1 : Etrunc),
             E0.data+iter, PA, deps, db, dmax 
          ); E4_.getCol(0,E4);

          if (!gotTR) {
             if (NK(iter,0)<int(E4.len)) gotTR=iter+1; else {
                double x=1.1*(iter+1>=ED.len ? Etrunc : ED[iter+1]);
                if (E4.last()>x && (ED.len || Etrunc1>x)) { gotTR=iter+2; }
             }
          }

          #ifdef __WB_MEM_CHECK__
             Wb::MemCheck(FL,"info");
          #endif

          sw.Switch(&nrgTime_3); 

          if (cgflag>0) { widx_t nx;
             n=AK.getDim(AK.rank(FL)-1, &nx); NK(iter,1)=nx;
             NK(iter,2)=D4.colSum(0);
             NK(iter,3)=D4.colSum(1);

             A4.contract(2,AK,1,AK,"132"); if (!AD.isEmpty())
             A4.contract(2,AD,1,AD,"132"); 
          }
          else {
             AK.Permute("2,3,1"); if (!AD.isEmpty())
             AD.Permute("2,3,1"); 
             NK(iter,1)=D4.colSum(0);
          }
       }

       AK.init_itags(FL,"nrg:AK",iter); 
       AD.init_itags(FL,"nrg:AD",iter);
       HK.init_itags(FL,"nrg:HK",iter);
       HD.init_itags(FL,"nrg:HD",iter);

       if (iter) { HKALL[iter]=HK; }
       else {
          HKALL[iter]=HK_.getReal(); 
       }

       if (iter) { HK_.initT(HK); }

       wf=(char(iter%2)==(zflag-2));

       sw.Switch(&nrgTime_4); 

       updateFOps(F1K, AK,AK,
          (zflag<=1 || wf) ? F1 : F2, 
          iter+1>=N
       );

       sw.stop();

       l=MIN( size_t(NK(iter,0)), EE.dim2 ); 
       memcpy(EE.rec(iter), E4.data, l*sizeof(double));

       if (iter>18 && Estop>0) {
          unsigned nk=MIN(MIN(int(EE.dim2),NK(iter-2,0))/4,NK(iter,0)/4);
          double x=sqrt(Wb::rangeNormDiff2(EE.rec(iter-2), EE.rec(iter), nk));
          if (x<Estop && iter+2<N) {
             wblog(FL,"%N==> "
              "NB! energy flow converged at %.3g (N=%d->%d)",x,N,iter+2);
             N=iter+2;
             Estop=-Estop; 

             EE.dim1=N; EK.dim1=N; NK.dim1=N;
             E0.len=N; HKALL.len=N;
             ff.dim1=N; if (gg.dim1>N) { gg.dim1=N; }
          }
       }

       n=NK(iter,0);
       EK(iter,0)=(n && n<E4.len ? E4[n-1] : E4.last()); 
       EK(iter,1)=(     n<E4.len ? E4[n  ] : NAN);       
       EK(iter,2)=E4.last();

       dbl = !ISNAN(EK(iter,1)) ? EK(iter,1) :  EK(iter,0);
       nrgDispIter(iter,N,E4, dbl, HK, vflag ?
         (vflag &14 ? vflag : (Etrunc<=0 ? 0 : (dbl>=0.8*Etrunc ? -1 : +1)))
          : -1
       );

       if (!nostore) { k = toFile ? 0 : iter;
          char qcl=(toFile ? '!':0); 

          mxReplaceField(FL, S, k, idAK, AK.toMx(qcl));
          mxReplaceField(FL, S, k, idAD, AD.toMx(qcl));
          mxReplaceField(FL, S, k, idHK, iter? HK.toMx(qcl) : H0.toMx(qcl));
          mxReplaceField(FL, S, k, idHD, HD.toMx(qcl));

          mxReplaceField(FL, S, k, idE0, numtoMx(E0[iter]));
          mxReplaceField(FL, S, k, idES, numtoMx(nrgScale(iter)));

          if (toFile) {
             Wb::UseClock NT5(&nrgTime_5);

             sprintf(str, "%s_%02d.mat", fout.data, iter);
             if (vflag &4) wblog(FL,"\nI/O writing %s",Wb::basename(str));

             Wb::matFile F(FL,str,"w");
             for (n=mxGetNumberOfFields(S), i=0; i<n; ++i) {
                F.put(FL,
                   mxGetFieldNameByNumber(S,i),
                   mxGetFieldByNumber(S,0,i)
                );

             }
             F.close();
          }
       }
       doflush(); 

       #ifdef __WB_MEM_CHECK__
          Wb::MemCheck(FL,"info");
          wblog(FL,"--- %50R","-");
       #endif
    }

    if (vflag) printf("\n");

    m=NK.colMax(0,i);
    n=NK.colMax(NK.dim2>3? 2:1,j);
    l=sprintf(str,"NK=%s / %s",I2STR(m),I2STR(n));
    if (NK.dim2>3) {
       l+=sprintf(str+l," (%s / %s)",I2STR(NK(i,1)),I2STR(NK(j,3))); }

    for (dbl=1E99, i=0; i<EK.dim1; ++i) {
       if (!ISNAN(EK(i,1)) && EK(i,1)>EK(i,0)) {
       if (dbl>EK(i,1)) dbl=EK(i,1); }
    }

    l+=sprintf(str+l," @ Etr=%.4g",dbl);
    if (Etrunc>0) l+=sprintf(str+l," / %g",Etrunc);

    if (dbl>0.9*Etrunc)
         wblog(FL,"==> %s",str);
    else wblog(FL,"WRN %s",str);

    if (!toFile || (toFile && nargout>1))
    argout[0] = S; 

    sprintf(str,"NRG data obtained using %s",myname);
    wbstring ver(str);

    MXPut Iout(0,0,"Inrg"); sw0.stop();
    Iout.add(ver,"istr").add(Wb::TimeStamp(),"stamp");
    Iout.addP(MXPut(0,0)
         .add (wbstring().time_sys(tstart),"started")
         .add (wbstring().time_sys(),"finished")
         .addP(nrgTime_0.toMx(),"time_all")
         .addP(nrgTime_1.toMx(),"time_1")
         .addP(nrgTime_2.toMx(),"time_2")
         .addP(nrgTime_3.toMx(),"time_3")
         .addP(nrgTime_4.toMx(),"time_4")
         .addP(nrgTime_5.toMx(),"time_5")
     .toMx(),"usage");

#ifdef __WB_MEM_CHECK__
    Iout.addP(Wb::gML.totStr('l'),"MEM");
#endif

    MXPut Iops(0,0);
      Iops.add(ID,"ID").add(H0,"H0").add(A0,"A0").add(FC,"FC").add(Z,"Z")
        .add(ff,"ff").add(F1,"F1").add(F2,"F2").add(gg,"gg");
      if (!FG.isEmpty()) Iops.add(FG,"FL"); else Iops.add(dd.init(),"FL");
      Iops.add(FCC,"FCC").add(zflag,"zflag");

    dd.init(N); for (i=0; i<N; ++i) dd[i]=nrgScale(i);

    Iout.addP(EE.toMx('r'),"EE").add(HKALL,"HK").add(E0,"E0").add(EK,"EK")
     .add(dd,"EScale").add(getPhysE0(E0),"phE0")
     .addP(Iops.toMx(),"ops");

    dd.init(4);
      dd[0]=CG_EPS1; dd[2]=CG_SKIP_EPS1;
      dd[1]=CG_EPS2; dd[3]=CG_SKIP_EPS2;
    Iout.add(Lambda,"Lambda")
     .addP(MXPut(0,0)
        .add(dd,"cg_eps").add(deps,"deps").add(db,"db")
        .add(dmax,"dmax").add(-Estop,"Estop").add(vflag,"vflag")
        .add(toFile,"toFile").addP(STR(A0.qtype),"sym")
        .add(WbUtil<TD>::isComplex(),"complex") 
     .toMx(),"paras");

    for (i=0; i<npass; ++i) if (apass[i]) {
       Iout.addP(apass[i],vpass[i]); apass[i]=NULL;
    }

    Iout.addP( MXPut(0,0)
       .add(Etrunc,"Etrunc").add(Etrunc1,"Etrunc1")
       .add(gotTR,"itrunc").add(ED,"ETRUNC").add(FNfac,"FNfac")
     .toMx(),"Itr")
     .add(Nkeep,"Nkeep").add(NK,"NK").add(N,"N");

    Iout.save2(S);

    if (toFile) {
       i=0; n=mxGetNumberOfFields(S);
       sprintf(str,"%s_info.mat",fout.data);

       if (vflag &14) wblog(FL,
          "I/O saving info data to `%s'",Wb::basename(str));
       Wb::matFile F(FL,str,"w");

       for (; i<n; ++i) { matPutVariable(F.mfp,
          mxGetFieldNameByNumber(S,i), mxGetFieldByNumber(S,0,i));
       }

    }

    if (nargout>=2)        { argout[1]=S; } else
    if (toFile && nargout) { argout[0]=S; } else { mxDestroyArray(S); }

    if (vflag) { wblog(FL,
       "NRG I/O time usage: %s",SEC2STR(nrgTime_5.gettime()));
       myCleanUp(); 
       printf("\n");
    }

    if (!(vflag &12) && nrgTime_0.gettime()<3600 && CG_VERBOSE<6) {
       nrgTime_1.reset(); nrgTime_2.reset();
       nrgTime_3.reset(); nrgTime_4.reset();
    }
    else {
       nrgTime_0.info(); nrgTime_0.reset();
       nrgTime_5.info(); nrgTime_5.reset();
    }

#ifdef __WB_MEM_CHECK__

    wblog(FL,"%NMTR %s() check memory management (before cleanup)",myname);

    Wb::MemCheck(FL,"info");

    opts.init();

    ff.init(); gg.init(); dd.init();
    E4.init(); EE.init(); E0.init(); EK.init();
    D4.init(); QS.init(); NK.init(); fout.init();

    FC.init(); FG.init(); FCC.init(); FX.init();
    F1.init(); F2.init(); F1K.init(); Z.init();

    AK.init(); AD.init(); A0.init(); A4.init(); ID.init();
    HK.init(); HD.init(); H0.init(); H4.init(); HKALL.init();

    gCS.BUF.clear(); gCS.buf3.clear();
    gRS.buf.clear(); gXS.clear();

    Wb::MemCheck(FL,"info");

    wblog(FL,"MTR %s() check memory management (after cleanup)%N",myname);

#else
#endif

};

template <class TQ, class TD>
void nrgBuildH4_abelian(
    QSpace<TQ,TD> &H4,
    const QSpace<TQ,TD> &HK,
    const QSpace<TQ,double> Es, 
    const wbvector<TD> &f,
    const wbvector< QSpace<TQ,TD> > &F1K,
    const wbvector< QSpace<TQ,TD> > &F2,
    const double *g,
    const wbvector< QSpace<TQ,TD> > &FG
){
    unsigned i,j,l=0; char xflag=0, dflag=1;

    QSpace<TQ,TD> HX,HY, EL, FX;
    wbMatrix<TQ> QQ;

    EL.initIdentity(HK, NRG_ITER>1 ? 'r' : 0);

    if (!H4.isHConj()) wblog(FL,"XXX %s() H4 not h.conj. !?",myname);

    Es.TensorProd(HK,HX).save2(H4); 

    if (HK.QIDX.dim1==1 || FG.len) {
       H4.ExpandDiagonal(2,4); xflag=1;
    }

    if (!xflag) { H4.ExpandDiagonal(2,4); xflag=1; } 
    if (!H4.isHConj()) {
       MXPut(FL,"q0").add(H4,"H4").add(EL,"EL").add(HK,"HK").add(Es,"Es");
       wblog(FL,"ERR %s() H4 not h.conj. !?",myname);
    }

    for (i=0; i<FG.len; ++i) { if (g[i]==0) continue;
       FX=FG[i]; FX*=g[i];
       FX.TensorProd(EL,HX).Append2AndDestroy(FL,H4);
    }  

    if (f.len!=F1K.len) { dflag=0;
       if (f.len!=F1K.len*F2.len) wblog(FL,
       "ERR %s() length mismatch (%d/%d/%d)",myname,f.len,F1K.len,F2.len);
    }

    for (j=0; j<F2.len; ++j) {
    for (i=0; i<F1K.len; ++i, ++l) {
       if (dflag) { if (i!=j) continue; else l=i; }
       if (!f[l]) { continue; }

       FX=F2[j]; FX*=f[l]; FX.TensorProd(F1K[i],HX);
       HX.hconj(HY); 
       HX.Append2AndDestroy(FL,H4);
       HY.Append2AndDestroy(FL,H4);                
    }}
    H4.MakeUnique();

    if (!xflag)
    H4.ExpandDiagonal(2,4); 

    if (!H4.isHConj()) { 
       MXPut(FL,"q").add(H4,"H4").add(EL,"EL").add(HX,"HX").add(Es,"Es");
       wblog(FL,"ERR %s() H4 not h.conj. !?",myname);
    }
};

template <class TQ, class TD>
void nrgBuildH4_cg(
   QSpace<TQ,TD> &H4,
   QSpace<TQ,TD> &A4,
   const QSpace<TQ,TD> &HK,
   const wbvector<TD> &f,
   const wbvector< QSpace<TQ,TD> > &F1K,
   const wbvector< QSpace<TQ,TD> > &F2,
   char addHC,
   const double *g,
   const wbvector< QSpace<TQ,TD> > &FG
){
   unsigned i,j,l=0; char dflag=1;
   QSpace<TQ,TD> Q,FX,HX,AX;

   for (i=0; i<F1K.len; ++i) if (F1K[i].isEmpty()) wblog(FL,
       "WRN got empty space F1K[%d]",i+1);
   for (i=0; i<F2.len; ++i) if (F2[i].isEmpty()) wblog(FL,
       "WRN got empty space F2[%d]",i+1);
   for (i=0; i<FG.len; ++i) if (FG[i].isEmpty()) wblog(FL,
       "WRN got empty space FG[%d]",i+1);

   A4.initIdentityCG(F1K,F2,"1 3 2"); 

   Q=HK; Q.ExpandDiagonal();
   Q.contract(2,A4,1,AX);

   #ifdef __WB_MEM_CHECK__
      Wb::gML.printSize(FL,'l'); 
   #endif

   #ifdef LOAD_CGC_QSPACE
   #endif

   A4.contract("13*",AX,"13",H4);

   #ifdef LOAD_CGC_QSPACE
   #endif

   #ifdef __WB_MEM_CHECK__
      Wb::gML.printSize(FL,'l');
   #endif

   if (f.len!=F1K.len) { dflag=0;
      if (f.len!=F1K.len*F2.len) wblog(FL,
      "ERR %s() length mismatch (%d/%d/%d)",myname,f.len,F1K.len,F2.len);
   }

   for (j=0; j<F2.len; ++j) {
   for (i=0; i<F1K.len; ++i, ++l) {
      if (dflag) { if (i!=j) continue; else l=i; }
      if (!f[l]) { continue; }

      F2[j].times(f[l],Q); A4.contract(FL,3,Q,2,FX);

      if (F2[j].rank(FL)==2) {
         if (!addHC) wblog(FL,"WRN got addHC=%d !?",addHC);
         F1K[i].contract(FL,2,FX,1,AX);
      }
      else {
         if (addHC)
              F1K[i].contract(FL,"23", FX,"14",AX);
         else F1K[i].contract(FL,"13*",FX,"14",AX);
      }

      #ifdef __WB_MEM_CHECK__
         Wb::gML.printSize(FL,'l');
      #endif

      #ifdef LOAD_CGC_QSPACE
      #endif

      try {
         A4.contract(FL,"13*",AX,"13",HX); }
      catch (...) {
         MXPut(FL,"i4").add(A4,"A4").add(AX,"AX").add(HX,"HX")
          .add(FX,"FX").add(F1K[i],"F1K").add(i+1,"i");
         wblog(FL,"ERR %s()",FCT);
      }

      #ifdef __WB_MEM_CHECK__
         Wb::gML.printSize(FL,'l');
      #endif

      #ifdef LOAD_CGC_QSPACE
      #endif

      if (addHC) { QSpace<TQ,TD> HY;
      HX.hconj(HY).Append2AndDestroy(FL,H4); }
      HX.Append2AndDestroy(FL,H4);
   }}

   AX.init();

   for (i=0; i<FG.len; ++i) { if (g[i]==0) continue;
      FG[i].times(g[i],FX);
      A4.contract(FL,3,FX,2,Q);
      A4.contract(FL,"13*",Q,"13",HX);
      HX.Append2AndDestroy(FL,H4);
   }

   H4.MakeUnique();

   HX.init2DiffOp(A4,1,H4,0); 
   HX.Append2AndDestroy(FL,H4);

#ifdef LOAD_CGC_QSPACE
#endif

};

template <class TQ, class TD>
void updateFOps(
    wbvector< QSpace<TQ,TD> > &F12, 
    const QSpace<TQ,TD> &A1,
    const QSpace<TQ,TD> &A2,
    const wbvector< QSpace<TQ,TD> > &FC, 
    char lflag
){
    unsigned i; 
    QSpace<TQ,TD> Xk;

    if (A1.QDIM!=A2.QDIM) wblog(FL,
    "ERR %s() severe data inconsistency (%d,%d,%d)",FCT,A1.QDIM,A2.QDIM);

    F12.init(FC.len);

    for (i=0; i<FC.len; ++i) {

       A2.contract(FL,3,FC[i],2,Xk); 
       A1.contract(FL,"13*",Xk,"13",F12[i]); 
       F12[i].otype=FC[i].otype;

       if (!F12[i] && !lflag) wblog(FL,
          "WRN got empty space F12[%d] (%s)",i+1,IT2STR(A2));
    }

};

template <class TQ, class TD>
void nrgDispIter(
    const unsigned iter,
    const unsigned N,
    const wbvector<TD> &E4, const TD &EK,
    const QSpace<TQ,TD> &HK,
    const char vflag 
){
    static unsigned D4last=0;
    unsigned i,j, m=HK.QIDX.dim1, n=HK.QDIM;
    wbvector<TQ> q, qmin, qmax;

    wbvector<unsigned> S, I;
    wbarray<unsigned> nn;

    char fmt[8]="%3d ";

    if (!HK.isConsistent()) wblog(FL,"ERR %s",str);
    if (HK.isEmpty()) { m=n=0; }

    qmin.init(n); qmax.init(n); S.init(n);

    for (i=0; i<n; ++i) {
       HK.QIDX.getCol(i,q);
       qmin[i]=q.min(); qmax[i]=q.max();
       S[i]=unsigned(qmax[i]-qmin[i]+1);
    }

    nn.init(S); I.init(n);

    for (i=0; i<m; ++i) {
       for (j=0; j<n; ++j) I[j]=unsigned(HK.QIDX(i,j)-qmin[j]);
       nn(I)+=HK.DATA[i]->SIZE.max();
    }

#ifdef __WBDEBUG__
    i=1; 
#else
    if (vflag>=0 && (
        (E4.last()==EK) || 
        (iter<3 || iter+3>N) ||
        fabs((double(E4.len)-D4last)/MAX(unsigned(E4.len),D4last)) > 0.20
    )) i=1; else i=0;
#endif

    n=nn.sum();
    if (n!=E4.len)
         sprintf(str,"NK=%d/%ld, EK=%.2f",n, E4.len, EK);
    else sprintf(str,"NK=%d (EK=%.2f)",n, EK);

    wblog(FL,"NRG %02d: Q=[%s : %s]; %s %s",
       iter, STR(qmin), STR(qmax), str, i? "":"\r\\");
    D4last=MAX(1U,unsigned(E4.len));

    if (vflag<=10 || HK.isEmpty()) return;
    if (iter>0) checkGSDeg(E4);
    if (nn.isEmpty() || nn.SIZE.len!=2) return;

    m=nn.SIZE[0];
    n=nn.SIZE[1];

    for (i=0; i<m; ++i) { printf("\n   ");
    for (j=0; j<n; ++j) {
        if (nn(i,j)) printf(fmt,nn(i,j));
        else printf("    ");
    }}

    printf("\n\n"); fflush(0);
};

void checkGSDeg(const wbvector<double> &E4) {
   unsigned i,n;
   double dE=E4.aMin(1);

   for (n=i=0; i<E4.len; ++i) if (ABS(E4[i])<1E-10) ++n;

   if (n==1)
            wblog(FL, " *  ground state is unique (%g).", dE); else
   if (n>1) wblog(FL, " *  ground state is not unique (%d; %g).", n,dE);
   else     wblog(FL, "ERR no ground state found (%d) !?", n);
}

double getPhysE0(const wbvector<double> &DE){

   double E=0.;
   for (unsigned i=0; i<DE.len; ++i) E+=(nrgScale(i)*DE[i]);

   return E;
}

template <class TQ, class TD>
void setupPOP(
   QSpace<TQ,TD> &P,
   const QSpace<TQ,TD> &AK,
   const QSpace<TQ,TD> &UK
){
   if (P.isEmpty()) return;
   wblog(FL,"ERR");
}

