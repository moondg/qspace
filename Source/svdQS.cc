
/* CHANGE LOG ======================================================== *

   added "Nkmin" (as follow up to discussion with Wei Li) // Wb,Sep16,17

 * =================================================================== */

char USAGE[]=""; // outsourced to svdQS.m // Wb,Jan12,19

#define PROG_TAG "svd"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "svdQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TD>
void MPS_ORTHO_1(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpace<gTQ,TD> &PSI
);

template<class T1, class T2, class T3>
void MPS_ORTHO_2(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpace<gTQ,T1> &A1,
   const QSpace<gTQ,T2> &A2,
   QSpace<gTQ,T3> &PSI
);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned got12, r=-1;

   MX_CHECK_HELPER_NARGS(2,-1,4); 

   mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
      "invalid input QSpace arg #1");

   r=-1; 
   got12 = (mxIsQSpace(0,0,argin[1],r,'c')>0);

   if (got12) wblog(FL,"WRN %s() using deprecated usage #2",myname);
   if (got12 && nargin<4) wblog(FL,
      "ERR too few input args (%d/4)",nargin);

   if (got12) {
      const unsigned
      isra=(mxIsQSpace(argin[0])>0), 
      isrb=(mxIsQSpace(argin[1])>0); 

      if (isra && isrb) {
         const QSpace<gTQ,double> A1(argin[0],'r');
         const QSpace<gTQ,double> A2(argin[1],'r');
         QSpace<gTQ,double> PSI;
         MPS_ORTHO_2(nargout,argout,nargin,argin,A1,A2,PSI);
      }
      else {
         const QSpace<gTQ,wbcomplex> A1(argin[0]);
         const QSpace<gTQ,wbcomplex> A2(argin[1]);
         QSpace<gTQ,wbcomplex> PSI;
         MPS_ORTHO_2(nargout,argout,nargin,argin,A1,A2,PSI);
      }
   }
   else {
      const unsigned
      isra=(mxIsQSpace(argin[0])>0); 

      if (isra) {
         QSpace<gTQ,double> PSI(argin[0],'r');
         MPS_ORTHO_1(nargout,argout,nargin,argin,PSI);
      }
      else {
         QSpace<gTQ,wbcomplex> PSI(argin[0],'r');
         MPS_ORTHO_1(nargout,argout,nargin,argin,PSI);
      }
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in svdQS"); }
   aclu.Check();
};

template<class T1, class T2, class T3>
void MPS_ORTHO_2(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   const QSpace<gTQ,T1> &A1,
   const QSpace<gTQ,T2> &A2,
   QSpace<gTQ,T3> &PSI
){
   unsigned i,ic1,ic2,r1,r2, K, vflag=0, Nkmin=0; int Nkeep=-1;
   double stol=WB_STOL_SVD;
   mxArray *Sout;
   itag_ tx;

   QSpace<gTQ,T3> UQ,SQ,VQ;
   OPTS opts;

   if (A1.isEmpty()) wblog(FL,"ERR arg #1 (A1) is empty");
   if (A2.isEmpty()) wblog(FL,"ERR arg #2 (A2) is empty");

   r1=A1.rank(); r2=A2.rank();

   for (i=2; i<4; ++i) if (!Mx::IsNumber(0,0,argin[i]))
      wblog(FL,"ERR input arg #%d must be a valid index.", i+1);

   if (mxGetNumber(argin[2], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[3], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
      --ic1; --ic2; 
   }
   else wblog(FL,
      "ERR arguments #3,#4 - index out of bounds (%d,%d)",ic1,ic2);

   opts.init(argin+4,nargin-4);

      opts.getOpt("Nkmin",Nkmin); 
      opts.getOpt("Nkeep",Nkeep); if (Nkeep>0) { stol=0; }
      opts.getOpt("stol", stol ); 
      opts.getOpt("itag", tx   ); 
      vflag=opts.getOpt("-v");

   opts.checkAnyLeft();

   twoSiteInit(A1,A2,PSI, K, ic2, ic1, 0,-1);

   Sout=getSVD(
      PSI, UQ, SQ, VQ, K,
      Nkmin,Nkeep,stol,tx,vflag ? "dS":"s");

   twoSiteFinal(A2, A1, VQ, UQ, ic2, ic1, 0,-1);

   argout[0]=UQ.toMx(); if (nargout>1) {
   argout[1]=SQ.toMx(); if (nargout>2) {
   argout[2]=VQ.toMx(); if (nargout>3) { argout[3]=Sout; Sout=0; }}}

   if (Sout) mxDestroyArray(Sout);
};

template<class TD>
void MPS_ORTHO_1(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   QSpace<gTQ,TD> &PSI
){
   unsigned i,r,K, vflag, Nkmin=0; int Nkeep=-1;
   double stol=WB_STOL_SVD;
   itag_ tx;

   wbindex I; wbperm P,iP,pa,pb;
   wbvector<unsigned> D;

   QSpace<gTQ,TD> UQ,SQ,VQ;
   const mxArray *a; mxArray *b, *Sout;
   OPTS opts;

   if (!PSI) wblog(FL,"ERR got empty PSI");
   r=PSI.rank(FL);

   a=argin[1]; {
      unsigned m=mxGetM(a), n=mxGetN(a);

      if (!Mx::IsDblMat(0,0,a) || (m>1 && n>1)) wblog(FL,
         "ERR (%d,%d)\n%s",m,n,str);

      I.init(FL,argin[1],"invalid index (arg #2)",'T'); 
      if (!I.len || I.len>=r) wblog(FL,
         "ERR invalid index set (arg #2) (len=%d/%d)",I.len,r);
      for (i=0; i<I.len; ++i) { if (I[i]<1 || I[i]>r) wblog(FL,
         "ERR index out of bounds (%d/%d)",I[i],r); }
      I-=1; 
   }

   opts.init(argin+2,nargin-2);

      opts.getOpt("Nkmin",Nkmin); 
      opts.getOpt("Nkeep",Nkeep); if (Nkeep>0) { stol=0; }
      opts.getOpt("stol", stol ); 
      opts.getOpt("itag", tx   ); 
      vflag=opts.getOpt("-v");

      opts.getOpt("permA",b); if (b) pa.init(FL,b,1);
      opts.getOpt("permB",b); if (b) pb.init(FL,b,1);

   opts.checkAnyLeft();

   if (pa.len>r || pb.len>r || (pa.len && pb.len && pa.len+pb.len!=r+2))
      wblog(FL,"ERR invalid perms [%s, %s; %d+2]",pa.len,pb.len);

   if (I.extend2Perm(r,P,'E')) 
      wblog(FL,"ERR invalid index set%N%N%s%N",str);
   K=r-I.len; iP.init(P,'i');

   PSI.Permute(P,0,'c');

   Sout=getSVD(
      PSI, UQ, SQ, VQ, K,  
      Nkmin,Nkeep,stol,tx,vflag ? "dS":"s");

   P.initLastTo(0,r-K+1); {
      if (pb.len) { P.Permute(pb); }
      P.save2(pb);
   }

   if (K+1==r) {
      if (pa.len) iP.Permute(pa);
      iP.save2(pa);
   }
   else if (K==1) {
      if (pb.len) iP.Permute(pb);
      iP.save2(pb);
   }

   if (pa) { UQ.Permute(pa); }
   if (pb) { VQ.Permute(pb); }

   argout[0]=UQ.toMx(); if (nargout>1) {
   argout[1]=SQ.toMx(); if (nargout>2) {
   argout[2]=VQ.toMx(); if (nargout>3) { argout[3]=Sout;  Sout=0; }}}

   if (Sout) mxDestroyArray(Sout);
};

