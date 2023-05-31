
/* CHANGE LOG ======================================================== *

   added "Nkmin" (as follow up to discussion with Wei Li) // Wb,Sep16,17

 * =================================================================== */

char USAGE[]=""; // outsourced to orthoQS.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "orthoQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

template<class TQ, class TD>
void MPS_ORTHO_1(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI, const ctrIdx &I    
);

template<class TQ, class TD>
void MPS_ORTHO_2(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI, ctrIdx &I, char ldir 
);

template<class TQ, class TD>
void MPS_ORTHO_3(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   const QSpace<TQ,TD> &Psi1,
   const QSpace<TQ,TD> &Psi2, QSpace<TQ,TD> &PSI 
);

template<class TQ, class TD>
int get_Nkeep_estimate(QSpace<TQ,TD> &PSI, unsigned K) {

   int Nkeep=-1, d1=1, d2=1;
   wbvector<widx_t> dd;

   PSI.getDim(dd); {
      if (K>dd.len) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",myname,K,dd.len);
      d1=Wb::prodRange(dd.data,K);
      d2=Wb::prodRange(dd.data+K,dd.len-K);
   }

   if (PSI.allAbelian())
        Nkeep=(d1<d2 ? d1 : d2); 
   else Nkeep=(d1<d2 ? d2 : d1); 

   return Nkeep;
};

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   char got2=0, ldir=0;
   unsigned r=-1, isra;
   ctrIdx I;  

   MX_CHECK_HELPER_NARGS(2,-1,3); 

   isra=(mxIsQSpace(argin[0])>0); if (!isra) 
   mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
      "argument #1 must be a valid QSpace");

   r=-1; 
   if (mxIsQSpace(0,0,argin[1],r,'c')>0) { got2=1; }
   else {
      I.init(FL,argin[1]);
      if (nargin>2 && mxIsChar(argin[2])) {
         ldir=mxGetLRDir(0,0,argin[2],3); 
      }
   }

   if (!got2 && !ldir) { 
      if (I.len!=1) wblog(FL,
         "ERR %s() invalid usage #1 (idx: '%s')",myname,STR(I));
      if (isra) {
         QSpace<gTQ,double> PSI(argin[0],'r');
         MPS_ORTHO_1(nargout,argout,nargin-2,argin+2,PSI,I);
      }
      else {
         QSpace<gTQ,wbcomplex> PSI(argin[0]);
         MPS_ORTHO_1(nargout,argout,nargin-2,argin+2,PSI,I);
      }
   }
   else if (!got2) { 
      if (isra) { 
         QSpace<gTQ,double> PSI(argin[0],'r');
         MPS_ORTHO_2(nargout,argout,nargin-3,argin+3,PSI,I,ldir);
      }
      else {
         QSpace<gTQ,wbcomplex> PSI(argin[0]);
         MPS_ORTHO_2(nargout,argout,nargin-3,argin+3,PSI,I,ldir);
      }
   }
   else { 
      unsigned isrb=(mxIsQSpace(argin[1])>0); 

      if (nargin<5) wblog(FL,
         "ERR usage #3: invalid number of input args (%d)",nargin);
      if (qsGotCGS(FL,argin[0])) wblog(FL, 
         "WRN %s() usage #3 for abelian symmetries only",myname);

      if (isra && isrb) {
         const QSpace<gTQ,double> Psi1(argin[0],'r');
         const QSpace<gTQ,double> Psi2(argin[1],'r');
         QSpace<gTQ,double> PSI;
         MPS_ORTHO_3(nargout,argout,nargin-2,argin+2,Psi1,Psi2,PSI);
      }
      else {
         const QSpace<gTQ,wbcomplex> Psi1(argin[0]);
         const QSpace<gTQ,wbcomplex> Psi2(argin[1]);
         QSpace<gTQ,wbcomplex> PSI;
         MPS_ORTHO_3(nargout,argout,nargin-2,argin+2,Psi1,Psi2,PSI);
      }
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in orthoQS"); }
   aclu.Check();
};

template<class TQ, class TD>
void MPS_ORTHO_1(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI,
   const ctrIdx &I
){
   unsigned K, vflag=1, r=PSI.rank(FL), Nkmin=0; int Nkeep=-1;
   double stol=WB_STOL_SVD;
   mxArray *S;
   itag_ tx;

   QSpace<TQ,TD> A1, A2;
   wbvector<widx_t> D;
   wbperm P,iP,pb;
   OPTS opts;

   if (PSI.isEmpty()) wblog(FL,"ERR input arg #1 (PSI) is empty");
   if (I.len!=1 || I.data[0]>r) wblog(FL,
      "ERR invalid index set arg #2 (%s/%d)",STR(I),r);

   opts.init(argin,nargin);

      opts.getOpt("Nkmin",Nkmin); 
      opts.getOpt("Nkeep",Nkeep); if (Nkeep>0) { stol=0; }
      opts.getOpt("stol", stol );
      opts.getOpt("itag", tx   ); 

      if (opts.getOpt("-v")) { vflag=2; } else
      if (opts.getOpt("-q")) { vflag=0; } 

   opts.checkAnyLeft();

   if (I.extend2Perm(r,P)) wblog(FL,
      "ERR invalid index set%N%N%s%N",str);
   iP.init(P,'i');

   PSI.Permute(P,0,'i'); K=I.len;

   if (Nkeep<0) {
   Nkeep=get_Nkeep_estimate(PSI,K); }

   S=orthoQS(PSI,A1,A2, K, Nkmin,Nkeep,stol,tx,vflag>1 ? "dS":"s");

   pb.initLastTo(0,r); pb.Permute(iP);
   A2.Permute(pb);

   argout[0]=A1.toMx();
   argout[1]=A2.toMx(); if (nargout>2) { argout[2]=S; S=0; }

   if (S) mxDestroyArray(S);
};

template<class TQ, class TD>
void MPS_ORTHO_2(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   QSpace<TQ,TD> &PSI, ctrIdx &I, char ldir
){
   unsigned r,K, vflag=1, Nkmin=0; int Nkeep=-1;
   wbperm P,iP,pa,pb;
   wbvector<widx_t> D;
   double stol=WB_STOL_SVD;
   itag_ tx;

   QSpace<TQ,TD> A1, A2;
   mxArray *S;
   OPTS opts;

   if (PSI.isEmpty()) wblog(FL,"ERR input arg #1 (PSI) is empty");
   r=PSI.rank(FL);

   opts.init(argin,nargin);

      opts.getOpt("Nkmin",Nkmin); 
      opts.getOpt("Nkeep",Nkeep); if (Nkeep>0) { stol=0; }
      opts.getOpt("stol", stol ); 
      opts.getOpt("itag", tx   ); 

      if (opts.getOpt("-v")) { vflag=2; } else
      if (opts.getOpt("-q")) { vflag=0; } 

      opts.getOpt("permA",S); if (S) pa.init(FL,S,1);
      opts.getOpt("permB",S); if (S) pb.init(FL,S,1);

      if (pa.len>r || pb.len>r || (pa.len && pb.len && pa.len+pb.len!=r+2))
      wblog(FL,"ERR invalid permutations perm[AB]");

   opts.checkAnyLeft();

   if (!I.len || I.len>=r) wblog(FL,
      "ERR invalid index set arg #2 (%d/%d)",I.len,r);

   if (abs(ldir)==3) { 
   if (I.len==1 || I.len+1==r) {
      I.Invert(r); ldir=-ldir;
   }}

   if (ldir>=0) { 
      if (I.extend2Perm(r,P)) wblog(FL,
         "ERR invalid index set%N%N%s%N",str);
      iP.init(P,'i');

      PSI.Permute(P,0,'i'); K=I.len; 

      if (Nkeep<0) {
      Nkeep=get_Nkeep_estimate(PSI,K); }

      S=orthoQS(PSI,A1,A2, K, Nkmin,Nkeep,stol,tx,vflag>1 ? "dS":"s");

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

      if (!pa.isIdentityPerm()) A1.Permute(pa); 
      if (!pb.isIdentityPerm()) A2.Permute(pb);

      argout[0]=A1.toMx();
      argout[1]=A2.toMx();
   }
   else {
      char ipb=0; 

      if (I.extend2Perm(r,P,'E')) wblog(FL, 
         "ERR invalid index set%N%N%s%N",str);
      iP.init(P,'i');

      PSI.Permute(P,0,'i'); K=r-I.len; 

      if (Nkeep<0) {
      Nkeep=get_Nkeep_estimate(PSI,K); }

      S=orthoQS(PSI,A2,A1, K, Nkmin,Nkeep,stol,tx,vflag>1 ? "dS":"s");

      if (I.len+1==r) {
         P.initLastTo(0,I.len+1);
         P.Permute(iP); if (pa.len) { P.Permute(pa); }
         P.save2(pa); 
      }
      else if (I.len==1) {
         if (pb.len) iP.Permute(pb);
         iP.save2(pb); ipb=1;
      }

      if (!ipb) {
         P.initLastTo(0,K+1); if (pb.len) { P.Permute(pb); }
         P.save2(pb);
      }

      if (!pa.isIdentityPerm()) A1.Permute(pa); 
      if (!pb.isIdentityPerm()) A2.Permute(pb);

      argout[0]=A1.toMx();
      argout[1]=A2.toMx();
   }

   if (nargout>2) { argout[2]=S; S=0; }
   if (S) mxDestroyArray(S);
};

template<class TQ, class TD>
void MPS_ORTHO_3(
   int nargout, mxArray *argout[], int nargin, const mxArray *argin[],
   const QSpace<TQ,TD> &Psi1,
   const QSpace<TQ,TD> &Psi2,
   QSpace<TQ,TD> &PSI
){
   unsigned i,ic1,ic2,r1,r2, K, vflag=1, Nkmin=0; int Nkeep=-1;
   double stol=WB_STOL_SVD;
   char ldir=0;
   mxArray *S;
   itag_ tx;

   QSpace<TQ,TD> A1, A2;
   OPTS opts;

   if (Psi1.isEmpty()) wblog(FL,"ERR input arg #1 is empty (Psi1)");
   if (Psi2.isEmpty()) wblog(FL,"ERR input arg #2 is empty (Psi2)");

   r1=Psi1.rank(); r2=Psi2.rank();

   for (i=0; i<2; ++i) {
      if (!Mx::IsNumber(0,0,argin[i]))
      wblog(FL,"ERR input arg #%d must be a valid index", i+1);
   }

   if (mxGetNumber(argin[0], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[1], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
       ic1--; ic2--; 
   }
   else wblog(FL,"ERR args #3 and #4\n"
      "index out of bounds (%d/%d, %d/%d)",ic1,r1,ic2,r2);

   ldir=mxGetLRDir(FL,argin[2],4); 

   opts.init(argin+3,nargin-3); {
      opts.getOpt("Nkmin",Nkmin); 
      opts.getOpt("Nkeep",Nkeep); if (Nkeep>0) { stol=0; }
      opts.getOpt("stol", stol ); 
      opts.getOpt("itag", tx   ); 

      if (opts.getOpt("-v")) { vflag=2; } else
      if (opts.getOpt("-q")) { vflag=0; } 
   }
   opts.checkAnyLeft();

   twoSiteInit(Psi1,Psi2,PSI, K,ic1,ic2, 0,ldir); 

   if (Nkeep<0) {
   Nkeep=get_Nkeep_estimate(PSI,K); }

   S=orthoQS(PSI,A1,A2,K, Nkmin,Nkeep,stol,tx,vflag>1 ? "dS":"s");

   twoSiteFinal(Psi1,Psi2, A1,A2, ic1,ic2, 0,ldir); 

   argout[0]=A1.toMx(); if (nargout>1) {
   argout[1]=A2.toMx(); if (nargout>2) { argout[2]=S; S=0; }}

   if (S) mxDestroyArray(S);
};

