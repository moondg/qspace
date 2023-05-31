
/* CHANGE LOG ======================================================== *

 * =================================================================== */

char USAGE[]=""; // outsourced to compactQS.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "compactQS"
#endif

#define LOAD_CGC_QSPACE

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(5,-1,2); 

   unsigned i,k=0;

   if (Mx::IsEqual(argin[k],"--stat")) { 
      argout[0]=gCS.toMx(); 
      return;
   }

   if (nargin>6 || !mxIsChar(argin[k]) || 
      !Mx::IsDblMat(0,0,argin[k+1]) || !Mx::IsDblMat(0,0,argin[k+2]) ||
      !Mx::IsDblMat(0,0,argin[k+3]) || !Mx::IsDblArr(0,0,argin[k+4])
   ){
       if (nargin || k || nargout) wblog(FL,"ERR invalid usage");
       else { usage(); return; }
   }

   QVec qvec(FL,argin[k++]);
   wbMatrix<gTQ> Q1(FL,argin[k++]), Q2(FL,argin[k++]), Q(FL,argin[k++]);
   wbarray<double> D3(FL,argin[k++]); D3.appendSingletons(3);

   QSpace<gTQ,double> A;
   wbperm P; if (nargin>5) P.init(FL,argin[k++],1); 
   wbvector< wbMatrix<gTQ>* > qq(3);
   qset<gTQ> qz;

   char isa;

   qq[0]=&Q1; qq[1]=&Q2; qq[2]=&Q;
   if (!P.isEmpty()) {
      if (P.len!=3) wblog(FL,"ERR invalid permutation length (%d/3)",P.len);
      qq.Permute(P); D3.Permute(P);
   }

   for (i=0; i<3; ++i)
   Wb::FixRational(FL,qq[i]->data,qq[i]->numel(),'r',4,1024);

   A.initOpZ_WET(FL,qvec,*(qq.data[0]),*(qq.data[1]),*(qq.data[2]),D3);

   isa=A.qtype.allAbelian();

   (A.qtype.len ? A.qtype : qvec).getScalar(qz,'z'); 

   if (qz.len!=Q.dim2) wblog(FL,
      "ERR %s() q-label length mismatch (%d/%d)",myname,Q.dim2,qz.len);

   if (Q.recAllEqual(qz.data)) {

      wbMatrix<gTQ> QQ(A.QIDX); QQ.getBlocks(0,1,A.QDIM,A.QIDX);
      wbarray<double> **dd=A.DATA.data;
      CRef<gTQ> *cg=A.CGR.data;
      unsigned n=A.CGR.numel(), ok=1;

      for (i=0; i<A.DATA.len; ++i) {
         if (dd[i]->SIZE.len!=3 || dd[i]->SIZE.data[2]!=1) {
            wblog(FL,
              "WRN failed to reduce to scalar op (%d)", dd[i]->SIZE.len,
               dd[i]->SIZE.len==3 ? dd[i]->SIZE.data[2] : -1);
            ok=0; break;
         }
      }

      if (ok) {
      for (i=0; i<n; ++i) { const CRef<gTQ> &ci=cg[i];
         if (!ci.cgb) { wblog(FL,
            "WRN %s() got cgb=null (%s)",FCT,SSTR(ci.cgw));
            continue;
         }

         if (ci.wnumel()!=1 || ci.cgb->qdir.len!=3) wblog(FL,
            "ERR %s() invalid CRef data\n%s",FCT,STR(ci));
         if (ci.Size(2)!=1) { wblog(FL,
            "ERR failed to reduce to scalar op\n%s having P=[%s]",
            STR(ci), STR(ci.cgp));
         }
      }

      for (i=0; i<A.DATA.len; ++i) { dd[i]->SIZE.Shorten2(2); }
      for (i=0; i<n; ++i) {
         cg[i].Reduce2Identity(isa? 0:'!');
      }

      A.otype=QS_NONE;
      A.itags.len=2;
   }}

   A.NormCGW();

   argout[0]=A.toMx(); 

#ifdef __WBDEBUG__
   Wb::MemCheck(FL);

   qvec.init(); P.init(); qq.init();
   Q1.init(); Q2.init(); Q.init(); D3.init();
   gCS.init(); A.init();

   Wb::MemCheck(FL);
#endif

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in compactQS"); }
   aclu.Check();
};

