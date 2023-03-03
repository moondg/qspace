
/* CHANGE LOG ======================================================== *
 * =================================================================== */

// see also CG::getSymmetryStates.m (MatLab script)

char USAGE[]=""; 

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "setupRCStore"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
) { Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(1,-1,1); 
   QType t(FL,argin[0]);

   if (nargin==2 && Mx::IsEqual(argin[1],"info")) { 
      genRG_struct<gTQ,RTD> &B=gRS.buf[t];
      if (!B.q.isKnown()) B.checkInit(FL,t);
      argout[0]=B.toMx(0); 
      if (nargout<2) { return; } 
   }
   if (nargout) { 
   MX_CHECK_HELPER_NARGS(2,-1,0); }

   unsigned n=t.qlen(), dmax=0, npass=3;
   char sdig=1; 

   if (nargin==3 && mxIsNumeric(argin[1]) && mxIsNumeric(argin[2])) {
      qset<gTQ> q1(FL,argin[1]), q2(FL,argin[2]);

      gRS.buf[t].checkInit(FL,t);

      if (q1.len!=n || q2.len!=n) wblog(FL,
         "ERR got invalid qsets for sym=%s (len=%d,%d/%d)",
         t.toStr().data,q1.len,q2.len,n);
      if (gStore.load_RSet(0,1,t,q1)<=0) wblog(FL,
         "ERR 1st multiplet %s (%s) does not yet exist",
         t.toStr().data,q1.toStr().data);
      if (gStore.load_RSet(0,1,t,q2)<=0) wblog(FL,
         "ERR 2nd multiplet %s (%s) does not yet exist",
         t.toStr().data,q2.toStr().data);

      wblog(FL,"%N==> "
        "%s() tensor product decomposition\n==> for %s [%s] x [%s]%N",
         myname,t.toStr().data, q1.toStr().data,q2.toStr().data);
      gRS.buf[t].getTensorProdReps_gen(q1,q2,TP3_TST);

      return;
   }

   if (nargin>1) {
      OPTS opts; opts.init(argin+1, nargin-1);
      opts.getOpt("dmax",dmax);
      opts.getOpt("npass",npass);
      if (opts.getOpt("-a")) { sdig=0; }
      opts.checkAnyLeft();
   }

   if (CG_VERBOSE>0) printf("\n   setup of (defining) "
      "symmetry multiplets for %s ...\n",t.toStr().data);
   Wb::SigHandler SIG(FL);

   qset<gTQ> q0;
   gRS.getRSet(FL,t,&q0);

   printf(
   "   defining multiplet: t=[%s]\n\n", q0.toStr().data);

   unsigned n3, ip=0;
   if (int(dmax)<=0) dmax=10*t.qlen();

   for (n=0; ip<npass; ++ip) { if (CG_VERBOSE>1) {
      wblog(FL,"=== pass %d/%d (dmax=%d) %30R",ip+1,npass,dmax,"=");  }
      n3=gRS.buf[t].getTensorProdReps(dmax,sdig);
      wblog(FL,"%d multiplets generated or loaded (dmax=%d)",n3,dmax);
      SIG.check911();
   }

   printf("\n"
   "   number of CGCs: %d (%d)\n",gCS.map3[t].size(),gCS.BUF.size());
   printf(
   "   number of multiplets: %d\n",gRS.buf[t].RSet.size());

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in setupRCStore"); }
};

