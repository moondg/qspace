
/* CHANGE LOG ======================================================== *

   Wb,Dec18,18 : adapted to QSpace v3.

 * =================================================================== */

char USAGE[]=""; // outsourced to traceQS.m // Wb,Jan12,19

#define PROG_TAG "tr"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "traceQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TQ, class TD>
mxArray* traceQS(const QSpace<TQ,TD> &A, ctrIdx &i1, ctrIdx &i2);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned r=-1; char isr;
   ctrIdx i1,i2;

   MX_CHECK_HELPER_NARGS(1,-1,1); 
   if ((nargin!=1 && nargin!=3)) usage(FL,"ERR invalid usage");

   int q=mxIsQSpace(FL,argin[0],r,'c',-2);
   if (q<0  ) wblog(FL,"ERR %s() invalid usage (arg #1 not a QSpace",FCT);
   if (q&128) wblog(FL,"ERR %s() expecting single QSpace for arg #1",FCT);
   isr=(q&4 ? 0 : 1); 

   if (nargin>2) {
      i1.init(FL,argin[1]);
      i2.init(FL,argin[2]);
   }

   if (i1.len!=i2.len) wblog(FL,"ERR %s() "
      "invalid index sets i1=%s <> i2=%s",PROG,STR(i1),STR(i2));

   if (isr) { 
      const QSpace<gTQ,double> A(argin[0],'r');
      argout[0]=traceQS(A,i1,i2); 
   }
   else {
      const QSpace<gTQ,wbcomplex> A(argin[0]);
      argout[0]=traceQS(A,i1,i2); 
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in traceQS"); }
   aclu.Check();
};

template<class TQ, class TD>
mxArray* traceQS(const QSpace<TQ,TD> &A, ctrIdx &i1, ctrIdx &i2) {

   if (A.isEmpty()) {
      QSpace<TQ,TD> A2;
      return A2.initScalar(0).toMx();
   }

   unsigned r=A.rank();

   if (!i1.len && !i2.len) {
      if (r%2) wblog(FL,"ERR %s() invalid usage (got odd rank r=%d)",FCT,r);
      return numtoMx(A.trace());
   }
   else {
      QSpace<TQ,TD> B; B.mt=Wb::MEX_RETURN; 
      A.trace(FL,i1,i2,B);

      if (B.isScalar())
           { return B.DATA[0]->toMx(); }
      else { return B.toMx(); }
   }
};

