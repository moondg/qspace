
/* CHANGE LOG ======================================================== *

// Wb,Dec09,14 :: added optional conj as trailing arument or string.

   NB! about x2 faster than MatLab based permuteQ() for PSI =
      Q: {[197x2 double]  [197x2 double]  [197x2 double]  [197x2 double]}
   data: {197x1 cell}
   Wb,Aug07,06

 * =================================================================== */

char USAGE[]=""; // outsourced to permuteQS.m // Wb,Jan12,19

#define PROG_TAG "pmt" 

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "permuteQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TD>
mxArray* PERMUTE_QS(
   const char *F, int L, const QSpace<gTQ,TD> &A, wbperm &P){

   unsigned r=A.rank(F_L);

   if (!r || !P.relevant()) { return A.toMx(); }
   else {
      QSpace<gTQ,TD> B; B.mt=Wb::MEX_RETURN; 

      if (r>P.len) { P.Extend(r); } 
      else if (r<P.len) wblog(F,L,
         "ERR invalid permutation (%s; len=%d/%d)",STR(P),P.len,r);

      A.permute(B,P);
      B.NormCGW();  

      B.ctime=Wb::getTimeNow();

      return B.toMx();
   }
};

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   wbperm P; unsigned r=-1;
   char isr=1;

   MX_CHECK_HELPER_NARGS(2,3,1); 
   if (nargin>3) usage(FL,"ERR invalid number of I/O arguments");

   isr=(mxIsQSpace(argin[0])>0); 
   if (!isr && mxIsQSpace(FL,argin[0],r,'c')<=0) wblog(FL,
      "ERR invalid input QSpace arg #1");

   P.init(FL,argin[1]); 

   if (nargin>2) { 
      if (nargin==3 && Mx::IsEqual(argin[2],"conj")) {
         if (++P.conj>1) { wblog(FL,
         "ERR %s() invalid usage (`conj' specified twice)",FCT); }
      }
      else usage(FL,"invalid usage (got %d input arguments)",nargin);
   }

   if (isr) {
      const QSpace<gTQ,double> A(argin[0],'r');
      argout[0]=PERMUTE_QS(FL,A,P); 
   }
   else {
      const QSpace<gTQ,wbcomplex> A(argin[0],'r');
      argout[0]=PERMUTE_QS(FL,A,P); 
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in permuteQS"); }
   aclu.Check();
};

