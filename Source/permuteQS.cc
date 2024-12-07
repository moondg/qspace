
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
   const char *F, int L,
   const QSpace<gTQ,TD> &A, wbperm &P, char conj=0
){
   unsigned r=A.rank(F_L);
   QSpace<gTQ,TD> B; B.mt=Wb::MEX_RETURN; 

   if (P && r) { 
      if (r>P.len) { P.Extend(r); } 
      else if (r<P.len) { wblog(F,L,
         "ERR invalid permutation (%s; len=%d/%d)",STR(P),P.len,r); }
      A.permute(P,B);
   }
   else {
      if (!conj) { return A.toMx(); }
      B=A;
   }

   if (conj) { B.Conj(); B.SortDegQ(); }
   B.ctime=Wb::getTimeNow();

   return B.toMx();
};

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   wbperm P; unsigned r=-1;
   char isr=1, conj=0;

   MX_CHECK_HELPER_NARGS(2,-1,1); 
   if (nargin>3) usage(FL,"ERR invalid number of I/O arguments");

   isr=(mxIsQSpace(argin[0])>0); 
   if (!isr && mxIsQSpace(FL,argin[0],r,'c')<=0) wblog(FL,
      "ERR invalid input QSpace arg #1");

   if (mxIsChar(argin[1])) { 
      ctrIdx q(0,0,argin[1]); 
      P.init(FL,q); conj=q.conj;

      if (nargin>2) usage(FL,"ERR invalid usage "
        "(got too many input arguments for char permutation)"
      );
   }
   else {
      P.init(FL,argin[1], 1); 
      if (!P.isValidPerm()) { P.print("P");
         wblog(FL,"ERR invalid permutation (2nd arugment)");
      }
      if (nargin>2) {
         if (Mx::IsEqual(argin[2],"conj")) { conj=1; }
         else usage(FL,"ERR invalid last argument (expecting 'conj')");
      }
   }

   if (isr) {
      const QSpace<gTQ,double> A(argin[0],'r');
      argout[0]=PERMUTE_QS(FL,A,P,conj); 
   }
   else {
      const QSpace<gTQ,wbcomplex> A(argin[0],'r');
      argout[0]=PERMUTE_QS(FL,A,P,conj); 
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in permuteQS"); }
   aclu.Check();
};

