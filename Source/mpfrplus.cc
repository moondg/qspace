
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // Wb,Oct10,25

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpfrplus"
#endif

   #define LD_CLEBSCH_QS

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin,  const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   char dflag=0;
   double a,b;

   MX_CHECK_HELPER_NARGS(4,5,1); 

   if (!Mx::IsNumber(0,0,argin[0]) ||
       !Mx::IsNumber(0,0,argin[2])) { usage(FL,"invalid usgae"); }

   if (nargin==5) { char s[8]; s[0]=0;
      if (!mxIsChar(argin[4]) || mxGetString(argin[4],s,8) ||
           strcmp(s,"-d")) { usage(FL,"invalid usgae (s=%s)",s); }
      dflag=1;
   }

   wbsparray<Wb::quad> A(argin[1]), B(argin[3]); 
   mxGetNumber(argin[0],a);
   mxGetNumber(argin[2],b);

   A.Plus(FL,B,b,a);

   if (dflag) {
      if (A.SIZE.len<=2)
           { wbsparray<double> X; X.init(A); argout[0]=X.toMx(); }
      else { wbarray<double> X; A.toFull(X); argout[0]=X.toMx(); }
   }
   else { argout[0]=A.toMx(); }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in mpfrplus"); }
   aclu.Check();
};

