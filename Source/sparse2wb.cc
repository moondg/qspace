
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to sparse2wb.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "sparse2wb"
#endif

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(1,-1,1); 

   if (nargin!=1 || !mxIsDouble(argin[0])) wblog(FL,
      "ERR %s() invalid usage",myname);

   wbsparray<double> A(FL,argin[0]);
   argout[0]=A.toMx(); 

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in sparse2wb"); }
   aclu.Check();
};

