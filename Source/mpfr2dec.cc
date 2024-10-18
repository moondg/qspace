
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to mpfr2dec.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpfr2dec"
#endif

   #define LOAD_CGC_QSPACE

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin,  const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned base=0; char tflag=0;
   wbvector<Wb::quad> q;

   MX_CHECK_HELPER_NARGS(1,-1,1); 
   if (nargin>3) usage(FL,"invalid number of I/O arguments");

   for (unsigned i=1; int(i)<nargin; ++i) {
      if (mxIsChar(argin[i]) && !tflag) {
         mxGetString(argin[i],str,128);
         if (!strcmp(str,"-t")) tflag=1;
         else { 
            snprintf(str+128,128,"invalid option '%.64s'",str);
            usage(FL,str);
         }
      }
      else if (!base) {
        if (mxGetNumber(argin[i],base)) wblog(FL,"ERR %s()",FCT);
        if (base<2 || base>62) wblog(FL,
           "ERR %s() invalid base=%d",FCT,base);
      }
      else usage(FL,"invalid usage");
   }

   if (mxIsChar(argin[0]) || mxIsInt8(argin[0])) {
      q.init_mpfr(FL,argin[0],base);
      argout[0]=wbvector<double>(q).toMx(tflag); 
   }
   else if (mxIsNumeric(argin[0])) {
      wbarray<double> A(FL,argin[0]); if (tflag) { A.Transpose(); }
      argout[0]=A.toMx();
   }
   else if (mxIsCell(argin[0])) {
      const mxArray *a=argin[0];
      const size_t *sp=mxGetDimensions(a);
      unsigned i=0, r=mxGetNumberOfDimensions(a), n=mxGetNumberOfElements(a);

      mxArray *C=mxCreateCellArray(r,sp);

      for (; i<n; ++i) {
          try {
             q.init_mpfr(FL,mxGetCell(a,i),base);
          }
          catch (...) {
             wblog(FL,"ERR %s() invalid cell %d/%d",FCT,i+1,n);
          }
          mxSetCell(C,i, wbvector<double>(q).toMx(tflag));
      }
      argout[0]=C;
   }
   else {
      sprintf_str("expecting (cell of) "
         "strings as input #1 (%s)", mxGetClassName(argin[0]));
      usage(FL,str);
   }

   return;

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in mpfr2dec"); }
   aclu.Check();
};

