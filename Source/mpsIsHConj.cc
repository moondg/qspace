
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to mpsIsHConj.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsIsHConj"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    char vflag=0; int i; double eps=1E-12;

    MX_CHECK_HELPER_NARGS(1,-1,2); 

    if (mxIsQSpace(FL,argin[0],'c')<=0) wblog(FL,
       "ERR invalid QSpace arg #1");

    if (nargin>1) {
       wbstring mark(nargin);
       for (i=1; i<nargin; ++i) {
          if (mxIsChar(argin[i])) {
             wbstring s(argin[i]);
             if (s=="-v") vflag='v';
             else wbdie(FL,str);
          }
          else if (Mx::IsNumber(0,0,argin[i])) {
             if ((++mark[i])>1 || mxGetNumber(argin[1], eps))
             wbdie(FL,str);
          }
       }
    }

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,double> A(argin[0],'r');
       i=A.isHConj(0,0,eps,vflag);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       i=A.isHConj(0,0,eps,vflag);
    }

    argout[0]=numtoMx(i); 
    if (nargout>1) argout[1]=mxCreateString(i ? "" : str);

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in mpsIsHConj"); }
   aclu.Check();
};

