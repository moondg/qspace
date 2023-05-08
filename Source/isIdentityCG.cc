
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to isIdentityCG.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "isIdentityCG"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned i,r=-1; double eps=1E-12;

    MX_CHECK_HELPER_NARGS(1,-1,-1); 

    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
    "first argument requires valid QSpace");

    if (nargin>1)
    if (mxGetNumber(argin[1], eps)) wbdie(FL,str);

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,double> A(argin[0],'r');
       i=A.hasIdentityCGS(eps) ? 1 : 0;
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       i=A.hasIdentityCGS(eps) ? 1 : 0;
    }

    argout[0]=numtoMx(i); 

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in isIdentityCG"); }
};

