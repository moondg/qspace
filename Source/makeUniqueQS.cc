
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to makeUniqueQS.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mpsMakeUnique"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned r=-1; 

    wbMatrix<gTQ> Q;
    OPTS opts;

    str[0]=0;

    MX_CHECK_HELPER_NARGS(1,-1,1); 

   opts.init(argin+1,nargin-1);

   opts.checkAnyLeft();

    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
    "input not valid QSpace object");

    if (mxIsQSpace(argin[0])>0) { 
       QSpace<gTQ,double> A(argin[0]);
       A.MakeUnique();
       argout[0]=A.toMx(); 
    }
    else {
       QSpace<gTQ,wbcomplex> A(argin[0]);
       A.MakeUnique();
       argout[0]=A.toMx(); 
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in makeUniqueQS"); }
   aclu.Check();
};

