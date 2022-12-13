
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to normQS.m // Wb,Jan12,19

#define PROG_TAG "nrm"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "normQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    double nrm2; unsigned r=-1;

    MX_CHECK_HELPER_NARGS(1,-1,1); 
    if (nargin!=1) usage(FL,"ERR invalid number of I/O arguments");

    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,"valid QSpace required");

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,double> A(argin[0],'r');
       nrm2=A.norm2();
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       wbcomplex z2=A.norm2();
       if (fabs(z2.i)>1E-14) wblog(FL,
          "ERR %s() got imaginary part (%s)",FCT,z2.toStr().data);
       nrm2=z2.r;
    }

    argout[0]=numtoMx(sqrt(nrm2));

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in normQS"); }
};

