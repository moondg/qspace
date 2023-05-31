
/* CHANGE LOG ======================================================== *

 * =================================================================== */

char USAGE[]=""; // outsourced to getDimQS.m // Wb,Jan12,19

#define PROG_TAG "dim"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getDimQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    wbvector<widx_t> D,D2;
    char cgflag=0; 

    MX_CHECK_HELPER_NARGS(1,1,2); 

    str[0]=0;

    try { mxIsQSpace(FL,argin[0],'c'); }
    catch (...) {
       if (str[0]) printf("\n%s\n\n",str);
       wblog(FL,"ERR input not a valid QSpace object.");
    }

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,gTD> A(argin[0],'r',0); 
       A.getDim(D,&D2); cgflag=A.gotCGS(FL);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0],'r',0); 
       A.getDim(D,&D2); cgflag=A.gotCGS(FL);
    }

    if (nargout==2) {
        argout[0]=D .toMx_base_d();
        argout[1]=D2.toMx_base_d();
    }
    else if (cgflag || D!=D2) {
       wbvector< wbvector<widx_t> const* > dd(2); dd[0]=&D; dd[1]=&D2;
       argout[0]=wbMatrix<widx_t>().CAT(1,dd).toMx_base_d();
    }
    else argout[0]=D.toMx_base_d(); 

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getDimQS"); }
   aclu.Check();
};

