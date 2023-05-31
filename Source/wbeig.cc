
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to wbeig.m // Wb,Feb14,19

#include "wblib.h"

template<class T>
void WB_EIG(wbarray<T> &H, mxArray *argout[], int nargout);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    MX_CHECK_HELPER_NARGS(1,-1,2); 
    if (nargin!=1) usage(FL,"ERR invalid number of I/O arguments");

    if (!mxIsComplex(argin[0])) {
       wbarray<double> H(argin[0]);
       WB_EIG(H,argout,nargout);
    }
    else {
       wbarray<wbcomplex> H(argin[0]);
       WB_EIG(H,argout,nargout);
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbeig"); }
   aclu.Check();
};

template<class T>
void WB_EIG(wbarray<T> &H, mxArray *argout[], int nargout) {

    wbvector<double> E;
    wbarray<T> U;

    wbEigenS(H,U,E);

    argout[0]=E.toMx(); if (nargout>1) {
    argout[1]=U.toMx(); }
};

