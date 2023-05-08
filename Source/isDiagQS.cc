char USAGE[] =
/* ================================================================== */

"   Usage: isdiag = isDiagQS(A)                                   \n\
                                                                  \n\
       check whether QSpace specified by A represents a           \n\
       block-diagonal matrix.                                     \n\
                                                                  \n\
   AW (C) May 2006                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    double eps=1E-14;
    unsigned i, r=-1;

    if (nargin<1) usage(FL,"ERR Invalid number of I/O arguments.");

    mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
    "first argument requires valid QSpace");

    if (nargin>1)
    if (mxGetNumber(argin[1], eps)) wbdie(FL,str);

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<double,double> A(argin[0],'r');
       i=A.isDiagMatrix(eps) ? 1 : 0;
    }
    else {
       const QSpace<double,wbcomplex> A(argin[0]);
       i=A.isDiagMatrix(eps) ? 1 : 0;
    }

    argout[0]=numtoMx(i); 
}

