char USAGE[] =
/* ==================================================================
 * getDRangeQS.cc */

"   Usage: rr = getDRangeQS(A)                                    \n\
                                                                  \n\
       get data range for given QSpace object                     \n\
       where rr=[dmin,dmax]                                       \n\
                                                                  \n\
   AW (C) May 2006                                                \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    if (nargin<1 || nargout>1)
    usage(FL,"ERR Invalid number of I/O arguments.");

    try { mxIsQSpace(FL,argin[0]); }
    catch (...) { wblog(FL,"ERR Input not valid QSpace object."); }

    const QSpace<double,double> A(argin[0],'r');
    wbvector<double> dd(2);

    A.getDRange(dd[0],dd[1]);
    argout[0]=dd.toMx(); 
}

