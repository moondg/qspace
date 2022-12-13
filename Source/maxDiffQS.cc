char USAGE[] =
/* ================================================================== */

"   Usage: Aout = maxDiffQS(A, B)                                 \n\
                                                                  \n\
       maximum difference of elements in QSpaces A and B.         \n\
       if QSpaces are incompatible, NaN is returned.              \n\
                                                                  \n\
   AW (C) Aug 2006                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,disp,r=-1;
    double dmax;

    OPTS opts;

    if (nargin<2) usage(FL,"ERR Invalid number of I/O arguments.");

   opts.init(argin+2,nargin-2);
     disp=opts.getOpt("disp");
   opts.checkAnyLeft();

    for (i=0; i<2; i++) {
       try { mxIsQSpace(FL,argin[i],r,'c'); }
       catch (...) {
          wblog(FL,"ERR invalid QSpace arg #%d", i+1);
       }
    }

    if (mxIsQSpace(argin[0])>0 && mxIsQSpace(argin[1])>0) { 
       const QSpace<double,double> A(argin[0],'r'), B(argin[1],'r');

       if (A.QDIM==B.QDIM && A.QIDX.dim2==B.QIDX.dim2) {
          dmax=A.maxDiff(B);
          if (isnan(dmax) && disp) printf("\n%s\n\n",str);
       }
       else {
          if (disp) wblog(FL,"Dimension mismatch (%d,%d; %d,%d)",
          FL, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
          dmax=NAN;
       }
    }
    else {
       const QSpace<double,wbcomplex> A(argin[0]), B(argin[1]);

       if (A.QDIM==B.QDIM && A.QIDX.dim2==B.QIDX.dim2) {
          dmax=A.maxDiff(B);
          if (isnan(dmax) && disp) printf("\n%s\n\n",str);
       }
       else {
          if (disp) wblog(FL,"Dimension mismatch (%d,%d; %d,%d)",
          FL, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
          dmax=NAN;
       }
    }

    argout[0]=numtoMx(dmax);
}

