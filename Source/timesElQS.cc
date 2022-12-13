char USAGE[] =
/* ================================================================== */

"   Usage: C = timesElQS(A,B);                                    \n\
                                                                  \n\
       calculates A .* B in QSpace data.                          \n\
                                                                  \n\
   AW (C) Feb 2007                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================ */

#include "wblib.h"

#define _TQ double

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,disp,r=-1;
    mxArray *a=NULL;

    OPTS opts;

    if (nargin!=2) usage(FL,"ERR Invalid number of I/O arguments.");

    opts.init(argin+2,nargin-2);
      disp=opts.getOpt("disp");
    opts.checkAnyLeft();

    for (i=0; i<2; i++) {
       try { mxIsQSpace(FL,argin[0],r,'c'); }
       catch (...) {
          wblog(FL,"ERR arg #%d must be valid QSpace", i+1);
       }
    }

    if (mxIsQSpace(argin[0])>0 && mxIsQSpace(argin[1])>0) { 
       QSpace<double,double> A(argin[0]), B(argin[1],'r');

       if (A.QDIM==B.QDIM && A.QIDX.dim2==B.QIDX.dim2) {
          A.TimesEl(B); a=A.toMx();
       }
       else {
          if (disp) wblog(FL,"Dimension mismatch (%d,%d; %d,%d)",
          FL, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
       }
    }
    else {
       QSpace<double,wbcomplex> A(argin[0]), B(argin[1]);

       if (A.QDIM==B.QDIM && A.QIDX.dim2==B.QIDX.dim2) {
          A.TimesEl(B); a=A.toMx();
       }
       else {
          if (disp) wblog(FL,"Dimension mismatch (%d,%d; %d,%d)",
          FL, A.QDIM, B.QDIM, A.QIDX.dim2, B.QIDX.dim2);
       }
    }

    if (!a) { QSpace<int,double> Q; a=Q.toMx(); }

    argout[0]=a; 
}

