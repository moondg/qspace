char USAGE[] =
/* ================================================================== */
"   Usage: [A,B] = mpsTransferSpace(A,B,dim[,W,eps])              \n\
                                                                  \n\
       shift space between two A-tensors such as AK and AD        \n\
       with respect to dimenion dim.                              \n\
                                                                  \n\
       Without [W,eps], all space of A is transferred to B        \n\
       returning the result as single output argument.            \n\
                                                                  \n\
       With [W,eps] specicified, W must be diagonal QSpace        \n\
       describing weights, where only w<eps will be transferred   \n\
       from QSpace A to QSpace B. If B is empty upon input,       \n\
       then A is split into kept (A) / discarded (B).             \n\
       Both, the updated A and B are returned.                    \n\
                                                                  \n\
   (C) Wb,Feb24,11                                                \n";

/* ================================================================== */

// This is a MEX-file for MATLAB
// tags: shiftspace movespace mergespace combinespace

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    unsigned i,dim=0,r=-1;
    double eps;
    mxArray *a=NULL;

    QSpace<double,double> W;

    if ((nargin!=3 && nargin!=5) ||
        (nargin==3 && nargout>1) || (nargin==5 && nargout>2)) {
       if ((!nargin && !nargout) || (nargin && isHelpIndicator(argin[0]))) {
          usage(); return; }
       else wblog(FL,
         "ERR invalid number of I/O arguments (%d/%d)",nargin,nargout);
    }

    if (mxGetNumber(argin[2],dim)) wblog(FL,"ERR invalid dim (arg #3)");
    if (nargin==5) {
       try { mxIsQSpace(FL,argin[3],1); } catch (...) {
       wblog(FL,"ERR arg #4 must be valid real 1D QSpace"); }
       W.init(FL,argin[3]);
       if (mxGetNumber(argin[4],eps)) wblog(FL,"ERR invalid eps (arg #5)");
    }

    if (W.isEmpty() && nargout>1) wblog(FL,
       "ERR only one output argument (full transfer)");

    for (i=0; i<2; i++) {
       try { mxIsQSpace(FL,argin[i],r,'c'); } catch (...) {
       wblog(FL,"ERR arg #%d must be valid QSpace", i+1); }
    }
    if (dim<1 || dim>r) wblog(FL,"ERR invalid dim=%d",dim);

    if (mxIsQSpace(argin[0])>0 && mxIsQSpace(argin[1])>0) { 
       QSpace<double,double> A(argin[0]), B(argin[1]);

       A.TransferSpace(FL,B,dim,W,eps);

       if (W.isEmpty()) { a=B.toMx(); }
       else { a=A.toMx(); if (nargout>1) argout[1]=B.toMx(); }
    }
    else {
       QSpace<double,wbcomplex> A(argin[0]), B(argin[1]);
       A.TransferSpace(FL,B,dim,W,eps);

       if (W.isEmpty()) { a=B.toMx(); }
       else { a=A.toMx(); if (nargout>1) argout[1]=B.toMx(); }
    }

    argout[0]=a; 
}

