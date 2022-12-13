char USAGE[] =
/* ==================================================================
 * mpsExpandQS.cc */

"   Usage: A_out [,Iold, mQ] = mpsExpandQS(A_in, QK, d, k [,nrm]) \n\
                                                                  \n\
       expand given qidx at index k with QK whos dimensions       \n\
       must be specified in vector d. the overall norm of newly   \n\
       added space can be determined by nrm (1E-12)               \n\
                                                                  \n\
   Output                                                         \n\
                                                                  \n\
      A_out  updated A_in                                         \n\
      Iold   index of where to find A_in in A_out                 \n\
      mQ     maker space whether given symmetry already existed in A_in\n\
                                                                  \n\
   AW (C) Dec 2006                                                 \n";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    wbvector<unsigned> d;
    wbMatrix<double> QK;
    wbvector<char> mQ, *mQp;
    wbindex kk;
    double nrm=1E-12;
    unsigned i,r;
    mxArray *a;
    wbindex I, *Ip;
    int rk=-1;

    if (nargin) if (checkHelpVersion(argin[0])) { return; }

    if (nargin<4 || nargin>5 || nargout>3)
    usage(FL,"ERR Invalid number of I/O arguments.");

    if (mxIsQSpace(argin[0],rk,'C')<=0) wblog(FL,
       "ERR invalid QSpace input");
    r=(unsigned)rk;

    QK.init(argin[1]);
    d .init(FL,argin[2]);
    kk.init(FL,argin[3]);

    if (nargin>4)
    mxGetNumber(argin[4],nrm);

    if (kk.len==0 || !kk.isUnique(r,1)) 
    wblog(FL,"ERR Invalid index set k !?");

    Ip  = nargout>1 ? &I  : NULL;
    mQp = nargout>2 ? &mQ : NULL;

    if (mxIsQSpace(argin[0])>0) { 
       QSpaceRM<double,double> A(argin[0],'r');
       for (i=0; i<kk.len; i++) A.ExpandQ(QK,d,r-kk[i],nrm, Ip, mQp);
       a=A.toMx('r'); 
    }
    else {
       QSpaceRM<double,wbcomplex> A(argin[0],'r');
       for (i=0; i<kk.len; i++) A.ExpandQ(QK,d,r-kk[i],nrm, Ip, mQp);
       a=A.toMx('r'); 
    }

    argout[0]=a; 

    if (nargout>1) { I+=1; argout[1]=I.toMx('r'); }
    if (nargout>2) argout[2]=mQ.toMx('r');
}

