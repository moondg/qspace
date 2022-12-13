char USAGE[]=
/* ================================================================== */
"[A, I] = mpsFullQS(A)              \n\
                                    \n\
   convert QSpace to full array.    \n\
   I provides additional info.      \n\
                                    \n\
Wb,Aug20,06                         \n\
";

// This is a MEX wrapper routine for MATLAB.
/* ================================================================ */

#include "wblib.h"

#define _TQ double
#define _TD wbcomplex

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){

    QSpace<_TQ,_TD> Ain,C;
    wbarray<_TD> A;
    wbvector< wbMatrix<_TQ> > QB;
    wbvector< wbvector<unsigned> > SB;

    if (nargin) if (checkHelpVersion(argin[0])) { return; }

    if (nargin!=1) wblog(FL,
    "ERR Need single input QSpace.");

    Ain.init(FL,argin[0]);

    Ain.toFull(A,QB,SB);
    argout[0]=A.toMx(); 

    if (nargout>1) {
       mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

       mxAddField2Scalar(FL,S, "QB", vecMat2Mx(QB));
       mxAddField2Scalar(FL,S, "SB", vecVec2Mx(SB));

       argout[1]=S;
    }
}

