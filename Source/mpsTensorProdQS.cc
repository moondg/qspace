/* ==================================================================
 * tstmex.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]="\
   C = mpsTensorProdQS(A,B)                                       \n\
                                                                  \n\
      Tensor product of two equal-rank QSpaces with A being       \n\
      the fast index (on the data level!) consistent with         \n\
      MatLab's column-major convention.                           \n\
                                                                  \n\
   Wb,Aug10,06                                                    \n\
";

#include "wblib.h"

template<class TQ, class TD>
inline mxArray* MPS_TENSOR_PROD_QS(
   QSpace<TQ,TD> &C, const mxArray *argin[], char ref=0
){
   const QSpace<TQ,TD> A(argin[0],ref);
   const QSpace<TQ,TD> B(argin[1],ref);

   A.TensorProd(B,C); return C.toMx();
};

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
    char isra, isrb; 

    if (nargin==0 || checkHelpVersion(argin[0])) { return; }
    if (nargin!=2) wblog(FL,
    "ERR mpsTensorProdQS requires two input QSpaces A and B.");

    isra=(mxIsQSpace(argin[0])>0); 
    isrb=(mxIsQSpace(argin[1])>0);

    if (isra && isrb) {
       QSpace<double,double> C;
       argout[0]=MPS_TENSOR_PROD_QS(C,argin,'r'); 
    }
    else {
       QSpace<double,wbcomplex> C;
       argout[0]=MPS_TENSOR_PROD_QS(C,argin); 
    }
};

