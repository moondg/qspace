
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to wbsvd.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "wbsvd"
#endif

#include "wblib.h"

template<class T>
int WB_SVD(
   wbarray<T> &A, wbvector<unsigned> &I,
   mxArray *argout[], int nargout
);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    wbvector<unsigned> I;

    MX_CHECK_HELPER_NARGS(1,-1,3); 
    if (nargin>2) usage(FL,"ERR invalid number of I/O arguments");

    if (nargin>1) { I.init(FL,argin[1]); }

    if (!mxIsComplex(argin[0])) {
       wbarray<double> A(argin[0]);
       WB_SVD(A,I,argout,nargout);
    }
    else {
       wbarray<wbcomplex> A(argin[0]);
       WB_SVD(A,I,argout,nargout);
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbsvd"); }
   aclu.Check();
};

template<class T>
int WB_SVD(
   wbarray<T> &A, wbvector<unsigned> &I,
   mxArray *argout[], int nargout
){
   wbarray<T> U,Vd,V;
   wbvector<double> S;
   wbperm P;

   if (I.len) {
      unsigned r=A.SIZE.len;
      if (I.len>r || I.anyLT(1) || I.anyGT(r)) wblog(FL,
         "ERR %s() got invalid index [%s]",myname,STR(I));
      I-=1; 
   }

   wbSVD(A,U,S,Vd,I);

   if (nargout<=1) argout[0]=S.toMx(); 
   else {
      argout[0]=U.toMx();
      argout[1]=S.toMx(); if (nargout>2) {
      argout[2]=V.toMx(); }
   }

   return 0;
};

