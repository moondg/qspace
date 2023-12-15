
#include <mex.h>

// using helloworld without wblib.h
// to include wblib.h -> see mexworld.cc // Wb,Jun27,23

#define PP_STR__(a) #a
#define PP_STRFY(a) PP_STR__(a)

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargout) { mxArray *a;
      for (unsigned i=0; i<nargout; ++i) {
         if (i<nargin) { a=mxDuplicateArray(argin[i]); }
         else {
            a=mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetDoubles(a)[0]=double(i+1);
         }
         argout[i]=a;
      }
   }

   printf("\n   Hello world from MEX! (.%s: nargout=%d, nargin=%d)\n\n",
      PP_STRFY(MEX_EXT), 
      nargout, nargin
   );
}

