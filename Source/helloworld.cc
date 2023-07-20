
#include <mex.h>

// using helloworld without wblib.h
// to include wblib.h -> see mexworld.cc // Wb,Jun27,23

#define PP_STR__(a) #a
#define PP_STRFY(a) PP_STR__(a)

void mexFunction(
   int nargin, mxArray *argin[], int nargout, const mxArray *argout[]
){
   printf(
   "\n   hello world!\n"
   "\n   This is %s() [.%s] (nargout=%d, nargin=%d)\n\n",
   mexFunctionName(), PP_STRFY(MEX_EXT), nargout, nargin);
}

