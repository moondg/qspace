
#include <mex.h>

void mexFunction(
   int nargin, mxArray *argin[], int nargout, const mxArray *argout[]
){
   printf("\n   hello world!  (%d/%d)\n\n",nargout,nargin);
}

