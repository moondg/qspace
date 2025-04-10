
#include <string.h>
#include <mex.h>

// using helloworld without wblib.h
// to include wblib.h -> see mexworld.cc // Wb,Jun27,23

#define PP_STR__(a) #a
#define PP_STRFY(a) PP_STR__(a)

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   int i=0, l=0; char vflag=1;

   char mark[nargin]; memset(mark,0,nargin);
   for (; i<nargin; ++i) {
      if (mxIsChar(argin[i])) { char str[8]; str[0]=0;
         if (!mxGetString(argin[i],str,7)) { 
            if (!strcmp(str,"-q")) { mark[i]=(++l); vflag=0; }
         }
      }
   }

   if (nargout) { mxArray *a;
      for (i=0; i<nargout; ++i) {
         if (i<nargin && !mark[i]) {
            a=mxDuplicateArray(argin[i]);
         }
         else {
            a=mxCreateDoubleMatrix(1,1,mxREAL);
            mxGetDoubles(a)[0]=double(i+1);
         }
         argout[i]=a;
      }
   }

   if (vflag) { printf("\n   "
     "Hello world from MEX! (.%s: nargout=%d, nargin=%d)\n\n",
      PP_STRFY(MEX_EXT), 
      nargout, nargin);
   }
}

