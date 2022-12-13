/*------------------------------------------------------------------- *
 * mexworld.cc                                                        *
 * substiture *.mex[glx|a64] files in git-repository with *this       *
 * as a friendly reminder to compile locally                          *

   NB! rather use separate m-files for help info
    * this is easier to change the help without recompiling the binary!
    * if kept with the mex-file, this indicates in git repos
      in which path to the put the mex file!
    * easier to search for other people

   help mexworld will not accept / look into the mex-file itself.
   hence it looks on its path for mexworld.m; if it finds it, it prints
   the header help section; any commands in it are ignored.

   NB! for this work, the *.m file needs to be in the same
   directory as the mex-file, or *LATER* on the matlab path;
   so if the command is called, it finds the mex-file first.
   Only if it looks for the help, it looks for m-file.

   Wb,Jan12,19
 * ------------------------------------------------------------------ */

char USAGE[]="";

/*
>> help mexworld % prints header help in mexworld.m

>> mexworld -?   % since USAGE is empty, and mexworld.m exists
                 % this calls 'help mexworld' in matlab!
                 % i.e. becomes EQUIVALENT to 'help mexworld'!
                 % Wb,Jan12,19
*/

#include <string.h>
#include <mex.h>

int isHelpIndicator(const mxArray *a) {

   char istr[8];

   if (!a || mxGetString(a,istr,8)) { return 0; }

   if (istr[0]=='-') {
      if (!strcmp(istr, "-?") || !strcmp(istr, "-h")) return 1;
      if (!strcmp(istr, "--help")) return 1;
   }
   return 0;
}

void usage(const char *F=0, int L=0, const char* estr=0) {
    if (!USAGE[0] && (!estr || !estr[0])) {
       char cmd[256], f[64]; snprintf(f,64,"%s",mexFunctionName());
       if (f[0]) { snprintf(cmd,256,
            "if exist('%s.m','file')==2, help %s; end",f,f);
          mexPrintf("\n"); 
          mexEvalString(cmd); return;
       }
    }

    printf("\n%s\n", USAGE);
    if (estr && estr[0]) mexErrMsgIdAndTxt("Wb:MEX:wblog",estr);
};

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   if (nargin && isHelpIndicator(argin[0])) { usage(__FILE__,__LINE__); return; }

   mexPrintf("\n"
   "   Hello world!\n\n"
   "   This is a friendly reminder via %s()\n"
   "   that mex-files need to be recompiled locally.\n"
   "   Talk to A. Weichselbaum for details.\n\n"
   "   AW (C) 2019\n\n",mexFunctionName()
   );
};

