/*------------------------------------------------------------------- *
 * mexworld.cc                                                        *
 * substitute *.mex[glx|a64] files in git-repository with *this       *
 * as a friendly reminder to compile locally                          *

   NB! rather use separate m-files for help info
    * this is easier to change the help without recompiling the binary!
    * if kept with the mex-file, this indicates in git repos
      in which path to the put the mex file!
    * easier to search for other people

  `help mexworld' will not accept / cannot run the mex-file itself.
   hence it looks on its path for mexworld.m; if it finds it,
   it prints the header help section; any commands in it are ignored.

   NB! for this work, the *.m file needs to be in the same
   directory as the mex-file, or *LATER* on the matlab path;
   so if the command is called, it finds the mex-file first.
   Only if it looks for the help, it looks for m-file.

   Wb,Jan12,19
 * ------------------------------------------------------------------ */

char USAGE[]=""; // outsourced to wbhist.m // Wb,Feb14,19

/* see Archive/mexworld_230627.cc
// for old setup without wblib.h testing `help mexworld' etc.

   prints header help in mexworld.m
   >> help mexworld

   >> mexworld -?
      since USAGE is empty, and mexworld.m exists
      this calls 'help mexworld' in matlab!
      i.e. becomes EQUIVALENT to 'help mexworld'! // Wb,Jan12,19
*/

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "mexworld"
#endif

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(0,-1,0);

   if (!nargin) { mexPrintf("\n"
   "   Hello world, this is the %s() [.%s] !\n\n"
   "   This C++ mex-file needs to be (re)compiled locally.\n"
   "   It is complimentary to %s.m which contains the usage information (help) only.\n\n"
   "   Type `help %s', `%s -h', or `%s --version'\n"
   "   for typical help (usage, documentation, and compilation info) of mex files.\n\n"
   "   AW (2019)\n\n",myname,PP_STRFY(MEX_EXT),myname,myname,myname,myname);
       return;
   }

   for (int i=0; i<nargin; ++i) {
      if (!mxIsChar(argin[i])) { wblog(FL,
         "WRN ignoring input argument #%d (%s) expecting strings",
          i+1,mxGetClassName(argin[i]));
         continue;
      }
      if (mxGetString(argin[i],str,24)) { wblog(FL,
         "WRN invalid option %d (string out of bounds)",i+1);
         continue;
      }
      if (!strcmp(str,"--check-hl")) { if (i) { PRINTF("\n"); }
         wblog(FL,"TST useCol = %d",Wb::useCol);
         wblog(FL,"NB! some message");
         wblog(FL,"WRN some message");
         wblog(FL,"ENV some message");
         wblog(FL,"ok. some message");
         wblog(FL,"*   some message");
         wblog(FL," *  some message");
         wblog(FL,"  * some message");
         wblog(FL,"hdr some message");
         wblog(FL,"DBG some message");
         wblog(FL,"XXE some message");
         wblog(FL,"XXW some message");
         wblog(FL,"XXG some message");
         wblog(FL,"XXY some message");
         wblog(FL,"XXB some message");
      }
      else if (!strcmp(str,"--check-env")) { if (i) { PRINTF("\n"); }
         wblog(FL," *  envVRB = %d = %s",Wb::envVRB,BITS(Wb::envVRB));
         wblog(FL," *  envDKT = %d",Wb::envDKT);
         wblog(FL," *  useCol = %d",Wb::useCol);
         wblog(FL," *  envDBG = %d",Wb::envDBG);
      }
      else { wblog(FL,"TST ignoring input '%s'",str); }
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbhist"); }
   aclu.Check();
};

