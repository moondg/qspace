
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to getpid.m // Wb,Feb14,19

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    int pid;

    MX_CHECK_HELPER_NARGS(0,-1,1);
    if (nargin) { ExitMsg("no input argument expected"); }

    pid=getpid();

    argout[0]=numtoMx(pid); 

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getpid"); }
   aclu.Check();
};

