
/* COMMENTS / CHANGE LOG ============================================= *
   e.g. that also allows to use backspace on command line output
   for i=1:999, printfc([zeros(1,10)+8, sprintf('i=%03d',i)]); end
 * =================================================================== */

char USAGE[]=""; // outsourced to getchar.m // Wb,Feb14,19

#include "wblib.h"

void got_alarm(int sig) {
   fprintf(stderr,"Got signal %d\n",sig);
}

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin,  const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    wbvector<int> cc(16);
    unsigned i=0;

#if 0
    Wb::SigHandler SIG(FL);
    alarm(1); 
    signal(SIGALRM,got_alarm);
#endif

    MX_CHECK_HELPER_NARGS(0,-1,1); 
    if (nargin) usage(FL,"no input arguments expected");

    cc[i++]=getc(stdin);

    if (cc[0]==27) { 

       cc[i++]=getc(stdin);
       cc[i++]=getc(stdin);

       if (cc[1]=='[') {
          while (1) {
             cc[i++]=getc(stdin);
             if (cc[i-1]<0 || cc[i-1]==126) break;
          }
       }
    }

    cc.len=i;
    argout[0]=cc.toMx(); 

    return;

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getchar"); }
   aclu.Check();
};

