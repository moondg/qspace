
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to skipZerosQS.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "skipZerosQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned k=1,r=-1; char all=0, isr;
    double eps=1E-14;

    MX_CHECK_HELPER_NARGS(1,-1,-1); 

    for (; int(k)<nargin; ++k) {
       if (mxIsNumeric(argin[k])) { int e; str[0]=0;
          try { e=mxGetNumber(argin[k],eps); }
          catch (...) { e=99; }
          if (e) wblog(FL,"ERR invalid arg #%d (eps, e=%d)\n%s",k+1,e,str);
       }
       else if (mxIsChar(argin[k])) { char s[8]; s[0]=0;
          if (mxGetString(argin[k],s,8)) wblog(FL,
             "ERR invalid arg #%d (string)",k+1); 
          if (!strcmp(s,"--all") || !strcmp(s,"-f")) { all=1; }
          else { 
             wblog(FL,"ERR invalid usage (arg #%d '%s')",k+1,s);
          }
       }
       else { 
         wblog(FL,"ERR invalid usage (arg #%d) !?",k+1);
       }
    }

    isr=(mxIsQSpace(argin[0])>0); str[0]=0;
    if (!isr) {
       mxIsQSpace(FL,argin[0],r,'c',-1,NULL,NULL,
       "invalid input QSpace arg #1");
    }

    if (isr) {
       QSpace<gTQ,double> A;
       A.init(FL,argin[0],'r'); k=A.SkipZeroData(eps,0,1,all);
       argout[0]=A.toMx(); 
    }
    else {
       QSpace<gTQ,wbcomplex> A;
       A.init(FL,argin[0]); k=A.SkipZeroData(eps,0,1,all);
       argout[0]=A.toMx(); 
    }

    if (nargout>1) argout[1]=numtoMx(k);

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in skipZerosQS"); }
   aclu.Check();
};

