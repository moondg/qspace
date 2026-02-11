
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to plusQS.m // Wb,Jan12,19

#define PROG_TAG "pls"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "plusQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TQ, class TD>
inline QSpace<TQ,TD>& MPS_PLUS_QS(
   QSpace<TQ,TD> &C, const mxArray *argin[], int nargin
){
   double bfac=1;

   unsigned l=nargin-1;
   char vflag=0, xflag=0, str[8]; str[0]=0;

   for (; l>=2 && mxIsChar(argin[l]); --l) {
      mxGetString(argin[l],str,8);
      if (!strcmp(str,"-v")) { ++vflag; } else
      if (!strcmp(str,"-x")) { ++xflag; }
      else wblog(FL,"ERR %s() invalid input option '%s' [%d]",FCT,str,l+1);
   }

   if (l>2) wblog(FL,"ERR invalid usage (nargin=%d)",nargin);
   if (l==2) {
      if (mxGetNumber(argin[l], bfac))
      wblog(FL,"ERR %s() invalid option '%s'",FCT,str);
   }

   const QSpace<TQ,TD> A(FL,argin[0],'r',1, xflag? 0:1);
   const QSpace<TQ,TD> B(FL,argin[1],'r',1, xflag? 0:1);

   if (A.gotCGS(FL)<2) {
      C.mt=Wb::MEX_RETURN; 
      A.plus_plain(B,C,bfac); 
   }
   else {
      C.Cat(FL,A,B,TD(1),TD(bfac),vflag? vflag:1);
   }

   C.SkipZeroData();   
   C.NormCGW();

   C.SetFDir(FL,A,B);  
   C.ctime=Wb::getTimeNow();

   return C;
};

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    char isra, isrb; 

    MX_CHECK_HELPER_NARGS(2,-1,1); 
    if (nargin>4) wblog(FL,"ERR invalid usage (nargin=%d)",nargin);

    isra=(mxIsQSpace(argin[0])>0); 
    isrb=(mxIsQSpace(argin[1])>0);

    if (isra && isrb) {
       QSpace<gTQ,double> C;
       MPS_PLUS_QS(C,argin,nargin); 
       argout[0]=C.toMx();
    }
    else {
       QSpace<gTQ,wbcomplex> C;
       MPS_PLUS_QS(C,argin,nargin); 
       argout[0]=C.toMx();
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in plusQS"); }
   aclu.Check();
};

