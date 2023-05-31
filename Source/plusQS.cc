
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to plusQS.m // Wb,Jan12,19

#define PROG_TAG "pls"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "plusQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

template<class TQ, class TD>
inline QSpace<TQ,TD>& MPS_PLUS_QS(
   QSpace<TQ,TD> &C, const mxArray *argin[], int nargin
){
   const QSpace<TQ,TD> A(FL,argin[0],'r');
   const QSpace<TQ,TD> B(FL,argin[1],'r');

   double bfac=1; str[0]=0;
   unsigned l=nargin-1; char vflag=0;
   if (l>=2 && mxIsChar(argin[l])) { str[0]=0;
      mxGetString(argin[l],str,128);
      if (!strcmp(str,"-v")) { vflag='v'; --l; --nargin; }
      else wblog(FL,
         "ERR %s() invalid input option '%s' [%d]",str,FCT,l+1
      );
   }
   if (nargin>2) if (mxGetNumber(argin[2], bfac)) wblog(FL,str);

   if (A.gotCGS(FL)<2) {
      C.mt=Wb::MEX_RETURN; 
      A.plus_plain(B,C,bfac); 
   }
   else {
      C.Cat(FL,A,B,TD(1),TD(bfac),vflag? vflag:1);
   }

   C.SkipZeroData();   

   C.NormCGW();
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

