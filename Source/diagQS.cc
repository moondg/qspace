
/* CHANGE LOG ======================================================== *

   NB! not yet implemented for QSpace 3.0
   see QSpace/trace.m (!) // Wb,Apr10,15

 * =================================================================== */

char USAGE[]=""; // outsourced to diagQS.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "diagQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TQ, class TD>
mxArray* diagQS(const QSpace<TQ,TD> &A);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned r=-1; char sflag=0;

    MX_CHECK_HELPER_NARGS(1,-1,1); 

    if (nargin>2 || nargout>1) usage(FL,
       "ERR invalid number of I/O arguments");
    if (nargin==2) {
       if (Mx::IsEqual(argin[1],"-s")) sflag=1;
       else wblog(FL,"ERR %s() invalid usage",myname);
    }

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,double> A(argin[0],'r');
       wbMatrix<double> xd; A.diag(xd);
       if (sflag) { xd.SortRecs(); }
       argout[0]=xd.toMx(); 
    }
    else if (mxIsQSpace(0,0,argin[0],r,'c')>0) {
       const QSpace<gTQ,wbcomplex> A(argin[0]);
       wbMatrix<wbcomplex> xd; A.diag(xd);
       if (sflag) { xd.SortRecs(); }

       double nr=xd.normReal(), ni=xd.normImag();
       if (ni<nr && ni/nr<1E-12) {
          wbMatrix<double> xr; xd.getReal(xr);
          argout[0]=xr.toMx(); 
       }
       else {
          argout[0]=xd.toMx(); 
       }
    }
    else wblog(FL,"ERR %s() invalid input QSpace",myname);

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in diagQS"); }
   aclu.Check();
};

