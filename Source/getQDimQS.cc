
/* CHANGE LOG ======================================================== *

 * =================================================================== */

char USAGE[]=""; // outsourced to getQDimQS.m // Wb,Jan12,19

#define PROG_TAG "qdm"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getQDimQS"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template <class TQ, class TD>
mxArray* get_Qinfo(const QSpace<TQ,TD> &A) { 
   wbvector<unsigned> dd; A.qtype.Qlen(dd);
   return MXPut(0,0)
      .add(A.itags,"itags")
      .add(QDir(A.itags),"qdir")
      .add(A.qtype,"qtype")
      .add(dd,"qlen")
      .add(A.gotCGS(),"cgflag")
   .toMx();
};

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned k=0, r=-1; char isop=0; int i=0;

    wbMatrix<gTQ> Q;
    wbvector<widx_t> dd;
    wbMatrix<widx_t> dc; str[0]=0;

    MX_CHECK_HELPER_NARGS(2,2,4); 

    try { i=mxIsQSpace(FL,argin[0],r,'c'); }
    catch (...) { wblog(FL,"ERR invalid QSpace argument"); }
    if (i<=0) { wblog(FL,"ERR invalid QSpace argument"); }

    if (mxGetNumber(argin[1],k,'q')==0) {
       if (k==0 || (int)k>(int)r) wblog(FL,
          "ERR 2nd argument out of bounds (%d/%d)",k,r);
       k--; 
    }
    else {
       if (!mxIsChar(argin[1])) wbdie(FL,str);
       if (mxGetString(argin[1],str,12))
       wblog(FL,"ERR failed to read string (arg #2) ???");

       if (strcmp(str,"op")) wbdie(FL,"invalid 2nd argument");
       isop=1;
    }

    if (mxIsQSpace(argin[0])>0) { 
       const QSpace<gTQ,double> A(argin[0],'r',0);
       if (isop)     
            A.getQDim(  Q,dd,&dc);
       else A.getQDim(k,Q,dd,&dc);

       if (nargout>2 && !dc.data) dc.init2val(dd.len,A.nsym(), 1);
       if (nargout>3) argout[3]=get_Qinfo(A);
    }
    else {
       const QSpace<gTQ,wbcomplex> A(argin[0],'r',0);
       if (isop)     
            A.getQDim(  Q,dd,&dc);
       else A.getQDim(k,Q,dd,&dc);

       if (nargout>2 && !dc.data) dc.init2val(dd.len,A.nsym(), 1);
       if (nargout>3) argout[3]=get_Qinfo(A);
    }

    argout[0]=Q.toMx(); 

    if (nargout>1) argout[1]=dd.toMx_base_d();
    if (nargout>2) argout[2]=dc.toMx_base_d();

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getQDimQS"); }
   aclu.Check();
};

