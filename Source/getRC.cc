
/* CHANGE LOG ======================================================== *

   formerly called: getRCData() // Wb,Nov06,16

 * =================================================================== */

char USAGE[]=""; // outsourced to getRC.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "getRC"
#endif

#define LD_CLEBSCH_QS
#include "wblib.h"

template<class TD>
mxArray* load_CData(
   wbMatrix< CData<gTQ,TD> > &C, const mxArray *a0, char full=0);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    mxArray *a=0;

    MX_CHECK_HELPER_NARGS(1,-1,-1);

    if (CG::isCRef(argin[0],0)>0) {
       char full=0;

       if (nargin==2 && mxIsChar(argin[1])) { char s[8];
          if (!mxGetString(argin[1],s,8)) { 
             if (!strcmp(s,"-f")) { full=1; --nargin; } else
             if (!strcmp(s,"-F")) { full=2; --nargin; }
          }
       }

       if (nargin>1 || nargout>1) usage(FL,
          "ERR invalid number of I/O arguments");

       if (full<=1) {
          wbMatrix< CData<gTQ,double> > C;
          a=load_CData(C,argin[0],full);
       }
       else {
          wbMatrix< CData<gTQ,RTD> > C;
          a=load_CData(C,argin[0],full);
       }
    }
    else if (CG::isQSet(argin[0],0)>0) {
       QSet<gTQ> Q(FL,argin[0]); CData<gTQ,RTD> C(Q); CRef<gTQ> R(C);
       a=R.toMx();
    }
    else if (mxIsChar(argin[0])) {
       if (nargin==1) {
          mxGetString(argin[0],str,128);
          if (!strcmp(str,"--empty")) {
             CRef<gTQ> R; a=R.toMx();
          }
          else {
             QSet<gTQ> Q(FL,str); CData<gTQ,RTD> C(Q); CRef<gTQ> R(C);
             a=R.toMx();
          }
       }
       else if (mxIsChar(argin[1])) {
          wbstring ts(argin[0]);
          QType t(ts.data);

          genRG_struct<gTQ,RTD> &B=gRS.buf[t];
          B.checkInit(FL,t);

          if (nargin<2 || nargout>2) usage(FL,
             "ERR invalid number of I/O arguments");

          if (Mx::IsEqual(argin[1],"--ping")) { return; } 
          if (Mx::IsEqual(argin[1],"--info")) { 
             a=gRS.toMx(); if (nargout>1) {
             argout[1]=gCS.toMx(); }
          }
          else { 
             qset<gTQ> qs(FL,argin[1]);
             genRG_struct<gTQ,RTD> &B=gRS.Buf(t);
             const genRG_base<gTQ,RTD> &R=B.RSet[qs];

             a=R.toMx();
          }
       }
    }
    if (!a)
         { usage(); wblog(FL,"ERR %s() invalid usage",FCT); }
    else { argout[0]=a; }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getRC"); }
   aclu.Check();
};

template<class TD>
mxArray* load_CData(
   wbMatrix< CData<gTQ,TD> > &C, const mxArray *a0, char full) {

   unsigned i=0, m=mxGetM(a0), n=mxGetN(a0), nm=n*m;
   CRef<gTQ> R;

   C.init(m,n);

   for (; i<nm; ++i) {
      R.init(FL,a0,i);
      C[i].init(R); 
      C[i].cstat=R.cgb->cstat;
      if (full) {
         if (R.isRefInit()) { R.LoadRef(FL); }
         C[i].cgd.init(R.cgb->cgd);
      }
   }
   return C.toMx();
};

