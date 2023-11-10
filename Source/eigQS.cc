
/* CHANGE LOG ======================================================== *

   formerly mpsEigenSymmQS() // Wb,Aug24,15

 * =================================================================== */

char USAGE[]=""; // outsourced to eigQS.m // Wb,Jan12,19

#define PROG_TAG "eig"

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "eigQS"
#endif

#define LOAD_CGC_QSPACE
#include "wblib.h"

template<class TA>
void EIGEN_SYM_QS(
    const char *F, int L, QSpace<gTQ,TA> &A,
    int nargout, mxArray *argout[], int nargin, const mxArray *argin[]
);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    MX_CHECK_HELPER_NARGS(1,-1,-1);

    char isra=(mxIsQSpace(argin[0])>0);
    if (!isra) { mxIsQSpace(FL,argin[0],'c'); }

    if (isra) {
       QSpace<gTQ,double> A(argin[0]); 
       EIGEN_SYM_QS(FL,A,nargout,argout,nargin-1,argin+1);
    }
    else {
       QSpace<gTQ,wbcomplex> A(argin[0]);
       EIGEN_SYM_QS(FL,A,nargout,argout,nargin-1,argin+1);
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in eigQS"); }
   aclu.Check();
};

template<class TA>
void EIGEN_SYM_QS(
    const char *F, int L, QSpace<gTQ,TA> &A,
    int nargout, mxArray *argout[], int nargin, const mxArray *argin[]
){
    int Nk,Nkeep=-1; char vflag=0, Rflag=0, mKD=0; unsigned r;
    double Etrunc=0; 
    double deps=0;   
    char cgflag=(A.gotCGS(F_L)>0);

    wbMatrix<double> Etot;
    wbMatrix<unsigned> DD;

    QSpace<gTQ,double> Ek,Et;
    QSpace<gTQ,TA> Ak,At;

    if ((r=A.rank(FL))%2) wblog(FL,"ERR invalid rank-%d QSpace",r);
    A.skipZeroOffDiag(1E-14);

    if (nargin) {
       OPTS opts;
       opts.init(argin, nargin);

       opts.getOpt("Nkeep",Nkeep);
       opts.getOpt("deps",deps);

       if (!opts.getOpt("Etrunc",Etrunc)) {
       if ( opts.getOpt("Rtrunc",Etrunc)) Rflag=1; }

       vflag=opts.getOpt("-v");

       mKD=opts.getOpt("--mKD"); 

       opts.checkAnyLeft(); 
    }

    Nk=Nkeep;
    if (!Rflag) {
       if (deps) { 
          A.EigenSymmetric(Ak,At,Ek,Et,Etot,DD,Nk,Etrunc,NULL,mKD,
          wbperm(),deps,-1.,-1);
       }
       else {
          A.EigenSymmetric(Ak,At,Ek,Et,Etot,DD,Nk,Etrunc,NULL,mKD);
       }
    }
    else { A.EigenSymmetric(Ak,At,Ek,Et,Etot,DD,Nk,Etrunc,NULL,mKD,
         wbperm(),0,0,-1,"desc"); 
    }

    if (vflag && Nkeep>=0) wblog(FL,
       "<i> kept %d/%d %s", Nk,Etot.dim1, cgflag ? "multiplets" : "states");

    argout[0]=Etot.toMx(); 

    if (nargout>1) {
       mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

       mxAddField2Scalar(FL,S,"AK",Ak.toMx());
       mxAddField2Scalar(FL,S,"AD",At.toMx());
       mxAddField2Scalar(FL,S,"EK",Ek.toMx());
       mxAddField2Scalar(FL,S,"ED",Et.toMx());

       mxAddField2Scalar(FL,S,"DB",DD.toMx());

       mxAddField2Scalar(FL,S,"NK",numtoMx(Nk));
       if (Rflag==0)
            mxAddField2Scalar(FL,S,"Etrunc",numtoMx(Etrunc));
       else mxAddField2Scalar(FL,S,"Rtrunc",numtoMx(Etrunc));
       argout[1]=S;
    }
};

