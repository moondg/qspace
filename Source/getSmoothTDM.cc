
/* CHANGE LOG ======================================================== *
 * =================================================================== */

char USAGE[]=""; // outsourced to getSmoothTDM.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE 
   #define PROG mexFunctionName()
#else
   #define PROG "getSmoothTDM"
#endif

#include "wblib.h"
#include "spectral.hh" 

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

    unsigned i,k,K,M,N, disc, func, vflag=1, r=0, nlog=100; int e;
    double eps=1E-8, sigma=-1, alpha=0., Lambda=-1, emin=1E-8, emax=10.;
    TDSpectral<double> tdmData;
    char rawflag=0, istr[128]; istr[0]=0;

    wbvector<double> tt,om;
    wbvector<unsigned> D;
    wbMatrix<wbcomplex> az;
    wbMatrix<double> aa,a0;
    OPTS opts;

    wbarray<double> Om,Ar,AR;

    MX_CHECK_HELPER_NARGS(1,-1,-1); 

    if (!Mx::IsDblVector(FL,argin[0]))
    wberror(FL,"invalid first argument (real data required).");
    tt.init(FL,argin[0]);
    if (!tt.len) wblog(FL,"ERR No time data specified.");

    if (mxIsStruct(argin[1])) {
       k=2; if (nargout>2 || nargin<(int)k)
       usage(FL,"invalid number of I/O arguments.");

       if ((e=tdmData.init(argin[1],emin,emax,nlog)))
       wblog(FL,"ERR invalid input RAW data (%d).",e);

       tdmData.tt=tt;
    }
    else {
       k=3; if (nargout>2 || nargin<(int)k)
       usage(FL,"invalid number of I/O arguments.");

       if (!Mx::IsDblVector(0,__LINE__,argin[0]) ||
           !Mx::IsDblVector(0,__LINE__,argin[1]) ||
            Mx::IsNumArray(0,__LINE__,argin[2],-1)<=0 ||
           (r=mxGetNumberOfDimensions(argin[2]))>3 || r<2)
       wberror(FL,"invalid input arguments (%s)",str);

       om.init(FL,argin[1]);
       AR.init((mxArray*)argin[2]);

       if (om.len!=AR.SIZE[0]) wblog(FL,
          "ERR omega does not match Araw (%d; %d)", om.len, AR.SIZE[0]);

       rawflag=1;
    }

    opts.init(argin+k, nargin-k);

    if (!rawflag) {
       double dbl=-1; unsigned u=0; i=0;
       i+=opts.getOpt("emin",dbl);
       i+=opts.getOpt("emax",dbl);
       i+=opts.getOpt("nlog",u  );
       if (i) wblog(FL,"WRN emin/emax/nlog ignored on input for alpha>0");
    }
    else {
       opts.getOpt("emin",emin); eps=emin*10;
       opts.getOpt("emax",emax);
       opts.getOpt("nlog",nlog);
    }

    opts.getOpt("alpha", alpha);
    opts.getOpt("sigma", sigma);
    if (sigma<=0) {
       if (alpha>0)
            opts.getOpt("Lambda",Lambda);
       else sigma=0.1; 
    }

    if (!opts.getOpt("eps",eps)) wblog(FL, 
       "ERR %s() option 'eps' required (pick ~T/2)",FCT);

    disc=opts.getOpt("disc");
    func=opts.getOpt("func");

    if (disc && func) wblog(FL,
       "ERR specify either 'disc' or 'func' (%d/%d)",disc,func);

    if (opts.getOpt("nolog") || opts.getOpt("-q")) { vflag=0; } else
    if (opts.getOpt("-v")) { vflag+=1; }

    opts.checkAnyLeft(); 

    if (rawflag) {
       if (AR.SIZE.len<3) AR.appendSingletons(3);
       N=AR.dim(1); M=AR.dim(2); K=AR.dim(3);

       if (K==1) { if (vflag) wblog(FL,
          "<i> %s() %d data set(s) with %d pts",myname,M,N); }
       else if (AR.SIZE.len==3) wblog(FL,
          "<i> %s() %d data set(s) with %d pts (x%d)",myname,M,N,K);
       else wblog(FL,"ERR invalid Araw data (size %s)",SSTR(AR));

       Wb::replRange(FL,AR.data,AR.numel(),double(NAN),0.); 

       tdmData.init(FL,"TDM",tt,K,M, emin,emax,nlog,vflag);

       D.init(2); D[0]=1; D[1]=N;
       Om.init2ref(om.data,D);

       for (k=0; k<K; k++)
       for (i=0; i<M; i++) {
          Ar.init2ref(AR.ref(i,k),D);
          tdmData.Add(Om,Ar,k,i);
       }
    }

    tdmData.getSpecData(om,a0);

    if (sigma>0) {
       sprintf(istr,"frequency broadened TDM data (sigma=%.4g)",sigma);
       if (vflag) wblog(FL,"<i> %s",istr);

       tdmData.getSmoothSpec_t(sigma,eps,LOG_GAUSS_BRD,vflag);
       tdmData.Fourier(om,aa,tt,az,alpha,vflag); 
    }
    else if (alpha>0) {
       if (Lambda<=1) wblog(FL,"ERR alpha requires Lambda (%.4g)",Lambda);
       sprintf(istr,"time-domain broadened TDM data (alpha=%.4g)",alpha);
       if (vflag) wblog(FL,"<i> %s",istr);

       tdmData.Fourier(om,az,alpha,Lambda,vflag);
    }
    else {
       sprintf(istr,"plain Fourier transformed data (no broadening)");
       if (vflag) wblog(FL,"<i> %s",istr);

       Wb::Fourier(om,a0,tt,az, func? 'f': 0); 
    }
    if (vflag) printf("\n");

    argout[0] = az.toMx('r');

    if (nargout>1) { mxArray *S;
       S=mxCreateStructMatrix(1,1,0,NULL);

       mxAddField2Scalar(FL,S,"istr", wbstring(istr).toMx());
       mxAddField2Scalar(FL,S,"sigma",numtoMx(sigma));
       mxAddField2Scalar(FL,S,"eps",  numtoMx(eps));
       mxAddField2Scalar(FL,S,"alpha",numtoMx(alpha));  if (alpha>0) {
       mxAddField2Scalar(FL,S,"Lambda",numtoMx(Lambda)); }

       mxAddField2Scalar(FL,S,"om",   om.toMx('t'));
       mxAddField2Scalar(FL,S,"aa",   aa.toMx('r'));
       mxAddField2Scalar(FL,S,"a0",   a0.toMx('r'));

       mxAddField2Scalar(FL,S,"emin", numtoMx(emin));
       mxAddField2Scalar(FL,S,"emax", numtoMx(emax));
       mxAddField2Scalar(FL,S,"nlog", numtoMx(nlog));

       argout[1]=S;
    }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getSmoothTDM"); }
};

