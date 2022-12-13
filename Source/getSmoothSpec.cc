
/* CHANGE LOG ======================================================== *

 * =================================================================== */

char USAGE[]=""; // outsourced to getSmoothSpec.m // Wb,Jan12,19

#ifdef MATLAB_MEX_FILE 
   #define PROG mexFunctionName()
#else
   #define PROG "getSmoothSpec"
#endif

#include "wblib.h"
#include "spectral.hh"

template <class TS>
void GET_SMOOTH_SPEC(
   int nargout, mxArray** &argout,
   int nargin, const mxArray** &argin, Spectral<TS> &ASpec);

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(2,-1,-1); 

   if (nargout<2) usage(FL,"invalid number of I/O arguments");

   if (!Mx::IsDblVector(FL,argin[0])) wberror(FL,
	  "invalid input arg#1 (real om data required)");

   if (Mx::IsDblMat(0,0,argin[1])) {
	  Spectral<double> ASpec;
	  GET_SMOOTH_SPEC(nargout, argout, nargin, argin,ASpec);
   }
   else {
	  if (!Mx::IsDblMat(FL,argin[1],'C')) wberror(FL,
		 "invalid input arg#2 (double or complex required)");
	  Spectral<wbcomplex> ASpec;
	  GET_SMOOTH_SPEC(nargout, argout, nargin, argin,ASpec);
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in getSmoothSpec"); }
};

template <class TS>
void GET_SMOOTH_SPEC(
   int nargout, mxArray** &argout,
   int nargin, const mxArray** &argin, Spectral<TS> &ASpec
) {
    unsigned i,j,n, N, vflag=1, isfunc, a0flag=1;

    wbvector<double> om;
    wbarray<double> Om;

    wbMatrix<TS> aa;
    wbarray<TS> Ar;

    wbvector<unsigned> D;
    OPTS opts;

    mxArray *a;

    om.init(FL,argin[0]);
    aa.init0(argin[1]);

    ASpec.Emin=om.aMin(1);
    ASpec.eps=ASpec.Emin*1E-3;

    opts.init(argin+2, nargin-2);

    isfunc=opts.getOpt("func");

    if (!opts.getOpt("eps",ASpec.eps)) wblog(FL, 
       "ERR %s() option 'eps' required (pick ~T/2)",FCT);

    opts.getOpt("sigma", ASpec.sigma);
    opts.getOpt("sigma2",ASpec.sigma2);
    opts.getOpt("alpha", ASpec.alpha); 

    opts.getOpt("nlog",  ASpec.nlog);
    opts.getOpt("emin",  ASpec.Emin);
    opts.getOpt("emax",  ASpec.Emax);

    opts.getOpt("om", a);
    if (opts.getOpt("-a0")) a0flag=0;

    if (opts.getOpt("-q")) { vflag=0; } else 
    if (opts.getOpt("-v")) { vflag=2; } 
    else { opts.getOpt("vflag",vflag); }

    opts.checkAnyLeft(FL); 

    N=aa.dim2;

    if (om.len!=N) wblog(FL,
       "ERR %s() size mismatch (a0 vs. om: %d/%d !?)",myname,om.len,N);
    if (vflag) {
       if (aa.dim1==1)
            { sprintf(str,"one data set"); }
       else { sprintf(str,"%ld data sets",aa.dim1); }
       wblog(FL,"<i> %s() %s @ %d pts",myname,str,N);
    }

    for (n=i=0; i<aa.dim1; ++i)
    for (  j=0; j<aa.dim2; ++j) { if (ISNAN(aa(i,j))) { aa(i,j)=0; n++; }}

    if (n) wblog(FL, "    skipping %d NaN's", n);

    if (isfunc) {
       wbvector<double> dom(N); unsigned m;
       if (vflag) wblog(FL,
          "--> transforming 'func' to binned data");

       for (i=1; i<N; i++) if (om[i]<=om[i-1]) wblog(FL,
          "ERR option `func' requires increasing set of omega's");

       dom[0]=om[1]-om[0];
       for (m=N-1, i=1; i<m; i++) dom[i]=0.5*(om[i+1]-om[i-1]);
       dom[i]=om[i]-om[i-1];

       for (i=0; i<aa.dim1; ++i)
       for (j=0; j<aa.dim2; ++j) aa(i,j)*=dom[j];
    }

    ASpec.init("A",aa.dim1);

    if (a0flag) {
       wbindex I; wbMatrix<TS> ax;
       wbvector<TS> a2;
       unsigned gotdelta=0;
       widx_t i; 

       wbMatrix<TS> a0(aa.dim1,2); 

       om.maxneg(&i);
       if (int(i)>=0) {
          om.findRange(2*om[i],om[i],I,"[["); 
          if (!I.isEmpty()) {
             aa.getCols(I,ax).recSumA(a2); a2.setMin(1E-14);
             for (j=0; j<a2.len; ++j) {
                if (ABS2(aa(j,i))>ABS2(a2[j])) {
                   a0(j,0)+=aa(j,i); aa(j,i)=0; ++gotdelta;
                }
             }
          }
       }

       om.minpos(&i);
       if (int(i)>=0) {
          om.findRange(om[i],2*om[i],I,"]]"); 
          if (!I.isEmpty()) {
             aa.getCols(I,ax).recSumA(a2); a2.setMin(1E-14);
             for (j=0; j<a2.len; ++j) {
                if (ABS2(aa(j,i))>ABS2(a2[j])) {
                    a0(j,1)+=aa(j,i); aa(j,i)=0; ++gotdelta;
                }
             }
          }
       }

       if (gotdelta) {
          if (vflag) {
             const char *istr="subtracted delta(0) =";
             wbstring as=a0.toStr("%.3g"); unsigned l=as.strlen_();
             if (l<30) wblog(FL,"NB! %s [%s]",istr,as.data);
             else { wblog(FL,"NB! %s",istr);
                if (l<50)
                     { wblog(FL,"--> [%s]",as.data); }
                else { printf("\n\n[ %s ]\n\n",as.data); }
             }
          }

          a0.save2(ASpec.A0); 
       }
    }

    D.init(2); D[0]=1; D[1]=N;
    Om.init2ref(om.data,D);

    for (i=0; i<aa.dim1; ++i) {
       Ar.init2ref(aa.rec(i),D);
       ASpec.Add(Om,Ar,i);
    }

    if (a) {
         om.init(FL,a);
         ASpec.getSmoothSpec(aa,om,vflag); } 
    else ASpec.getSmoothSpec(om,aa,vflag);   

    argout[0]=om.toMx('t');
    argout[1]=aa.toMx('r');

    if (nargout>2) { MXPut S(FL); 
       S.add(ASpec.Emin,"emin").add(ASpec.Emax,"emax")
        .add(ASpec.nlog,"nlog").add(ASpec.sigma,"sigma")
        .add(ASpec.sigma2,"sigma2").add(ASpec.alpha,"alpha")
        .add(ASpec.eps,"eps").add(ASpec.A0,"a0");
       S.save2(argout[2]);
    }

    if (vflag) printf("\n");
};

