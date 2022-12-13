
/* COMMENTS / CHANGE LOG ============================================= *
   adapted from Wb::FixRational()
 * =================================================================== */

char USAGE[]=""; // outsourced to wbrat.m // Wb,Feb14,19

#define PROG_TAG "rat"
#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(1,-1,2); 
   if (!Mx::IsDblMat(FL,argin[0])) usage(FL,"invalid usage");

   wbMatrix<double> A(FL,argin[0]);

   char rflag=1, aflag=0, vflag=1, sflag=1;
   unsigned i, n2, d1=A.dim1, d2=A.dim2, n=d1*d2, nmax=6;
   int e1, e2, nx; long p,q, p2,q2;
   double a,x,Ai,a2,x2, err2, r2=0, eps1=1E-6, eps2=1E-12;
   wbMatrix<double> P(d1,d2), Q(d1,d2), err(d1,d2);
   wbMatrix<unsigned> nn(d1,d2), ee(d1,d2);
   wbMatrix<int> gotr(d1,d2); 

   wbMatrix< wbvector<double> > aa; 
   wbvector<double> ax;

   wbMatrix<wbstring> ss; 

   if (nargin>1) {
      OPTS opts; opts.init(argin+1, nargin-1);

      if (opts.getOpt("-v")) { vflag|=2; } else
      if (opts.getOpt("-q")) { vflag =0; }

      if (opts.getOpt("~r"))
           { rflag=0; }
      else { opts.getOpt("-r"); } 

      if (opts.getOpt("~s")) { sflag=0; } 

      aflag=opts.getOpt("-a");
      opts.getOpt("nmax",nmax);
      opts.getOpt("eps1",eps1);
      opts.getOpt("eps2",eps2);

      opts.checkAnyLeft();
   }

   if (aflag) { aa.init(d1,d2); }
   if (sflag) { ss.init(d1,d2); }

   for (i=0; i<n; ++i) { a=Ai=A[i]; e2=-1;
       e1=Wb::Rational( 
          a,p,q,&(err[i]),&(nn[i]),
          aflag ? &(aa[i]) : NULL, nmax,-1,eps1,eps2,vflag>1
       );

       if (rflag && (e1 || abs(p)>999 || abs(q)>999)) {
          if ( (x2=a2=Ai*Ai) > eps1) {
          e2=Wb::Rational(
             a2,p2,q2, &err2, &n2, 
             aflag ? &ax : NULL, nmax,-1,eps1,eps2,vflag>1 
          );
       }}

       if (!e1 && !e2) {
          double q1=fabs(a*a-x2), q2=fabs(a2-x2);
          if (q1>q2) {
             if (vflag>1) PRINTF("  ==> picking "
                "sqrt() over plain (e=%.3g / %.3g)\n",sqrt(q2),sqrt(q1));
             e1=-2;
          }
       }

       if (!e1) {
          if (err[i]>1E-12 && vflag) wblog(FL,
             "WRN fixing number %8.5g => %8.5g (%.3g)",A[i],a,err[i]);
          r2+=(err[i]*err[i]);

          if (q<0) { p=-p; q=-q; } 
          A[i]=a; P[i]=p; Q[i]=q;
          if (sflag) {
             if (q==1)
                  { sprintf(str,"%ld",p); }
             else { sprintf(str,"%ld/%ld",p,q); }
             ss[i]=str;
          }
          continue;
       }
       if (e2) { 
          if (sflag) { sprintf(str,"%.6g",A[i]); ss[i]=str; }
          ee[i]=1; continue;
       }

       err[i]=err2; nn[i]=n2; if (aflag) { ax.save2(aa[i]); }

       a=sqrt(a2); if (Ai<0) { a=-a; }
       x=Ai-a; r2+=(x*x);

       if (::fabs(x)>1E-12  && vflag) { wblog(FL, 
          "WRN fixing number %8.5g => %8.5g (%.3g / %.3g)",A[i],a,x,err[i]);
       }
       A[i]=a; P[i]=p2; Q[i]=q2; gotr[i]=(Ai<0 ? -1:1);

       if (!sflag) { continue; }
       else { str[0]=0; }

       if (q2==1) { sprintf(str,"%ssqrt(%ld)",  a<0 ? "-":"",p2); } else
       if (p2==1) { sprintf(str,"%s1/sqrt(%ld)",a<0 ? "-":"",q2); }

       if (!str[0]) { x=sqrt(double(p2));
          if (fabs(x-round(x))<1E-14) {
          sprintf(str,"%s%ld/sqrt(%ld)",a<0 ? "-":"",long(x),q2);
       }}
       if (!str[0]) { x=sqrt(double(q2));
          if (fabs(x-round(x))<1E-14) {
          sprintf(str,"%ssqrt(%ld)/%ld",a<0 ? "-":"",p2,long(x));
       }}
       if (!str[0]) {
          sprintf(str,"%ssqrt(%ld/%ld)",a<0 ? "-":"",p2,q2);
       }
       ss[i]=str;
   }

   nx=ee.sum();

   if (vflag && nargout<2) {
      if (nx) wblog(FL,
         "WRN FixRational() %d/%d values unaltered (%d)",nx,n,nn.max());
      else if (sqrt(r2)>eps2) wblog(FL,
         "TST FixRational() skipped weight is %.3g (%d/%d; %d)",
          sqrt(r2),nx,n,nn.max()
      );
   }

   if (sflag)
        { argout[0]=ss.toMx(); }
   else { argout[0]= A.toMx(); } 

   if (nargout>1) {
      MXPut IO; IO.add(P,"P").add(Q,"Q");
      IO.add(gotr,"rflag").add(ee,"missed").add(nx,"nmiss");
      if (aflag) { IO.add(aa,"a"); }
      IO.add(err,"relerr").add(sqrt(r2),"rtot").add(nn,"niter")
        .add(nmax,"nmax").add(eps1,"eps1").add(eps2,"eps2")
        .save2(argout[1]);
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbrat"); }
};

