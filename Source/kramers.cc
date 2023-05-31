
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to kramers.m // Wb,Feb14,19

#include "wblib.h"

void KKreal_linint(
   const wbvector<double> &X,
   const wbMatrix<double> &Yi,
   wbMatrix<double> &Yr,
   char gflag=0, double eps=1E-10,
   char vflag=0
);

void KKreal_disc(
   const wbvector<double> &X,
   const wbMatrix<double> &Yi,
   wbMatrix<double> &Yr,
   char gflag=0, double faceta=1.,
   char vflag=0
);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned i,N;
   wbvector<double> xi;
   wbMatrix<double> yi,yr; 
   char isdisc, vflag, gflag, r2i, tflag=0;
   double eps;
   OPTS opts;

   MX_CHECK_HELPER_NARGS(2,-1,2); 

   if (!Mx::IsDblArr(FL,argin[0]) || !Mx::IsDblArr(FL,argin[1]) ||
       mxIsComplex(argin[0]) || mxIsComplex(argin[1])) wblog(FL,
   "ERR first two arguments must be real arrays (%s)",
    mxGetClassName(argin[1]));

   if (!Mx::IsVector(argin[0])) wblog(FL,
   "ERR first argument must be vector");

   xi.init(FL,argin[0]); N=xi.len;
   yi.init(argin[1]);

   if (yi.dim1!=N && yi.dim2!=N) wblog(FL,
   "ERR Dimension mismatch of arguments (%dx%d; %d)",yi.dim1,yi.dim2,N);

   for (i=1; i<N; i++) if (xi[i]<=xi[i-1]) wblog(FL,
   "ERR first argument must contain data \nin strictly ascending order");

   if (yi.dim2!=N) { tflag=1; yi.Transpose(); }

   opts.init(argin+2, nargin-2);

   r2i    = opts.getOpt("r2i");
   vflag  = opts.getOpt("-v");

   gflag = (r2i ? 1 : 2);

   if (opts.getOpt("g1flag")) gflag=1; else
   if (opts.getOpt("g2flag")) gflag=2; else
   opts.getOpt("gflag",gflag);

   isdisc = opts.getOpt("-disc");
   if (!isdisc) {
      eps=1E-8*xi.aMin(1); 
      opts.getOpt("eps",   eps);
   }
   else {
      eps=1;
      opts.getOpt("faceta",eps);
   }

   opts.checkAnyLeft(); 

   if (isdisc)
        KKreal_disc  (xi,yi,yr, gflag, eps, vflag);
   else KKreal_linint(xi,yi,yr, gflag, eps, vflag);

   if (tflag) yr.Transpose();
   yr *= (r2i ? -1 : 1)/M_PI; 

   argout[0]=yr.toMx(); 

   if (nargout>1) {
      mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

      mxAddField2Scalar(FL,S,"isdisc", numtoMx(isdisc));
      mxAddField2Scalar(FL,S,"gflag",  numtoMx(gflag));
      mxAddField2Scalar(FL,S,"eps",    numtoMx(eps));
      mxAddField2Scalar(FL,S,"vflag",  numtoMx(vflag));

      argout[1]=S;
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in kramers"); }
   aclu.Check();
};

void KKreal_disc(
   const wbvector<double> &X,
   const wbMatrix<double> &Yi, 
   wbMatrix<double> &Yr,       
   char gflag, double faceta,
   char vflag
){
   wbvector<double> X2(X.len+1), dX(X.len);

   unsigned i,j,k, N=X.len, m=Yi.dim1;
   const double *x=X.data;
   double dij,eta,x0,xN, *x2=X2.data, *dx=dX.data;

   Wb::SigHandler SIG(FL); 

   if (vflag)
      wblog(FL,"%s (%dx%d; gflag=%d, faceta=%g)",
    __FUNCTION__,Yi.dim2, Yi.dim1, gflag,faceta);

   Yr.init(Yi.dim1, Yi.dim2);
   if (N<2) {
      wblog(FL,"WRN length of xdata to short.");
      return;
   }
   if (faceta<=0)
   wblog(FL,"ERR invalid faceta=%g",faceta);

   for (i=1; i<N; i++) x2[i]=0.5*(x[i-1]+x[i]);
   x2[0]=x[0]-0.5*(x[1]-x[0]);
   x2[N]=x[N-1]+0.5*(x[N-1]-x[N-2]);

   x0=x2[0]; xN=x2[N];

   for (i=0; i<N; i++) {
      dx[i]=x2[i+1]-x2[i]; if (dx[i]<=0) wblog(FL,
      "ERR x-data must be ascending! (dx[%d]=%g)",i,dx[i]);
   }

   if (faceta!=1.) dX*=faceta;

   for (i=0; i<N; i++) {
      SIG.check911();

      for (j=0; j<N; j++) {
         eta=dx[i]+dx[j];
         dij=log(  
           ( wbcomplex(x[j]-x2[i],eta) /
             wbcomplex(x[j]-x2[i+1],eta) ).abs2()
         );

         for (k=0; k<m; k++) Yr(k,i)+=dij*Yi(k,j);
      }

      if (gflag==1) {
         double d0,dN;
         d0= x0/x[i]*log( wbcomplex(x0-x[i],dx[i]).abs2() / (x0*x0) );
         dN=-xN/x[i]*log( wbcomplex(xN-x[i],dx[i]).abs2() / (xN*xN) );

         d0*=(dx[i]/dx[0]);
         dN*=(dx[i]/dx[N-1]);

         for (k=0; k<m; k++) Yr(k,i)+=(d0*Yi(k,0) + dN*Yi(k,N-1));
      }
      else
      if (gflag) wblog(FL,"ERR gflag=%d not implemented yet",gflag);
   }

   Yr*=0.5; 
}

void KKreal_linint(
   const wbvector<double> &X,
   const wbMatrix<double> &Yi, 
   wbMatrix<double> &Yr,       
   char gflag, double eps,
   char vflag
){
   unsigned i,j,k, N=X.len, m=Yi.dim1;
   const double *x=X.data, *yi;
   double *yr, dx1,dx2, xj, g1, g2, d;

   Wb::SigHandler SIG(FL); 

   if (vflag) wblog(FL,
      "%s (%dx%d; gflag=%d, eps=%.4g)",FCT,Yi.dim2,Yi.dim1,gflag,eps);

   Yr.init(Yi.dim1, Yi.dim2);

   for (k=0; k<m; ++k) {
      yi=Yi.rec(k); yr=Yr.rec(k);

      for (i=1; i<N; ++i) {
         SIG.check911();

         d=x[i]-x[i-1]; if (d<=0) wblog(FL,
           "ERR %s() strictly increasing x-data required",FCT);

         g1=yi[i-1]/d; g2=yi[i]/d;

         for (j=0; j<N; ++j) {
            dx1 = (j!=i-1 ? x[j]-x[i-1] : eps);
            dx2 = (j!=i   ? x[i]-x[j  ] : eps);

            yr[j] += (dx1*g2 + dx2*g1) * log(fabs((dx2/dx1)));
         }
      }

      d=yi[N-1]-yi[0];

      if (gflag==0) {
         for (j=0; j<N; ++j) yr[j]+=d;
      }
      else {
         double x0=x[0], xN=x[N-1], y0=yi[0], yN=yi[N-1];
         if (x0==0 || xN==0) wblog(FL,
            "ERR %s() got range [%.4g %.4g]",FCT,x0,xN);

         if (gflag==1) {
            double q0=y0*x0, qN=yN*xN;
            for (j=0; j<N; ++j) { xj=x[j];
               dx1=(j      ? x0-xj : eps);
               dx2=(j!=N-1 ? xN-xj : eps); if (xj==0) xj=eps;

               yr[j] += (
                   d + (q0/xj)*log(fabs(dx1/x0))
                     - (qN/xj)*log(fabs(dx2/xN))
               );
            }
         }
         else if (gflag==2) {
            double q0,qN;
            for (j=0; j<N; ++j) {
               xj=x[j]; if (xj==0) xj=eps;

               q0=x0/xj;
               if (fabs(q0)<1E6) { 
                  dx1=(j      ? x0-xj : eps);
                  q0 = y0*q0*(1+q0*log(fabs(dx1/x0))); }
               else {
                  q0 = y0*(-0.5 - 1./(3*q0) - 1./(4*q0*q0));
               }

               qN=xN/xj;
               if (fabs(qN)<1E6) { 
                  dx2=(j!=N-1 ? xN-xj : eps);
                  qN = yN*qN*(1+qN*log(fabs(dx2/xN))); }
               else {
                  qN = yN*(-0.5 - 1./(3*qN) - 1./(4*qN*qN));
               }

               yr[j] += (d + q0 - qN);
            }
         }
         else wblog(FL,"ERR");

         for (j=0; j<N; ++j) { if (std::isnan(yr[j])) {
            wblog(FL,"WRN %s() encountered NaN-value(s) in yr",FCT);
            break;
         }}
      }
   }
};

