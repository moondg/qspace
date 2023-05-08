/* ==================================================================
 * curvature.cc
 *
 *    C++ test routine.
 *
 * This is a MEX-file for MATLAB.
 * ================================================================== */

char USAGE[]="\n   c=curvate(xd,yd,fac)\n\n";

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   unsigned i,j,i1,i2,m,n;

   wbvector<double> xd,yd,cc,x1,y1,x2,y2;
   wbperm P;
   double fac,d1,d2,dbl,eps;

   if (nargin!=3) { usage(); return; }

   xd.init(FL,argin[0],"zflags");
   yd.init(FL,argin[1],"zflags");

   if (mxGetNumber(argin[2],fac)) wbdie(FL,str);

   n=xd.len;

   for (i=1; i<n; i++) if (xd[i]!=0) { eps=fabs(xd[i]); break; }
   for (i++; i<n; i++) if (xd[i]!=0) { eps=MIN(eps,fabs(xd[i])); }
   eps*=5;

   if (nargin>3) {
      OPTS opts;
      opts.init(argin+3, nargin-3);

      opts.getOpt("eps", eps);

      opts.checkAnyLeft(); 
   }

   if (yd.len!=n || n<4) wblog(FL,
   "ERR size mismatch of x-y data (%d/%d)",xd.len,yd.len);
   if (fac<=1) wblog(FL,"ERR invalid fac=%g",fac);

   xd.Sort(P);
   yd.Permute(P);

   cc.init(n);
   x1.init(n-1); x2.init(n-2);
   y1.init(n-1); y2.init(n-2);

   for (m=n-1, i=0; i<m; i++) {
      if (xd[i+1]==xd[i])
      wblog(FL,"ERR degenerate x-data! (%d/%d)",i+1,n);

      y1[i]=(yd[i+1]-yd[i])/(xd[i+1]-xd[i]);
      x1[i]=0.5*(xd[i+1]+xd[i]);
   }

   for (m=n-2, i=0; i<m; i++) {
      y2[i]=(y1[i+1]-y1[i])/(x1[i+1]-x1[i]);
      x2[i]=0.5*(x1[i+1]+x1[i]);
   }

   wblog(FL,"fac=%g; eps=%g",fac,eps);

   for (i1=i2=i=0; i<n; i++) {
       d1=xd[i]/fac;
       d2=xd[i]*fac; if (xd[i]<0) SWAP(d1,d2);

       if (d1>xd[i]-eps) d1=xd[i]-eps;
       if (d2<xd[i]+eps) d2=xd[i]+eps;

       for (j=i1; j<i; j++) if (xd[j]>d1) break; i1=j;
       for (j=i2; j<n; j++) if (xd[j]>d2) break; i2=(j>0 ? j-1 : 0);

       for (;;) {
          m=i2-i1+1; if (m>2) break;
          if (i1) i1--; if (i2+1<n) i2++;
       }

       for (dbl=0.,j=i1+2; j<i2; j++) dbl+=y2[j-2];
       cc[i]=dbl/(double)(m-2);
   }

   if (nargout==1) argout[0]=cc.toMx(); else
   if (nargout==2) {
      argout[0]=xd.toMx();
      argout[1]=cc.toMx();
   }
};

