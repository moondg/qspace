char USAGE[]=
/* =====================================================================
 * wbeigs.cc */

"   psi=wbeigs(H, psi [, task, m])                                       \n\
                                                                         \n\
      H     symmetric matrix                                             \n\
      task  task for arpack ('SA'=smallest algebraic, i.e. ground state) \n\
      m     number of states to keep (1)                                 \n\
                                                                         \n\
   Wb,Nov11,06                                                           \n\
";

/* This is a MEX-file for MATLAB.
 * Reference: matlab-14.3/toolbox/matlab/sparfun/eigs.m
 * ===================================================================== */

#include "wblib.h"
#include "wbeigs.hh"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   int e, ido=0, info=0;
   unsigned i1, i2, N, ncv, nev=1, mode=1, maxit=300;
   wbMatrix<double> H;
   wbvector<double> psi,Hpsi;
   double sigma=0.;

   wbMatrix<double> X,V,W,EV;
   wbvector<double> E,x0;

   wbvector<int>  ipntr(15); 

   wbvector<int> iparam(11), select;

   char task[]="SA"; 
   double tol=-1;

   MX_CHECK_HELPER_NARGS(2,-1,3); 

   H  .init(argin[0]); N=H.dim2;
   psi.init(FL,argin[1]);

   if (psi.isEmpty()) psi.init(N);
   else {
      if (psi.len!=N) wblog(FL,
      "ERR Input vector does not match H (%d/%d)", psi.len, N);
      info=1;
   }

   nev=MIN(6U,N);

   if (nargin>2) {
      str[0]=0; e=mxGetString(argin[2],str,STRLEN);
      if (e || strlen(str)!=2) wblog(FL,
         "ERR Failed to read task `%s' (arg 3)",str);
      strcpy(task,str);
   }

   if (nargin>3) {
      OPTS opts;
      opts.init(argin+3, nargin-3);

      opts.getOpt("tol",   tol   );
      opts.getOpt("nev",   nev   ); ncv=MIN(N,3*nev);
      opts.getOpt("ncv",   ncv   );
      opts.getOpt("maxit", maxit );

      opts.checkAnyLeft(); 
   }
   else ncv=MIN(N,3*nev);

   if (ncv<6U) ncv=MIN(6U,N);

   if (ncv>N) { wblog(FL,
      "WRN NCV too large ­ reset %d -> %d", ncv, N);
       ncv=N;
   }

   E.init(nev);       
   x0=psi;            
   select.init(ncv);

   V.init(ncv,N);
   W.init(3,N);       
   X.init(ncv+8,ncv); 

   iparam[0]=1;       
   iparam[2]=maxit;
   iparam[6]=mode;    

   unsigned iter=0;

   while (ido!=99) {

      AUPDFUN( 
         ido, 'I', N, task, nev, tol, x0.data,
         ncv, V.data, N, iparam.data, ipntr.data, W.data, X.data, X.dim1*X.dim2,
         info
      );

      if (info) {
         printf("\nWRN %3d: N=%d, task='%s', nev=%d, ncv=%d, m*n=%d, tol=%g\n",
         iter, N,task,nev,ncv, X.dim1*X.dim2, tol);
      }

      if (info<0) wblog(FL,
         "ERR ARPACK routine %s returns info=%d (ido=%d)\nmsg: %s",
          AUPDSTR, info, ido, info2str(info));
      if (info>0) wblog(FL,
      "WRN ARPACK routine %s returns info=%d (ido=%d)",AUPDSTR,info,ido);

      i1=(ipntr[0]-1)/N; 
      i2=(ipntr[1]-1)/N;

      if (ipntr[0]%N!=1 || ipntr[1]%N!=1) wblog(FL, 
      "ERR ipntr(1:2) does not refer to start of column !?\n"
      "[%d %d] for workd = (%dx%d)",ipntr[0]-1,ipntr[1]-1,W.dim1,W.dim2);
      if (i1>=W.dim1 || i2>=W.dim1 || i1==i2) wblog(FL,
      "ERR invalid index into WORKD (%d,%d; %d) !?",i1,i2,W.dim1);

      if (ido==1 || ido==-1) {
         wbvector<double> x,Hx;

         x .init2ref(N, W.rec(i1));
         Hx.init2ref(N, W.rec(i2));

         Wb::MatProd(H,x,Hx);
      }
      else if (ido!=99)
      wblog(FL,"ERR %s returned ido=%d !?", AUPDSTR, ido);

   }

   EV.init(nev,N);
   EUPDFUN( 
      (nargout>=2), 'A', select.data,
      E.data, EV.data, N, sigma,
      'I', N, task, nev, tol, x0.data, ncv,
      V.data, N, iparam.data, ipntr.data, W.data, X.data, X.dim1*X.dim2, info
   );

   if (info<0) wblog(FL,
      "ERR ARPACK routine %s returns info=%d\n(%s)",
       EUPDSTR, info, info2str(info));
   if (info<0) wblog(FL,
   "ERR ARPACK routine %s returns %d !?", EUPDSTR, info);
   if (info>0) wblog(FL,
   "<i> ARPACK routine %s returns info=%d",EUPDSTR, info);

   if (nargout  ) argout[0]=E.toMx();
   if (nargout>1) argout[1]=EV.toMx('r');

#ifndef MATLAB_MEX_FILE
   if (ido==99 && info==0) {
      wbvector<double> x,Hx,ee(nev);

      printf("  %d/%d eigenvalues and -vectors:\n\n",EV.dim1,EV.dim2);
      for (unsigned i=0; i<nev; i++) {
         x .init2ref(N, EV.rec(i)); ee[i]=fabs(x.norm()-1);
         Hx.init2ref(N, W.rec(1));
         Wb::MatProd(H,x,Hx);
         x*=E[i]; ee[i]=MAX((x-Hx).norm(),ee[i]);
         printf("%6d  %8.4g  %12.4g\n",i,E[i],ee[i]);
      }

      if (ee.norm()<1E-12)
           printf("\n  SUC :)\n");
      else printf("\n  ERR ???\n");
   }
   else printf("\n  ERR ???\n");
#endif

   if (nargout>2) {
      mxArray *S=mxCreateStructMatrix(1,1,0,NULL);
      x0.init(ncv, X.data+ipntr[8]);

      sprintf_str("EIGS w/mode `%s' %sconverged %d/%d (%d)",
      task,iparam[4]!=(int)nev ? "NOT ":"",iparam[4],nev,N);

      mxAddField2Scalar(FL,S, "info",  mxCreateString(str));

      mxAddField2Scalar(FL,S, "niter", numtoMx(ipntr[14])); 

      mxAddField2Scalar(FL,S, "maxit", numtoMx(maxit));
      mxAddField2Scalar(FL,S, "ncv",   numtoMx(ncv) );
      mxAddField2Scalar(FL,S, "tol",   numtoMx(tol) );
      mxAddField2Scalar(FL,S, "iparam",iparam.toMx());
      mxAddField2Scalar(FL,S, "ipntr", ipntr.toMx() );
      mxAddField2Scalar(FL,S, "ido",   numtoMx(ido) );
      mxAddField2Scalar(FL,S, "err",   numtoMx(info));

      mxAddField2Scalar(FL,S, "res",   x0.toMx());
      mxAddField2Scalar(FL,S, "W",     W.toMx());
      mxAddField2Scalar(FL,S, "X",     X.toMx());

      argout[2]=S;
   }

}

