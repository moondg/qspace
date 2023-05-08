// #define __WBDEBUG__

char USAGE[] =
/* ========================================================================
 * bicg2SiteQS.cc */

"   USAGE: [Psi1, Psi2, [,info]] =                                      \n\
   bicg2SiteQS(Psi1, Psi2, xPsi1, xPsi2, ic1, ic2, lrdir, z, ...        \n\
     { A, ia, la [,opA] }, ...                       % H_{L/R, local}   \n\
     {{ A, ia, la [,opA] }, { B, ib, lb [,opB] }}... % H_coupling       \n\
     OPTS ...                                                           \n\
   );                                                                   \n\
                                                                        \n\
     C++ program that solves (H-z)|A>=|B> using Abelian                 \n\
     symmetries followed by an orthonormaliation step either            \n\
     for the `left' (Psi1, LR) or the `right' (Psi2, RL) QSpaceRM.        \n\
                                                                        \n\
   VARIABLES                                                            \n\
                                                                        \n\
     Psi[12]  nearest neighbor MPS blocks    (rank-3, real)             \n\
     xPsi[12]                                                           \n\
                                                                        \n\
         partial overalp <chi|Psi> of Psi with the                      \n\
         correction vector chi to the left(right) of Psi1(2).           \n\
         Note that this must have the same structure as Psi[12].        \n\
                                                                        \n\
     ic[12]   index that connects (x)Psi1 to (x)Psi2                    \n\
                                                                        \n\
     lrdir = 'LR' or 'RL'                                               \n\
                                                                        \n\
        direction (left to right, right to left) that specifies order   \n\
        of SVD on output.                                               \n\
                                                                        \n\
     z  = omega + i * eta                     (scalar, complex)         \n\
                                                                        \n\
     A, opA, ia, la=1,2                                                 \n\
                                                                        \n\
        QSpaceRM operator that acts on index ia of Psi[idx]               \n\
        A is a scalar or vector QSpaceRM        (rank-2, real)            \n\
                                                                        \n\
   Both Psi1 and Psi2 are assumed to have at least one (local) index    \n\
   with dimension 1 in its QSpaceRM; data will be rearranged accordingly. \n\
                                                                        \n\
   By default all operators are assumed to be contracted                \n\
   with their second index onto PSI, i.e. as for H in H|psi>;           \n\
   opA() can be used to change this in analogy to BLAS routine dgemm()  \n\
                                                                        \n\
   OPTIONS for biconjugate gradient                                     \n\
                                                                        \n\
    'maxit',.. maximum number of iteration      (32)                    \n\
    'btol',... tolerance btol*norm(b)           (1E-16)                 \n\
    'rtol',... tolerance rtol*initial residual  (1E-4)                  \n\
    'disp',... more detailed log to stdout (0:none ... 4:full; (1))     \n\
                                                                        \n\
     stop criteria: absolute residual < MAX(btol*norm(b), rtol*res0)    \n\
     relres = min abs. residual / norm(b)                               \n\
     resvec = sequence of abs. residual                                 \n\
                                                                        \n\
     Stop criteria on residual: MAX(btol*|b|, MIN(|b|, rtol*|r0|))      \n\
                                                                        \n\
   OPTIONS for subsequent orthonormalization                            \n\
                                                                        \n\
    'Nkeep',.. maximum dimension connecting Psi1 to Psi2 (64)           \n\
    'stol',... rel. tolerance on SVD (1E-4)                             \n\
    'twoSite'  2-site orthonormalization (only relevant if 1-site bicg) \n\
                                                                        \n\
   OUTPUT                                                               \n\
                                                                        \n\
     Psi1(2)   updated and properly orthonormalized Psi1(2) using lrflag\n\
     info      info structure containing                                \n\
      .resvec  vector of residuals                                      \n\
      .svd     vector of singular values                                \n\
                                                                        \n\
   AW (C) May 2006                                                       \n";

#define WB_CLOCK

#include "wblib.h"

#include "mpsortho.cc"
#include "cpat.hh" 

#define _TQ int
#define _TD wbcomplex

template<class TQ, class TD>
int checkStag(
   const QSpaceRM<TQ,TD> &A,
   const QSpaceRM<TQ,TD> &B,
   const double &aa
);

template<class TQ, class TD>
mxArray* RunBiCG(
   QSpaceRM<TQ,TD> &X,       
   const QSpaceRM<TQ,TD> &BZ,
   CPAT<TQ,TD> &CP,        
   IDAT &idat,             
   unsigned maxit,         
   double btol,            
   double rtol,            
   unsigned disp
){
   unsigned iter, stag, ibest=0;
   TD gi, alpha, beta, rho, rho1, ptq;
   double normr, minr, nb, tolr;

   QSpaceRM<TQ,TD> AX, XX, XBEST, Q, QT, R, RT, P, PT;
   char istr[64];

   static unsigned ncall=0;

   if (!maxit) wblog(FL,"ERR maxit=%d (%d)???", maxit);

   if (!X.isEmpty())
        { if (!X.hasSameQ(BZ)) {
            X.put("XX"); BZ.put("BB");
            wblog(FL,"WRN QIDX mismatch !?");
        }}
   else {
   X=BZ; X.setDATA(0); }

   XBEST=X;

   nb=sqrt(ABS(BZ.norm2()));
   if (nb==0) {
       X=BZ; idat.init(1); idat.niter=0;
       return mxCreateStructMatrix(1,1,0,NULL);
   }
   idat.init(maxit); idat.flag=1;

   BZ.minus( CP.HZTimes(X,AX), R);
   CP.nr_check_ctr=3; 

   normr=sqrt(ABS(R.norm2()));
   idat.resvec[0]=normr;

   tolr=MAX(btol*nb, MIN(nb, rtol*normr));

   if (disp<=1) { if (!(ncall++)) {
      wblog(FL,
      "%N%s tol setting\n"
      " ¤  |b|=%8.3g (btol=%g)\n"
      " ¤  |r|=%8.3g (rtol=%g)\n"
      "--> tol= %.3g = %.3g·|r|%N", __FUNCTION__,
      nb, btol, normr, rtol, tolr, tolr/normr);
   }}
   else {
      wblog(FL,"%s norm(b)=%.3g, norm(r)=%.3g --> res.tol: %.3g * res0",
      " · ", nb, normr, tolr/normr);
   }

   rho1=rho=1.; stag=0; minr=normr; istr[0]=0;

   if (normr<=tolr) {
      iter=ibest=0;
      idat.flag=0;
   }
   else {
      RT=R; P=R; PT=RT;
      for (iter=1; iter<maxit; iter++) {

         if (disp>3) {
            printf("%s  %.4s:%d  %6d/%d  ", iter>1 ? "\r" : "\n", FL,
            iter+1, maxit); fflush(0);
         }

         rho1=rho; rho=R.scalarProd(RT);
         if (rho==0. || rho.isinf()) {
            idat.flag=4; sprintf(istr,"rho=%s", rho.toStr().data);
            break;
         }

         if (iter>1) {
            beta=rho/rho1;
            if (beta==0. || beta.isinf()) {
               idat.flag=4; sprintf(istr,"beta=%s", beta.toStr().data);
               break;
            }

            P *=beta;        P +=R ;
            PT*=beta.conj(); PT+=RT;
         }

         CP.HZTimes(P, Q);
         CP.HZTimes(PT,QT,'C');

         ptq=Q.scalarProd(PT);
         if (ptq==0.) {
            idat.flag=4; sprintf(istr,"ptq=%s", ptq.toStr().data);
            break;
         }

         alpha=rho/ptq;
         if (alpha.isinf()) {
            idat.flag=4; sprintf(istr,"alpha=%s", alpha.toStr().data);
            break;
         }

         if (alpha==0.) stag=1;
         if (stag==0) stag=checkStag(P, X, alpha.abs());

         X.Plus(P,alpha);

         BZ.minus(CP.HZTimes(X,AX), XX);
         normr=sqrt(ABS(XX.norm2())); idat.resvec[iter]=normr;

         if (normr<=tolr) {
            idat.flag=0; minr=normr; ibest=iter;
            break;
         }
         if (stag) {
            idat.flag=3; sprintf(istr,
            "stag=%d (normr=%.3g, tolr=%.3g)", stag, normr, tolr);
            break;
         }

         if (normr<minr) { XBEST=X; minr=normr; ibest=iter; }

         R .Minus(Q,  alpha);
         RT.Minus(QT, alpha.conj());
      }
   }

   if (disp>1) printf("\n\n");

   if (iter<maxit)
   idat.resvec.Resize(iter+1);

   if (idat.flag) {
   X=XBEST;
   }

   idat.niter=ibest+1;
   idat.relres=minr/nb;

   gi = X.scalarProd(BZ);

   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL, S, "niter",  numtoMx(idat.niter));
   mxAddField2Scalar(FL, S, "maxit",  numtoMx(maxit));
   mxAddField2Scalar(FL, S, "flag",   numtoMx(double(idat.flag)));
   mxAddField2Scalar(FL, S, "bx",     numtoMx(gi)); 

   mxAddField2Scalar(FL, S, "norma2", numtoMx(X.scalarProd(X)));
   mxAddField2Scalar(FL, S, "normb",  numtoMx(nb)); if (disp) {
   mxAddField2Scalar(FL, S, "normr",  numtoMx(normr)); }
   mxAddField2Scalar(FL, S, "relres", numtoMx(idat.relres));
   mxAddField2Scalar(FL, S, "facres", numtoMx(minr/idat.resvec[0]));
   mxAddField2Scalar(FL, S, "resvec", idat.resvec.toMx());

   { MXPut T(FL,"tol");
     T.add(tolr,"tolres").add(rtol,"rtol").add(btol,"btol")
      .add(nb,"nb").add(normr,"nr");
     T.add2Struct(S);
   }

   if (disp>1) {
   mxAddField2Scalar(FL, S, "alpha",  numtoMx(alpha));
   mxAddField2Scalar(FL, S, "beta",   numtoMx(beta));
   mxAddField2Scalar(FL, S, "rho",    numtoMx(rho));
   mxAddField2Scalar(FL, S, "rho1",   numtoMx(rho1));
   mxAddField2Scalar(FL, S, "ptq",    numtoMx(ptq));
   }

   if (idat.flag) {
      sprintf(str, "ERR @ %d/%d iterations (%s)", iter, maxit, istr);
   }
   else str[0]=0;

   mxAddField2Scalar(FL, S, "msg", mxCreateString(str));

#if 0
   CP.HZTimes(X,XX); alpha=XX.scalarProd(X);
   wblog(FL,"TST aHa = %.8g %+.8gi", alpha.r, alpha.i);

   sprintf(str, "%s %3d/%d iterations @ res=%10.4g (%d), |G|=%8.3g",
         idat.flag ? "ERR" : "TST",
         iter, maxit, idat.relres, idat.flag, gi.abs());
   printf("%s:%d %s\n", FL, str);
   wblog(FL,"    normr = %.4g (%.4g; nb=%.4g)", normr, minr, nb);

   X.put("Xout","caller"); BZ.put("BZ","caller");
   CP.HZTimes(X,XX);     XX.put("XX","caller");
   CP.HZTimes(X,XX,'C'); XX.put("XC","caller");

   BZ.minus( CP.HZTimes(X,AX), XX); XX.put("R","caller");

#endif

   return S;
}

template<class T>
double nrmCheckStag(const wbArray<T> &A, const wbArray<T> &B) {
   unsigned i, s=A.SIZE.prod();
   T x; const T *a=A.data, *b=B.data;
   double nrm=0;

   if (!A.sameSize(B)) wblog(FL,
   "ERR cannot compare arrays of different size");

   for (i=0; i<s; i++) {
      if (b[i]!=0) { x=ABS(a[i]/b[i]); nrm=MAX(nrm,ABS(x)); }
      else if (a[i]!=0) return Inf;
   }
   return nrm;
}

template<class TQ, class TD>
int checkStag(
   const QSpaceRM<TQ,TD> &A,
   const QSpaceRM<TQ,TD> &B,
   const double &aa
){

   unsigned i,i1,i2,d,k,l, m=A.QIDX.dim1, n=A.QIDX.dim2;
   wbvector<unsigned> D;
   wbperm P;
   wbMatrix<TQ> QQ;
   double nrm=0;

   if (A.isEmpty() || B.isEmpty()) {
      wblog(FL,"WRN Empty QSpaces P or X (%d,%d) ???",
      A.isEmpty(), B.isEmpty()); return 1;
   }

   if (n!=B.QIDX.dim2 || A.rank()!=B.rank()) wblog(FL,
   "ERR Incompatible QSpaces (%d,%d;%d,%d).",
   n, B.QIDX.dim2, A.rank(), B.rank());

   QQ.Cat(1, A.QIDX, B.QIDX);
   QQ.groupRecs(P,D);

   for (k=l=i=0; i<QQ.dim1; i++, l+=d) {
       d=D[i];

       if (d==2) {
          i1=P[l]; i2=P[l+1]; k++;

          if (i1<m && i2>=m) { i2-=m; } else
          if (i2<m && i1>=m) { SWAP(i1,i2); i2-=m; } else
          wblog(FL,"ERR QSpaces have non-unique QIDX ???");

          nrm = MAX(nrm, nrmCheckStag(*A.DATA[i1], *B.DATA[i2]));
       }
       else if (d>3) wblog(FL,
       "ERR QSpaces have non-unique QIDX ???");
   }

   if (k==0) wblog(FL,
   "WRN checkStag - no overlap in QSpaceRM !?");

   return (nrm*aa < DBL_EPSILON);
}

char RFLAG='r'; 

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

#ifdef WB_CLOCK
   Wb::Clock clk("bicg:2Site:all",0);
   ARG_CHECK=1;
#else
   ARG_CHECK=0;
#endif

   int i,r=-1; char ldir;
   const char cflag=(typeid(_TD)==typeid(wbcomplex) ? 'C' : 0);
   wbcomplex z;
   double btol=1E-8, rtol=1E-4, stol=1E-4;
   unsigned nc, ic1, ic2, r1, r2, m, k, K,
      disp=1, Nkeep=64, maxit=32,
      twoSite=0,  
      on2Site;    

   mxArray *S=NULL, *a;

   QSpaceRM<_TQ,_TD> xPsi1, xPsi2, PSIB;
   QSpaceRM<_TQ,_TD> Psi1, Psi2, PSI, A1, A2;
   CPAT<_TQ,_TD> CP;
   IDAT idat;
   OPTS opts;

   if (nargin) if (checkHelpVersion(argin[0])) { return; }

   if (nargin) {
      const mxArray *ac=argin[nargin-1];
      if (mxIsChar(ac) && mxGetNumberOfElements(ac)<32) {
         str[0]=0; mxGetString(ac,str,32);
         if (!strcmp(str,"NO_RFLAG")) {
            RFLAG=0; nargin--;
            wblog(FL,"TST *** UNSET RFLAG ***");
         }
      }
   }

   if (nargin <9) wbdie(FL,"Too few input arguments."  ); else
   if (nargout<2) wbdie(FL,"Too few output arguments." ); else
   if (nargout>3) wbdie(FL,"Too many output arguments.");

   for (k=i=0; i<4; i++) {
      try { mxIsQSpace(FL,argin[i],r,cflag); }
      catch (...) {
         wblog(FL,"ERR argument #%d must be a valid QSpaceRM",i+1);
      }
      if (mxIsEmptyQSpace(argin[i])) {
         if (i<2 || k>1) {
            if (i<2)
            wblog(FL,"ERR argument #%d is empty QSpaceRM !?",i+1); else
            wblog(FL,"ERR argument #%d and #%d are empty !?", i, i+1);
         }
         else if (i>=2) k++;
      }
   }

   Psi1 .init(FL, argin[0],RFLAG); r1=Psi1.rank();
   Psi2 .init(FL, argin[1],RFLAG); r2=Psi2.rank();
   xPsi1.init(FL, argin[2],RFLAG);
   xPsi2.init(FL, argin[3],RFLAG);

   for (i=4; i<6; i++) if (!Mx::IsNumber(argin[i])) wblog(FL,
   "ERR argument #%d must be a valid index.", i+1);

   if (mxGetNumber(argin[4], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[5], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
       ic1--; ic2--;  
   }
   else wblog(FL,
   "ERR arguments #5 or #6 - index out of bounds (%d,%d)", ic1, ic2);

   if (RFLAG) {
      ic1=r1-ic1-1; ic2=r1-ic2-1;
   }

   ldir=mxGetLRDir(FL,argin[6],6);

   if (mxGetNumber(argin[7], z)) wblog(FL,
   "ERR z must be complex scalar (arg #7)!%N%N%s%N",str);

   if (isnan(z.r) || isnan(z.i)) wblog(FL,
   "ERR invalid z=%s (contains NaN; arg #7) !?",z.toStr().data);

   i=8; m=isCPAT(argin[i]);
   if (m) {
      const mxArray *a = (m>1 ? mxGetCell(argin[i],0) : argin[i]);
      twoSite=isCPattern(a)-1;

      for (i++; i<nargin; i++)
      if (!isCPAT(argin[i])) break;
   }

   m=i; nc=m-8;

   if (!nc) { 
   wblog(FL,"ERR CPAT required as argument #9."); }

   for (; i<nargin; i++) if (mxIsCell(argin[i]))
   wbdie(FL, 
   "there should be no further cells after contraction patterns!");

   opts.init(argin+m,nargin-m);

   opts.getOpt(0,0, "disp",  disp );
   wbSetLogLevel(disp>1 ? 1 : 0); 

   wblog(FL,"<%d> using %sSite method",twoSite+1, twoSite ? "two":"one");

   opts.getOpt(FL, "maxit", maxit);
   opts.getOpt(FL, "rtol",  rtol );
   opts.getOpt(FL, "btol",  btol );
   opts.getOpt(FL, "stol",  stol );
   opts.getOpt(FL, "Nkeep", Nkeep);

   on2Site=opts.getOpt(FL,"twoSite");

   opts.checkAnyLeft(FL);

   if (rtol<0) wblog(FL,"ERR Invalid rtol=%g", rtol);
   if (btol<0) wblog(FL,"ERR Invalid btol=%g", btol);
   if (stol<0) wblog(FL,"ERR Invalid stol=%g", stol);
   if (rtol==0 && btol==0)
   wblog(FL,"ERR rtol or btol must be set finite.");

   wbSetLogLevel(disp ? 1 : 0);

   CP.init(
     z, argin+8, nc,
     ic1, r1, ic2, r2, twoSite, ldir, RFLAG
   );

   twoSiteInit(xPsi1, xPsi2, PSIB, k, ic1, ic2, twoSite, ldir);
   twoSiteInit( Psi1,  Psi2, PSI,  K, ic1, ic2, twoSite, ldir);

   if (k!=K) wblog(FL,"ERR K=%d,%d ???", k, K); 

   if (isFlagGlobal("DEBUG")) {
      wblog(FL,"TST DEBUG is on.");
      CP.put("CP");
   }

   try {
     #ifdef WB_CLOCK
      Wb::Clock cl2("bicg:2Site:actual",0);
     #endif

     S=RunBiCG(
        PSI, PSIB, CP,
        idat, maxit, btol, rtol, disp
     );
   }

#ifdef MATLAB_MEX_FILE

   catch (char* istr) {
       CP.put("CP","caller"); PSI.put("PSI","caller");
       if (istr && !strcmp(istr,"\n"))
            wblog(FL, "ERR gotcha (RunBiCG failed)");
       else wblog(FL, "ERR gotcha (RunBiCG failed, %s)", istr);
   }
   catch (...) {
       CP.put("CP","caller"); PSI.put("PSI","caller");
       wblog(FL, "ERR gotcha (RunBiCG failed ...)");
   }

#else

   catch (char* istr) { wblog(FL, "ERR RunBiCG failed (%s).", istr); }
   catch (...)        { wblog(FL, "ERR RunBiCG failed."); }

#endif

   if (!twoSite && on2Site) {
      QSpaceRM<_TQ,_TD> PSIX;
      twoSite=on2Site;

      if (ldir==+1) {
         PSI.permuteLastTo(ic1,PSIX);
         twoSiteInit(PSIX, Psi2, PSI,  K, ic1,ic2,twoSite, ldir);
      }
      else {
         PSI.permuteLastTo(ic2,PSIX);
         twoSiteInit(Psi1, PSIX, PSI,  K, ic1,ic2,twoSite, ldir);
      }
   }

   a=mpsOrthoQS(
      PSI, A1, A2,
      K, Nkeep, stol, disp ? (disp==1 ? "ds" : "dS") : ""
   );

   if (a && S)
   mxAddField2Scalar(FL, S, "Iortho", a);

   twoSiteFinal(Psi1, Psi2, A1, A2, ic1, ic2, twoSite, ldir);

   argout[0]=A1.toMx(RFLAG);
   argout[1]=A2.toMx(RFLAG);

   if (nargout>2) argout[2]=S; else
   if (S) mxDestroyArray(S);
};

