// #define __WBDEBUG__

char USAGE[] =
/* ========================================================================
 * eigs2SiteQS.cc */

"   eigs2SiteQS                                                         \n\
                                                                        \n\
     C++ program that solves H|A>=E0|A> using Abelian symmetries        \n\
     followed by an orthonormaliation step either for the `left'        \n\
     (Psi1, LR) or the `right' (Psi2, RL) QSpace.                       \n\
                                                                        \n\
     The eigenvalue problem is solved iteratively                       \n\
     using ARPACK package (Sorensen et al.)                             \n\
                                                                        \n\
   ONE-SITE USAGE                                                       \n\
                                                                        \n\
     [Psi1, Psi2, [,info]] =                                            \n\
     eigs2SiteQS(Psi1, Psi2, ic1, ic2, lrdir, ...                       \n\
       { A, i1 [,opA] }, ...                           % H_{L/R, local} \n\
       {{ A, ia [,opA] }, { B, ib [,opB] }, ... }...   % H_coupling     \n\
       OPTS ...                                                         \n\
     );                                                                 \n\
                                                                        \n\
   TWO-SITE USAGE (NB! either i1 or i2 must be empty [])                \n\
                                                                        \n\
     [Psi1, Psi2, [,info]] =                                            \n\
     eigs2SiteQS(Psi1, Psi2, ic1, ic2, lrdir, ...                       \n\
       { A, i1, i2 [,opA] }, ...                       % H_{L/R, local} \n\
       {{A,i1,[] [,opA]}, {B,[],i2 [,opB]}, ... } ...  % H_coupling     \n\
       OPTS ...                                                         \n\
     );                                                                 \n\
                                                                        \n\
   VARIABLES                                                            \n\
                                                                        \n\
     Psi[12]  nearest neighbor MPS blocks    (rank-3, real)             \n\
     ic[12]   index that connects Psi1 to Psi2                          \n\
                                                                        \n\
     lrdir = 'LR' or 'RL' ('>>' or '<<')                                \n\
                                                                        \n\
        direction (left to right, right to left) that specifies order   \n\
        of SVD on output.                                               \n\
                                                                        \n\
     { A, i1 [,i2], opA }                                               \n\
                                                                        \n\
        QSpace operator contraction set that acts on index ia of Psi[la]\n\
        A is a scalar or vector QSpace        (rank-2, real)            \n\
                                                                        \n\
   By default all operators A are assumed to be contracted with         \n\
   their second index onto PSI=(Psi1,Psi2), i.e. as in A|psi>;          \n\
   opA() can be used to change this in analogy to BLAS routine dgemm()  \n\
                                                                        \n\
   OPTIONS for eigs() [ARPACK]                                          \n\
                                                                        \n\
    'xflag',.. which dimension of Hilbert space to expand to full       \n\
               randomizing newly acquired coefficient blocks looking up \n\
               input contraction pattern                                \n\
    'rnrm',..  overall amplitude of randomly initialized newly acquired \n\
               data blocks (1E-3)                                       \n\
                                                                        \n\
    'maxit',.. maximum number of iteration      (300)                   \n\
    'rtol', .. rel. tolerance (-1, ARPACK default (machine precision))  \n\
    'disp',... more detailed log to stdout (0:none ... 4:full; (1))     \n\
                                                                        \n\
   OPTIONS for subsequent orthonormalization                            \n\
                                                                        \n\
    'Nkeep',.. maximum dimension connecting Psi1 to Psi2 (64)           \n\
    'stol',... tolerance on SVD (i.e. on density matrix; 1E-4)          \n\
                                                                        \n\
   Test parameters                                                      \n\
                                                                        \n\
    'NCV'...   number of Lanczos vectors (arpack; max(20, 2*nev)).      \n\
    'nev'...   number of eigenvalues to compute (1).                    \n\
    'randall'  randomize input A matrizes.                              \n\
                                                                        \n\
   OUTPUT                                                               \n\
                                                                        \n\
     Psi1(2)   updated and properly orthonormalized Psi1(2) using lrflag\n\
     info      info structure containing                                \n\
      .info    info summary of eigs()                                   \n\
      .svd     vector of singular values                                \n\
                                                                        \n\
   AW (C) Dec 2006                                                     \n";

#define WB_CLOCK

#define QTYPE int
#define _TQ QTYPE
#define _TD double

#include "wblib.h"
#include "mpsortho.cc"
#include "cpat.hh" 

#include "wbeigs.hh"

unsigned NCV=20U;

template<class TQ, class TD>
mxArray* RunEigs(
   QSpace<TQ,TD> &PSI0,    
   CPAT<TQ,TD> &CP,        
   wbvector<unsigned> &xflag,
   char disp=0,
   double rnrm=1E-3,       
   unsigned maxit=300,     
   unsigned nev=1,         
   double tol=-1.,         
   char randall=0
){
   static unsigned check_Hconj_id=0;

   int ido=0,info=0,info2=0;
   char task[]="SA"; 
   double dbl,sigma=0.;

   unsigned i,j,m,r,i1,i2,ncv=0,N=0,count=0;
   wbvector<TD> A;
   wbvector<unsigned> sa;
   wbvector<double> E,res,ritz;
   wbMatrix<TD> V,W,X,EV;

   wbvector<int>  ipntr(15); 

   wbvector<int>  iparam(11), select;
   unsigned mode=1; 

   QSpace<TQ,TD> PSI, HPsi;

   Wb::SigHandler SIG(FL);

   PSI=PSI0;
   if (CP.CP.len) { 
      if (PSI.isEmpty()) {
         for (r=i=0; i<CP.CP.len; i++) { m=CP.CP[i].len;
         for (j=0; j<m; j++) r=MAX(r,CP.CP[i][j].ic); }

         PSI.init(0,r+1,CP.CP[0][0].A[0].QDIM);
         PSI.ExpandQ(CP, 1., xflag, disp);
      }
      else {
         PSI.ExpandQ(CP,rnrm,xflag,disp);

         if (rnrm!=0 && tol>=0) tol=-1;
      }
   }

   if (randall) PSI.setRand(1.);

#if 0
    wblog(FL,"TST %s",FCT);
    PSI0.put("PSI0"); PSI.put("PSIX");

    wbMatrix<double> dd(3,2); wbindex I(PSI0.rank(FL),'i');
    QSpace<TQ,TD> Q1,Q2,Q3;
    PSI0.contract(I,PSI0,I,Q1); dd(0,0)=Q1.DATA[0]->data[0]; dd(0,1)=sqrt(Q1.norm2());
    PSI0.contract(I,PSI, I,Q2); dd(1,0)=Q2.DATA[0]->data[0]; dd(1,1)=sqrt(Q2.norm2());
    PSI .contract(I,PSI, I,Q3); dd(2,0)=Q3.DATA[0]->data[0]; dd(2,1)=sqrt(Q3.norm2());

    dd.put("dd");

    for (i=1; i<6; i++) if (fabs(dd[i]-dd[0])>1E-8) {
        dd.print("dd"); Q1.put("Q1"); Q2.put("Q2"); Q3.put("Q3");
        wblog(FL,"ERR");
    }
#endif

   if (CP.CP.len) 
   if ((check_Hconj_id++)<12 || check_Hconj_id%37) {
      double h1,h2; QSpace<TQ,TD> RR(PSI); RR.setRand(1.);
      CP.HZTimes(RR,HPsi,'C'); h1=PSI.scalarProd(HPsi);
      CP.HZTimes(RR,HPsi    ); h2=PSI.scalarProd(HPsi);

      dbl=fabs(h1-h2)/MAX(fabs(h1),fabs(h2));
      if (dbl>1E-12) {  
         if (dbl<1E-8) wblog(FL,
            "WRN H only approx. Hermitian (%.3g) !?",dbl);
         else wblog(FL,
            "ERR H not Hermitian (%.3g, %.3g -> %.3g) !?",h1,h2,dbl);
      }
   }

   CP.HZTimes(PSI,HPsi);

   if (PSI.QIDX.dim1!=HPsi.QIDX.dim1) {
      PSI.put("PSI"); HPsi.put("HPSI");
      wblog(FL,"ERR block setting changes (%d->%d)",
      PSI.QIDX.dim1,HPsi.QIDX.dim1);
   }

   PSI.Sort(HPsi); 
   N=PSI.map2Vec(FL,A);

   if (maxit==0 || N<=1 || CP.CP.len==0) { str[0]=0; E.init(1);
      E[0]=PSI.scalarProd(HPsi);
   }
   else {

   if (nev>N) nev=N;

   ncv=MAX(2*nev,NCV); if (ncv>N) ncv=N;

   if (NCV<2*nev) wblog(FL,
   "WRN setting ncv=%d (nev=%d)",ncv,nev);

   E.init(nev);
   select.init(ncv);

   V.init(ncv,N);
   W.init(3,N);       
   X.init(ncv+8,ncv); 

   iparam[0]=1;     
   iparam[2]=maxit;
   iparam[6]=mode;  

   CP.nr_check_ctr=2;

   while (ido!=99) {
      SIG.check911();

      AUPDFUN(
         ido, 'I', N, task, nev, tol, A.data,
         ncv, V.data, N, iparam.data, ipntr.data, W.data, X.data, X.dim1*X.dim2,
         info
      );

      if (info) {
      if (info<0) wblog(FL,
         "ERR arpack::%s returns info=%d (ido=%d, i=%d)\n(%s)",
          AUPDSTR, info, ido, iparam[8], info2str(info));
      if (info==1) wblog(FL,
         "WRN arpack::%s returns info=%d (ido=%d, nconv=%d)\n"
         "numop: [%d %d %d %d; %d]", AUPDSTR, info, ido, iparam[4],
          iparam[8], iparam[9], iparam[10], count, iparam[2]);
      else if (info>0) wblog(FL,
         "WRN arpack::%s returns info=%d (ido=%d, i=%d)",
          AUPDSTR, info, ido, iparam[8]);
      }

      i1=(ipntr[0]-1)/N; 
      i2=(ipntr[1]-1)/N;

      if (ipntr[0]%N!=1 || ipntr[1]%N!=1) wblog(FL, 
      "ERR ipntr(1:2) does not refer to start of column !?\n"
      "[%d %d] for workd = (%dx%d)",ipntr[0]-1,ipntr[1]-1,W.dim1,W.dim2);
      if (i1>=W.dim1 || i2>=W.dim1 || i1==i2) wblog(FL,
      "ERR invalid index into WORKD (%d,%d; %d) !?",i1,i2,W.dim1);

      if (ido==1 || ido==-1) {
         wbvector<TD> xin,xout; count++;

         xin .init2ref(N, W.rec(i1)); PSI .map2VecI(FL,xin,'r');
         CP  .HZTimes(PSI,HPsi);      PSI .getDataSize(sa);
         xout.init2ref(N, W.rec(i2)); HPsi.map2Vec (FL,xout,sa);
      }
      else if (ido!=99)
      wblog(FL,"ERR %s returned ido=%d !?", AUPDSTR, ido);

   }

   EV.init(nev,N);
   EUPDFUN(1, 
     'A', select.data, 
      E.data, EV.data, N, sigma,
     'I', N, task, nev, tol, A.data, ncv,
      V.data, N, iparam.data, ipntr.data, W.data, X.data, X.dim1*X.dim2, info2
   );

   if (info2<0) wblog(FL,
      "ERR arpack::%s returns info=%d\n(%s)",
       EUPDSTR, info2, info2str(info2));
   if (info2>0) wblog(FL,
      "<i> arpack::%s returns info=%d", EUPDSTR, info2);

   if (nev>1) {
#ifdef MATLAB_MEX_FILE
      A.init2ref(N, EV.rec(1));
      PSI.map2VecI(FL,A,'r',sa); 
      PSI.put("PSI_2");
#else
      wblog(FL,
      "WRN all but the first eigenvectors are ignored (nev=%d)",nev);
#endif
   }

   A.init2ref(N, EV.rec(0));
   PSI.map2VecI(FL,A,'r',sa); 

   if (info || info2 || ncv==0) {
       wblog(FL,"WRN eigs() got info=%d, ncv=%d - take input PSI",info,ncv);
#ifdef MATLAB_MEX_FILE
       PSI0.put("psi_in","caller");
       PSI.put("psi_tmp","caller");
#endif
   }
   else {
       PSI0=PSI;
   }

   ritz.init(ncv, X.data+(ipntr[7]-1)); 
   res .init(ncv-1, X.data+(ipntr[8]-1));

   sprintf(str,"EIGS w/mode `%s' %sconverged %d/%d (%d)",
   task,iparam[4]!=(int)nev ? "NOT ":"",iparam[4],nev,N);

   } 

   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,S, "msg",   mxCreateString(str));
   mxAddField2Scalar(FL,S, "E0",    E.toMx());
   mxAddField2Scalar(FL,S, "ido",   numtoMx(ido) );
   mxAddField2Scalar(FL,S, "flag",  numtoMx(info ? info: info2));
   mxAddField2Scalar(FL,S, "N",     numtoMx(N));

   mxAddField2Scalar(FL,S, "niter", numtoMx(ipntr[14])); 

   mxAddField2Scalar(FL,S, "maxit", numtoMx(maxit));
   mxAddField2Scalar(FL,S, "ncv",   numtoMx(ncv));
   mxAddField2Scalar(FL,S, "rtol",  numtoMx(tol));
   mxAddField2Scalar(FL,S, "rnrm",  numtoMx(rnrm));
   mxAddField2Scalar(FL,S, "xflag", xflag.toMx());
   mxAddField2Scalar(FL,S, "ritz",  ritz.toMx()); 
   mxAddField2Scalar(FL,S, "relres",res.toMx());  
   mxAddField2Scalar(FL,S, "iparam",iparam.toMx());

#ifdef MATLAB_MEX_FILE
   if (info || info2 || ncv==0) {
       mxPutArray(FL,S,"Ieigs","caller");
   }
#endif

   return S;
}

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

#ifdef WB_CLOCK
   Wb::Clock clk("eigs:2Site:all",0); 
   ARG_CHECK=1;
#else
   ARG_CHECK=0;
#endif

   int i,r=-1; char ldir;
   const char isC=(typeid(_TD)==typeid(wbcomplex) ? 'C' : 0);
   double tol=-1, stol=1E-4, rnrm=1E-3;
   unsigned nc, ic1, ic2, r1, r2, k, K, nev=1, iw0, iw,
      disp=1, Nkeep=64, maxit=32,
      twoSite=0;  
   char randall;
   wbMatrix<_TQ> Qtot0, Qtot4, Qtot;
   wbvector<double> w;

   wbvector<unsigned> xflag; 

   mxArray *S=NULL, *a;

   QSpace<_TQ,_TD> Psi1, Psi2, PSI, A1, A2;
   CPAT<_TQ,_TD> CP;
   IDAT idat;
   OPTS opts;

   if (nargin) if (checkHelpVersion(argin[0])) { return; }

   if (nargin <6) wbdie(FL,"Too few input arguments."  ); else
   if (nargout<2) wbdie(FL,"Too few output arguments." ); else
   if (nargout>3) wbdie(FL,"Too many output arguments.");

   for (k=0; k<2; k++) { r=-1;
      try { mxIsQSpace(FL,argin[k],r,isC); }
      catch (...) {
         wblog(FL,"ERR argument #%d must be a valid QSpace",k+1,str);
      }
      if (mxIsEmptyQSpace(argin[k]))
      wblog(FL,"ERR argument #%d is empty QSpace !?",k+1);
   }

   Psi1.init(FL,argin[0],'r'); r1=Psi1.rank(FL);
   Psi2.init(FL,argin[1],'r'); r2=Psi2.rank(FL);

   for (; k<4; k++) if (!Mx::IsNumber(argin[k])) wblog(FL,
   "ERR argument #%d must be a valid index.", k+1);

   if (mxGetNumber(argin[2], ic1)) wblog(FL,"ERR %s",str);
   if (mxGetNumber(argin[3], ic2)) wblog(FL,"ERR %s",str);

   if (ic1 && ic1<=r1 && ic2 && ic2<=r2) {
       ic1--; ic2--;  
   }
   else wblog(FL,
   "ERR arguments #3 or #4 - dim out of bounds (%d,%d)", ic1, ic2);

   ldir=mxGetLRDir(FL,argin[k],k); k++;

   i=isCPAT(argin[k]); nc=0;
   if (i) {
      const mxArray *a = (i>1 ? mxGetCell(argin[k],0) : argin[k]);
      twoSite=isCPattern(a)-1;

      k++; nc++;

      for (; int(k)<nargin; k++)
      if (isCPAT(argin[k])) nc++; else break;
   }

   if (!nc) 
   wblog(FL,"ERR CPAT required as argument #%d.",k+1);
   for (i=k; i<nargin; i++) if (mxIsCell(argin[i]))
   wblog(FL,"invalid option as argument #%d.",i+1);

   opts.init(argin+k,nargin-k);

   opts.getOpt(FL,"NCV",NCV);

   opts.getOpt(0,0, "disp",  disp );
   wbSetLogLevel(disp>1 ? 1 : 0); 

   wblog(FL,"<%d> using %sSite method",twoSite+1, twoSite ? "two":"one");

   opts.getOpt(FL, "nev",   nev  ); 
   opts.getOpt(FL, "maxit", maxit);
   opts.getOpt(FL, "rtol",  tol  ); 
   opts.getOpt(FL, "stol",  stol );
   opts.getOpt(FL, "rnrm",  rnrm );
   opts.getOpt(FL, "Nkeep", Nkeep);
   opts.getOpt(FL, "xflag", a); if (a) xflag.init(FL,a);

   randall=opts.getOpt(FL,"rand");

   opts.checkAnyLeft(FL);

   if (stol<0) wblog(FL,"ERR Invalid stol=%g", stol);

   wbSetLogLevel(disp ? 1 : 0);

   CP.init(
     0., argin+5, nc,  
     ic1, r1, ic2, r2, twoSite, ldir
   );

   twoSiteInit(Psi1,Psi2,PSI, K,ic1,ic2, twoSite,ldir);

   if (isFlagGlobal("DEBUG")) {
      wblog(FL,"TST DEBUG is on.");
      CP.put("CP");
   }

   PSI.getDistQtot(Qtot0,w); w.max(iw0);

   try {
     #ifdef WB_CLOCK
      Wb::Clock cl2("eigs:2Site:actual",0); 
     #endif
     S=RunEigs(PSI,CP,xflag,disp,rnrm,maxit,nev,tol,randall);
   }

#ifdef MATLAB_MEX_FILE

   catch (char* istr) {
       CP.put("CP","caller"); PSI.put("PSI","caller");
       if (istr && !strcmp(istr,"\n"))
            wblog(FL, "ERR gotcha (RunEigs failed)");
       else wblog(FL, "ERR gotcha (RunEigs failed, %s)", istr);
   }
   catch (...) {
       CP.put("CP","caller"); PSI.put("PSI","caller");
       wblog(FL, "ERR gotcha (RunEigs failed ...)");
   }

#else

   catch (char* istr) { wblog(FL, "ERR RunEigs failed (%s).", istr); }
   catch (...)        { wblog(FL, "ERR RunEigs failed."); }

#endif

   a=mpsOrthoQS(
      PSI, A1, A2,
      K, Nkeep, stol, disp ? (disp==1 ? "ds" : "dS") : ""
   );

   PSI.getDistQtot(Qtot4,w); w.max(iw);
   A2.getQtot(Qtot,'u'); 

   if (S) {
   mxAddField2Scalar(FL, S, "Qtot", Qtot.toMx());
   mxAddField2Scalar(FL, S, "Qtot4",Qtot4.toMx());
   mxAddField2Scalar(FL, S, "Qtotw",w.toMx()); }

   if (disp>9 || disp>1 && !Qtot0.recEqual(iw0,Qtot4.rec(iw))) {
      defaultFmt(str, (_TQ)0, 2);

      if (Qtot0.dim1==1 && Qtot.dim1==1)
      wblog(FL," *  DQtot = [%s] -> [%s]  ",
         Qtot0.rec2Str(0,str).data,
         Qtot .rec2Str(0,str).data);
      else if (Qtot0.dim1>=1 && Qtot.dim1>=1)
      wblog(FL," *  DQtot = [%s] (%d) -> [%s] (%d) %12s\r\\",
         Qtot0.rec2Str(iw0,str).data, Qtot0.dim1,
         Qtot4.rec2Str(iw, str).data, Qtot.dim1,"");
      else wblog(FL,"WRN DQtot = %d -> %d ???  ", Qtot0.dim1, Qtot.dim1);
   }

   if (a && S)
   mxAddField2Scalar(FL, S, "Iortho", a);

   twoSiteFinal(Psi1,Psi2, A1,A2, ic1,ic2, twoSite,ldir);

   argout[0]=A1.toMx();
   argout[1]=A2.toMx();

   if (nargout>2) argout[2]=S; else
   if (S) mxDestroyArray(S);
}

