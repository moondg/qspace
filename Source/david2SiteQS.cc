// #define __WBDEBUG__

char USAGE[] =
/* ========================================================================
 * david2SiteQS.cc */

"   david2SiteQS                                                        \n\
                                                                        \n\
     C++ program that solves H|A>=E0|A> using Abelian symmetries        \n\
     followed by an orthonormaliation step either for the `left'        \n\
     (Psi1, LR) or the `right' (Psi2, RL) QSpace.                       \n\
                                                                        \n\
     The eigenvalue problem is solved iteratively using the             \n\
     Jacobi-Davidson algorithm for obtaining the ground state.          \n\
                                                                        \n\
   1-SITE USAGE                                                         \n\
                                                                        \n\
     [Psi1, Psi2, [,info]] =                                            \n\
     david2SiteQS(Psi1, Psi2, ic1, ic2, lrdir, ...                      \n\
       { A, i1 [,opA] }, ...                           % H_{L/R, local} \n\
       {{ A, ia [,opA] }, { B, ib [,opB] }, ... }...   % H_coupling     \n\
       OPTS ...                                                         \n\
     );                                                                 \n\
                                                                        \n\
   2-SITE USAGE (NB! either i1 or i2 must be empty [])                  \n\
                                                                        \n\
     [Psi1, Psi2, [,info]] =                                            \n\
     david2SiteQS(Psi1, Psi2, ic1, ic2, lrdir, ...                      \n\
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
   Options for Davidson algorithmus                                     \n\
                                                                        \n\
    'npass',.. number of davidson passes (1)                            \n\
    'ddav', .. dimension to build for each davidson pass (4)            \n\
    'rtol', .. relative tolerance on H-residual to call converged (1E-8)\n\
    'disp',... more detailed log to stdout (0:none ... 4:full; (1))     \n\
                                                                        \n\
   Options related to symmetry search                                   \n\
                                                                        \n\
    'xflag',.. which dimension of Hilbert space to expand to full       \n\
               randomizing newly acquired coefficient blocks looking up \n\
               input contraction pattern                                \n\
    'rnrm',..  overall amplitude of randomly initialized newly acquired \n\
               data blocks (1E-3)                                       \n\
                                                                        \n\
   Options for subsequent orthonormalization                            \n\
                                                                        \n\
    'Nkeep',.. maximum dimension connecting Psi1 to Psi2 (64)           \n\
    'stol',... tolerance on SVD (i.e. on density matrix; 1E-4)          \n\
                                                                        \n\
   Test parameters                                                      \n\
                                                                        \n\
    'nev'...   number of eigenvalues to compute (1).                    \n\
    'randall'  randomize input A matrizes.                              \n\
                                                                        \n\
   OUTPUT                                                               \n\
                                                                        \n\
     Psi1(2)   updated and properly orthonormalized Psi1(2) using lrflag\n\
     info      info structure containing                                \n\
      .info    info summary of david()                                   \n\
      .svd     vector of singular values                                \n\
                                                                        \n\
   AW (C) Oct 2009                                                     \n";

#define WB_CLOCK

#define QTYPE int
#define _TQ QTYPE
#define _TD double

#include "wblib.h"
#include "mpsortho.cc"
#include "cpat.hh" 

template<class TQ, class TD>
mxArray* RunDavidson(
   QSpace<TQ,TD> &PSI0, CPAT<TQ,TD> &CP, wbvector<unsigned> &xflag,
   char vflag, double rnrm, unsigned npass, unsigned ddav,
   unsigned nev, double rtol,  char randall
);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){

#ifdef WB_CLOCK
   Wb::Clock clk("davidson:2Site:all",0); 
   ARG_CHECK=1;
#else
   ARG_CHECK=0;
#endif

   static unsigned iparam[]={0,0,0};
   static double dparam[]={0,0,0};

   int i,r=-1; char ldir;
   const char isC=(typeid(_TD)==typeid(wbcomplex) ? 'C' : 0);
   double rtol=1E-8, stol=1E-6, rnrm=1E-3;
   unsigned nc, ic1, ic2, r1, r2, k, K, nev=1, iw0, iw,
      disp=1, Nkeep=64, npass=1, ddav=4,
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

   opts.getOpt(0,0, "disp",  disp );
   wbSetLogLevel(disp>1 ? 1 : 0); 

   wblog(FL,"<%d> using %sSite method  ",twoSite+1, twoSite ? "two":"one");

   opts.getOpt(FL, "nev",   nev  ); 
   opts.getOpt(FL, "rtol",  rtol );
   opts.getOpt(FL, "stol",  stol );
   opts.getOpt(FL, "npass", npass);
   opts.getOpt(FL, "ddav",  ddav );
   opts.getOpt(FL, "rnrm",  rnrm );
   opts.getOpt(FL, "Nkeep", Nkeep);
   opts.getOpt(FL, "xflag", a); if (a) xflag.init(FL,a);

   randall=opts.getOpt(FL,"rand");

   opts.checkAnyLeft(FL);

   if (rtol<0) wblog(FL,"ERR invalid rtol=%g", rtol);
   if (stol<0) wblog(FL,"ERR invalid stol=%g", stol);

   wbSetLogLevel(disp ? 1 : 0);

   if (disp) {
      if (iparam[0]==0) { iparam[0]=npass;
         wblog(FL,
         " *  npass=%d, ddav=%d, nev=%d\n"
         " *  rtol=%g, stol=%g, rnrm=%g%N",
         npass,ddav,nev, rtol,stol,rnrm);
      }
      else { str[0]=0;
         if (checkUpdate(iparam[0], npass))
            sprintf(str+strlen(str)," npass=%d",npass);
         if (checkUpdate(iparam[1], ddav ))
            sprintf(str+strlen(str)," ddav=%d", ddav );
         if (checkUpdate(iparam[2], nev  ))
            sprintf(str+strlen(str)," nev=%d",  nev  );
         if (checkUpdate(dparam[0], rtol ))
            sprintf(str+strlen(str)," rtol=%g", rtol );
         if (checkUpdate(dparam[1], stol ))
            sprintf(str+strlen(str)," stol=%g", stol );
         if (checkUpdate(dparam[2], rnrm ))
            sprintf(str+strlen(str)," rnrm=%g", rnrm );
      }
   }

   CP.init(
     0., argin+5, nc,  
     ic1, r1, ic2, r2, twoSite, ldir
   );

   twoSiteInit(Psi1,Psi2,PSI, K,ic1,ic2, twoSite,ldir);

#ifdef MATLAB_MEX_FILE
   if (isFlagGlobal("DEBUG")) {
      wblog(FL,"TST DEBUG is on.");
      CP.put("CP");
   }
#endif

   PSI.getDistQtot(Qtot0,w); w.max(iw0);

   try {
     #ifdef WB_CLOCK
      Wb::Clock cl2("davidson:2Site:actual",0); 
     #endif
      S=RunDavidson(PSI,CP,xflag,disp,rnrm,npass,ddav,nev,rtol,randall);
   }

#ifdef MATLAB_MEX_FILE

   catch (char* istr) {
       CP.put("CP_","caller"); PSI.put("PSI","caller");
       if (istr && !strcmp(istr,"\n"))
            wblog(FL, "ERR gotcha (%s failed)",__FUNCTION__);
       else wblog(FL, "ERR gotcha (%s failed, %s)",__FUNCTION__,istr);
   }
   catch (...) {
       CP.put("CP_","caller"); PSI.put("PSI","caller");
       wblog(FL, "ERR gotcha (%s failed ...)",__FUNCTION__);
   }

#else

   catch (char* istr) { wblog(FL, "ERR %s failed (%s).",__FUNCTION__,istr); }
   catch (...)        { wblog(FL, "ERR %s failed.",__FUNCTION__); }

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

   if (disp>9 || (disp>1 && !Qtot0.recEqual(iw0,Qtot4.rec(iw)))) {
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

}; 

template<class TQ, class TD>
mxArray* RunDavidson(
   QSpace<TQ,TD> &PSI0,       
   CPAT<TQ,TD> &CP,           
   wbvector<unsigned> &xflag, 
   char vflag,                
   double rnrm,               
   unsigned npass,            
   unsigned ddav,             
   unsigned nev,              
   double rtol,               
   char randall               
){

   static unsigned check_Hconj_id=0;

   double res,tol,dbl;
   char msg[256];
   double E00=0;

   unsigned i,j,l,m,r,ip, flag=0, id=0, ncv=0, N=0;
   unsigned w=0; 
   wbvector<TD> v,Hv;
   wbvector<unsigned> sa;
   wbvector<double> E(1),Elast(1);
   wbarray<TD> H,U,V,X,EV;

   QSpace<TQ,TD> PSI, HPsi;

   Wb::SigHandler SIG(FL);

   PSI=PSI0; msg[0]=0;

   if (CP.CP.len) { 
      if (PSI.isEmpty()) {
         for (r=i=0; i<CP.CP.len; i++) { m=CP.CP[i].len;
         for (j=0; j<m; j++) r=MAX(r,CP.CP[i][j].ic); }

         PSI.init(0,r+1,CP.CP[0][0].A[0].QDIM);
         PSI.ExpandQ(CP, 1., xflag,vflag);
      }
      else {
         PSI.ExpandQ(CP,rnrm,xflag,vflag);

         if (rnrm!=0 && rtol>=0) rtol=-1;
      }
   }

   if (randall) PSI.setRand(1.);
   CP.initCounter();

   if (CP.CP.len) 
   if ((check_Hconj_id++)<12 || check_Hconj_id%23) {
      double h1,h2; QSpace<TQ,TD> RR(PSI); RR.setRand(1.);
      CP.HZTimes(RR,HPsi,'C'); h1=PSI.scalarProd(HPsi);
      CP.HZTimes(RR,HPsi    ); h2=PSI.scalarProd(HPsi);

      dbl=fabs(h1-h2)/MAX(fabs(h1),fabs(h2));
      if (dbl>1E-11) {  
         if (dbl<1E-8) {
            sprintf(msg,"WRN H only approx. Hermitian (%.3g)",dbl);
            wblog(FL,"%N%N%s",msg); flag++;
         }
         else wblog(FL, "%N%N"
         "ERR H not Hermitian (%.3g, %.3g -> %.3g) !?",h1,h2,dbl);
      }
   }

   CP.HZTimes(PSI,HPsi);

   if (PSI.QIDX.dim1!=HPsi.QIDX.dim1) {
#ifdef MATLAB_MEX_FILE
      PSI.put("PSI"); HPsi.put("HPSI");
#endif
      wblog(FL,"ERR block setting changes (%d->%d)",
      PSI.QIDX.dim1,HPsi.QIDX.dim1);
   }

   PSI.Sort(HPsi); N=PSI.getDataSize(sa); Hv.init(N);

   if (nev>N) nev=N;
   if (nev<=0) wblog(FL,
   "ERR invalid number of eigenvectors (nev=%d)",nev);

   ncv=MAX(2*nev,ddav);

   if (ncv>N) ncv=N;
   if (ncv<=nev) wblog(FL,"WRN setting ncv=%d (nev=%d)",ncv,nev);

   H.init(ncv,ncv); Elast[0]=1E32;
   U.init(N,ncv); 
   U.SIZE[1]=2;   

   v.init2ref(N, U.ref(0));
   PSI.map2Vec(FL,v,sa); v.Normalize(FL); 
   v.unRef();

   if (ddav==0 || npass==0 || N<=1 || CP.CP.len==0) { str[0]=0;
      E[0]=PSI.scalarProd(HPsi); 
   }
   else {

   CP.nr_check_ctr=2;

   if (vflag>120) printf("\n");

   for (ip=0; ip<npass; ip++) { SIG.check911();
      if (vflag>120) printf("\n");

      for (id=(ip ? nev : 1); id<ncv; id++) {
         CP.HZTimes(PSI,HPsi); HPsi.map2Vec(FL,v,sa);

         for (l=id-1, j=0; j<id; j++)
         H(l,j)=H(j,l)=v.scalarProd(U.ref(j));

         if (ip==0 && id==1) E00=H(0,0); 

         if (id>1) {
            wbarray<TD> Htmp; unsigned sh[2]={id,id};
            H.resize(wbvector<unsigned>(2,sh),Htmp);

            wbEigenS(Htmp,V,E); Wb::MatProd(U,V,EV);
            v.init(N, EV.ref(w)); 
            if (id+1<=ncv) U.SIZE[1]=id+1;
         }
         else {
            E[0]=H(0,0);
            v.init(N,U.ref(w)); 
         }

         if (vflag>120) wblog(FL,"TST %d/%d %d/%d %d: %10.7f %10.7f",
         ip+1,npass, id+1,ncv,w+1,H(0,0),E[0]);

         if (id>ncv) break;

         PSI.map2VecI(FL,v,'r',sa);
         CP.HZTimes(PSI,HPsi); HPsi.map2Vec(FL,Hv,sa);
         Hv.Plus(v,-E[w]);

         res=Hv.norm(); 
         tol=rtol*ABS(E[w]); 
         if (res<rtol && fabs(Elast[w]-E[w])<rtol) w++; 
         if (w>=nev) break; 

         Elast=E;

         PSI.map2VecI(FL,Hv,'r',sa);
         CP.HZApplyDiag(PSI,HPsi,1,E[0],'I'); HPsi.map2Vec(FL,v,sa);

         for (j=0; j<id; j++) v.GSProject(U.ref(j),1); v.Normalize(FL);
         for (j=0; j<id; j++) v.GSProject(U.ref(j),1); v.Normalize(FL);

         U.SetRef(id,v.data);
         PSI.map2VecI(FL,v,'r',sa);
      }

      if (E.len<nev || E.len!=EV.SIZE[1]) wblog(FL,
      "ERR nev length inconsistency (%d,%d/%d)",E.len,EV.SIZE[1],nev);
      if (EV.SIZE.len!=2 || EV.SIZE[0]!=N) wblog(FL,
      "ERR size inconsistency (EV: %dx%d <> %s)",N,nev,EV.sizeStr().data);

      if (ip+1<npass) {
         memcpy(U.data,EV.data,N*nev*sizeof(TD));
         for (i=0; i<E.len; i++)
         for (j=0; j<E.len; j++) H(i,j)=(i!=j ? 0 : E[i]);

         U.SIZE[1]=MAX(nev,2U);

         v.init(N, EV.ref(MIN(w,nev-1)));
         PSI.map2VecI(FL,v,'r',sa);
      }

      if (w>=nev) break; 
   }

   if (w>=nev && !flag)
   sprintf(msg,"converged at ip=%d/%d (%d/%d)",ip+1,npass,id,ncv);

   if (vflag>120) {
      printf("\n\n");
      if (msg[0]) wblog(FL,"TST msg: %s",msg);
   };

   if (nev>1) {
#ifdef MATLAB_MEX_FILE
      v.init2ref(N, EV.ref(1));
      PSI.map2VecI(FL,v,'r',sa); 
      PSI.put(FL,"PSI_2");
#else
      wblog(FL,
      "WRN all but the first eigenvectors are ignored (nev=%d)",nev);
#endif
   }

   v.init2ref(N, EV.ref(0)); E.len=nev;
   PSI.map2VecI(FL,v,'r',sa);

   PSI0=PSI;  

   } 

   mxArray *S=mxCreateStructMatrix(1,1,0,NULL);

   mxAddField2Scalar(FL,S, "msg",   mxCreateString(str));
   mxAddField2Scalar(FL,S, "flag",  numtoMx(flag));
   mxAddField2Scalar(FL,S, "oldE0", numtoMx(E00));
   mxAddField2Scalar(FL,S, "E0",    E.toMx());
   mxAddField2Scalar(FL,S, "res",   numtoMx(res));
   mxAddField2Scalar(FL,S, "rtol",  numtoMx(rtol));
   mxAddField2Scalar(FL,S, "N",     numtoMx(N));
   mxAddField2Scalar(FL,S, "niter", numtoMx(CP.HPsi_counter));
   mxAddField2Scalar(FL,S, "npass", numtoMx(npass) );
   mxAddField2Scalar(FL,S, "nev",   numtoMx(nev));
   mxAddField2Scalar(FL,S, "ddav",  numtoMx(ddav));
   mxAddField2Scalar(FL,S, "ncv",   numtoMx(ncv));
   mxAddField2Scalar(FL,S, "xflag", xflag.toMx());
   mxAddField2Scalar(FL,S, "rnrm",  numtoMx(rnrm));

   return S;
};

