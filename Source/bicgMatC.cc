// #define __WBDEBUG__

char USAGE[] =
/* =============================================================
 * cell2matc.cc */

"   Usage: [A,tol] = bicgMatC(A, B, z, sx, ...                 \n\
                                                               \n\
      { Op1, i1 }, ...                                         \n\
      {  [], [], g, i3 }, ...                                  \n\
      { vOp, i1, f, i3 }, ...                                  \n\
      OPTS);                                                   \n\
                                                               \n\
      C++ program that does bicg() on A block to solve         \n\
      (H-z)|A>=|B> (without symmetries).                       \n\
                                                               \n\
   H|psi> is calculated in analogy to pseudo contraction scheme\n\
                                                               \n\
      -z|psi> + contractmat(A,Op1,i1) + ...                    \n\
              + sum_{p=1:3} contractmat(A,g{p}*sx{p},i3) + ... \n\
              + sum_{p=1:3} contract(A,vOp{p},i1, f{p}*sx{p},i3) + ... \n\
                                                               \n\
   OPTIONS                                                     \n\
                                                               \n\
       'tol', tol (1E-6; set min(tol) if rtol is set)          \n\
       'rtol', rtol  (relative tolerance)                      \n\
       'maxit', maxit (20)                                     \n\
                                                               \n\
   AW (C) Feb 2006                                              \n";

#include "wblib.h"

#ifndef Inf
#define Inf FP_INFINITE
#endif

class XMAP_CMat {

  public:

    XMAP_CMat () {
        memset(this, 0, sizeof(XMAP_CMat));
        z.r=1.;
    };

    ~XMAP_CMat () {};

    wbCMat *M;

    unsigned ic;           
    unsigned i1;           
    unsigned i2;           
    wbcomplex z;           

  protected:
  private:

};

   DATA sx;
   wbvector<wbCMat*> MALL;
   wbvector<XMAP_CMat> MAP;

void getCtrPattern(
   const mxArray *P,
   wbvector<XMAP_CMat> &map
){
   unsigned i,j,nz,l,ic,e=0,n,d=sx.D[0].R.dim1, N=mxGetNumberOfElements(P);

   wbcomplex z;
   wbCMat S,*M;
   DATA Sx;

   XMAP_CMat *m;

   if (N==1) {
       if (mxIsCell(P))
            if (mxGetNumber(mxGetCell(P,0), z)) wblog(FL,"ERR %s",str);
       else if (mxGetNumber(          P,    z)) wblog(FL,"ERR %s",str);

       map.init(d);

       for (i=0; i<d; i++) {
           m=map.data+i;

           m->z=z;
           m->i1=m->i2=i; 
       }

       return;
   }

   if (!mxIsCell(P) || mxGetM(P)!=1)
   wbdie(FL,"Contraction patterns must be grouped into cell vectors.");

   if (N==2) {
       M = new wbCMat[1]; M->init(mxGetCell(P,0));
       MALL.Append(M);

       ic=getFlag(mxGetCell(P,1));
       if (ic==0 || ic>2) wblog(FL,"ERR Invalid index ic (%d)", ic);

       map.initDef(d);

       for (i=0; i<d; i++) {
           m=map.data+i;

           m->M=M;
           m->ic=ic;
           m->i1=m->i2=i; 
       }

       return;
   }

   if (N==4) {
       char Sflag;
       mxArray *aS, *a3;

       aS=mxGetCell(P,0); Sflag=!mxIsEmpty(aS); 
       a3=mxGetCell(P,2);

       n=sx.D.len;

       e = !Mx::IsVector(a3) || n!=(unsigned)mxGetNumberOfElements(a3);
       i = !Mx::IsVector(aS) || n!=(unsigned)mxGetNumberOfElements(aS);
       if (Sflag) if (i || e)
       wblog(FL,"ERR Args 1 and 3 must be (cell) vectors of length %d.",n);
       if (e) wblog(FL,"ERR Arg 3 must be (cell) vector of length %d.", n);

       i=getFlag(mxGetCell(P,3));

       if (i!=3) wblog(FL,
       "ERR sx (arg 4) must be contracted onto index 3.");

       if (!Sflag) {
           wbvector<wbcomplex> gg;
           wbMatrix<wbcomplex> sz;

           if (!mxIsEmpty(mxGetCell(P,1)))
           wblog(FL,"ERR Need first two args to be [] for local term (g)");

           mxGetVector(FL,a3,gg);
           sx.sum(S,gg); S.toWbComplex(sz);

           if (sz.dim1!=d || sz.dim2!=d) wblog(FL, 
           "ERR Dimension mismatch (%dx%d; %d)",sz.dim1,sz.dim2,d);

           map.initDef(sz.nnz()); n=0;

           for (i=0; i<d; i++)
           for (j=0; j<d; j++) if (sz(i,j)!=double(0)) {
                if (n>=map.len) wblog(FL,"ERR Index out of bounds.");
                m=map.data+(n++);

                m->i1=i; m->i2=j;
                m->z=sz(i,j);
           }
       }
       else {
           wbvector<wbcomplex> ff, fs;

           ic=getFlag(mxGetCell(P,1));
           if (ic==0 || ic>2) wblog(FL,"ERR Invalid index ic (%d)", ic);

           mxCellVec2Data(aS, Sx);

           mxGetVector(FL,a3,ff);

           map.initDef(d*d); n=0;

           for (i=0; i<d; i++)
           for (j=0; j<d; j++) {

                fs=ff;
                for (nz=l=0; l<fs.len; l++) {
                    fs[l]*=sx.D[l](i,j);
                    if (fs[l]!=double(0)) nz++;
                }
                if (nz==0) continue; 

                Sx.sum(S,fs);

                M = new wbCMat[1]; (*M)=S;
                MALL.Append(M);

                m=map.data+(n++);
                m->M=M;
                m->ic=ic;
                m->i1=i; m->i2=j;  
           }
           map.Resize(n);
       }
       return;
   }

   wblog(FL,"ERR Invalid contraction pattern (%d).", N);
}

void makeUniqueMap(wbvector<XMAP_CMat> &map) {

   unsigned i,k,n=0,m=0,d,e=0,ik,kmin,imin;

   wbMatrix<unsigned> I3(map.len,3);
   wbvector<unsigned> dg,i1;
   wbperm is;
   wbcomplex z;
   XMAP_CMat *m1;
   wbCMat *M;

   for (i=0; i<map.len; i++) {
       I3(i,0)=map[i].ic;
       I3(i,1)=map[i].i1;
       I3(i,2)=map[i].i2;
   }

   I3.groupRecs(is,dg); i1.init(dg.len);

   for (m=i=0; i<dg.len; i++, m+=d) {
      d=dg[i]; i1[i]=is[m];
      if (d>1) {
         n+=(d-1); 

         m1=map.data+is[m];
         if ((*m1).M==NULL) {
            for (kmin=0, imin=is[m], k=1; k<d; k++) {
                ik=is[m+k];
                if (map.data[ik].M!=NULL) break;
                if (imin>ik) { imin=ik; kmin=k; }
            }
            if (k<d) { wblog(FL,"... %d, %d", is[m], is[m+k]); e++; break; }

            ik=i1[i]=is[m+kmin]; 
            m1=map.data+ik;
            for (k=0; k<d; k++) {
               if (k==kmin) continue;
               (m1->z)+=map.data[is[m+k]].z;
            }
         }
         else {
            for (kmin=0, imin=is[m], k=1; k<d; k++) {
                ik=is[m+k];
                if (map.data[ik].z!=double(1)) break;
                if (imin>ik) { imin=ik; kmin=k; }
            }
            if (k<d) { wblog(FL,"... %d, %d", is[m], is[m+k]); break; }

            ik=i1[i]=is[m+kmin]; 
            m1=map.data+ik;

            M = new wbCMat[1]; (*M)=*(m1->M);
            m1->M=M; MALL.Append(M);

            for (k=0; k<d; k++) {
                if (k==kmin) continue;
                M->plus(*(map.data[is[m+k]].M));
            }
         }
      }
      else
      m1=map.data+is[m];

      if (m1->M) if (m1->M->isdiag) { 
         if (!(m1->M->isDiag()))
         wblog(FL,"ERR isdiag inconsistency (%d).", m1->M->isdiag);

         m1->M->R.isdiag=1;
         if (m1->M->I.data) m1->M->I.isdiag=1;
#ifdef __WBDEBUG__
         wblog(FL,"TST got diag M entry.");
#endif
      }
   }

   if (e) { wblog(FL,
      "WRN Failed to reduce MAP to unique (%d/%d entries).", n,map.len);
       return;
   }

   map.selectSU(i1);

   for (i=0; i<MALL.len; i++) MALL[i]->user[0]=0;
   for (i=0; i<map.len;  i++) if (map[i].M) map[i].M->user[0]++;
   for (i=0; i<MALL.len; i++) if (MALL[i]->user[0]==0) MALL[i]->init();

   return;
}

DATA& ATIMES(const DATA &X, DATA &AX, const char tflag='N') {

   unsigned i, i1, i2;
   wbcomplex zt;
   XMAP_CMat *m;

   if (!Wb::strchri("NTC",tflag))
   wblog(FL,"ERR Invalid flag %c<%d>", tflag, tflag);

   AX.initSize(X);

   for (i=0; i<MAP.len; i++) {
      m=MAP.data+i;
      i1=m->i1; i2=m->i2; if (tflag!='N') SWAP(i1,i2);

      if (m->M==NULL) {
         zt=m->z; if (tflag=='C') zt.Conj();
         AX.D[i2].plus(X.D[i1], zt);
      }
      else {
         if (m->z!=1.) wblog(FL,
         "ERR z must be contained in M - will be ignored.");

         X.D[i1].Contract(m->ic, *(m->M), AX.D[i2], tflag);
      }
   }

   return AX;
}

double nrmCheckStag(const DATA &A, const DATA &B) {

    wbcomplex z1, z2;
    unsigned i,ib,s;
    double *r1,*r2,*i1,*i2, norm_inf=0.;

    if (!A.sameSize(B))
    wblog(FL,"ERR Incompatible data structures - return.");

    for (ib=0; ib<A.D.len; ib++) {
       r1=A.D[ib].R.data; i1=A.D[ib].I.data;
       r2=B.D[ib].R.data; i2=B.D[ib].I.data;

       s = A.D[ib].R.dim1 * A.D[ib].R.dim2;

       if (i1==NULL && i2==NULL) {
          for (i=0; i<s; i++) {
             if (r2[i]==0.) if (r1[i]!=0.) { norm_inf=Inf; break; }
             norm_inf=MAX(norm_inf, fabs(r1[i]/r2[i]));
          }
          if (isinf(norm_inf)) break;
       }
       else {
          for (i=0; i<s; i++) {
             z1.set(r1[i], i1 ? i1[i] : 0.);
             z2.set(r2[i], i2 ? i2[i] : 0.);

             if (z2==0.) if (z1!=0.) { norm_inf=Inf; break; }

             norm_inf=MAX(norm_inf, z1.abs()/z2.abs());
          }
          if (isinf(norm_inf)) break;
       }
    }

    return norm_inf;
}

void outPut_CtrPattern (const char* Qname="Q", const char* Mname="M") {

    unsigned i;
    const char *fnames[] = {"M","ic","i1","i2","z"};
    mxArray *C=mxCreateStructMatrix(MAP.len,1, 5,fnames);
    const char* ws="base";

    for (i=0; i<MAP.len; i++) {
       mxSetFieldToNumber(C,i,0, (double)(unsigned long)MAP[i].M);
       mxSetFieldToNumber(C,i,1, (double)MAP[i].ic);
       mxSetFieldToNumber(C,i,2, (double)MAP[i].i1);
       mxSetFieldToNumber(C,i,3, (double)MAP[i].i2);
       mxSetFieldToNumber(C,i,4, MAP[i].z);
    }

    i=mexPutVariable(ws, Qname, C); if (i==1) wblog(FL,
    "ERR Could not put variable %s into workspace %s.",Qname,ws);
    mxDestroyArray(C);

    C=mxCreateCellMatrix(MALL.len,1);
    for (i=0; i<MALL.len; i++) MALL[i]->mat2mxc(C,i);

    i=mexPutVariable(ws, Mname, C); if (i==1) wblog(FL,
    "ERR Could not put variable %s into workspace %s.",Mname,ws);
    mxDestroyArray(C);
}

void dispIfComplex( const char *file, int line,
   const char *vname,
   const wbcomplex &z
){
   if (z.i!=0.) wblog(file,line,
      "TST NB! %-5s is complex (%10.4g %+10.4gi)", vname, z.r, z.i);
}

char tstCMatProd(const wbCMat &A0, const wbCMat &B0) {

    char flag[]="NTCc",e=0;
    unsigned i,j;
    wbCMat A1,A2,C1,C2;
    double dbl,diffmax=0.;
    wbcomplex z1(0.5967, 0.0846), z2(0.5990, 0.1894);

    for (i=0; i<4; i++)
    for (j=0; j<4; j++) {
       Wb::CMatProd(A0,B0,C1,flag[i],flag[j]);

       A1=A0; A2=B0;

       if (flag[i]=='c' || flag[i]=='C') A1.conj();
       if (flag[j]=='c' || flag[j]=='C') A2.conj();

       if (flag[i]=='T' || flag[i]=='C') A1.transp();
       if (flag[j]=='T' || flag[j]=='C') A2.transp();

       Wb::CMatProd(A1,A2,C2);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in CMatProd() %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }

       Wb::CMatProd(A0,B0,C1,flag[i],flag[j],M_PI,1.);
       C2*=(1+M_PI);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in CMatProd() %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }

       Wb::CMatProd(A0,B0,C1,flag[i],flag[j],z1,z2);
       Wb::CMatProd(A1,A2,C2,'N','N',z1,z2);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in CMatProd() %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }
    }

    if (e==0)
         wblog(FL, "CMatTest passed (max rel.diff=%.4g) - ok.", diffmax);
    else wblog(FL, "ERR CMatTest failed (%d).", e);

    return e;
}

char tstDiagMatProd() {

    unsigned i,j;
    char flag[]="NTC",e=0;

    wbMatrix<double> A, B, C1, C2, A0(3,7), B0(7,11), C0(3,11);
    double afac=M_PI, cfac=sqrt(2.), dbl, diffmax=0.;

    C0.setRand();

    A0.setDiagRand();
    B0.setRand();

    A=A0; B=B0; C1=C0; A.isdiag=B.isdiag=0;
    Wb::MatProd(A, B, C1, 'N', 'N', afac, cfac);

    for (i=0; i<3; i++)
    for (j=0; j<3; j++) {
       if (flag[i]=='N') A=A0; else A0.transpose(A);
       if (flag[j]=='N') B=B0; else B0.transpose(B);

       C2=C0;
       Wb::MatProd(A, B, C2, flag[i], flag[j], afac, cfac);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in MatProd() using isdiag %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }
    }

    A0.setRand();
    B0.setDiagRand();

    A=A0; B=B0; C1=C0; A.isdiag=B.isdiag=0;
    Wb::MatProd(A, B, C1, 'N', 'N', afac, cfac);

    for (i=0; i<3; i++)
    for (j=0; j<3; j++) {
       if (flag[i]=='N') A=A0; else A0.transpose(A);
       if (flag[j]=='N') B=B0; else B0.transpose(B);

       C2=C0;
       Wb::MatProd(A, B, C2, flag[i], flag[j], afac, cfac);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in MatProd() using isdiag %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }
    }

    A0.setDiagRand();
    B0.setDiagRand();

    A=A0; B=B0; C1=C0; A.isdiag=B.isdiag=0;
    Wb::MatProd(A, B, C1, 'N', 'N', afac, cfac);

    for (i=0; i<3; i++)
    for (j=0; j<3; j++) {
       if (flag[i]=='N') A=A0; else A0.transpose(A);
       if (flag[j]=='N') B=B0; else B0.transpose(B);

       C2=C0;
       Wb::MatProd(A, B, C2, flag[i], flag[j], afac, cfac);

       dbl = C1.maxRelDiff(C2); diffmax=MAX(diffmax,dbl);
       if (dbl>1E-12) {
          wblog(FL, "WRN Inconsistency in MatProd() using isdiag %c%c:: %g",
          flag[i], flag[j], dbl); e++;
       }
    }

    if (e==0)
    wblog(FL, "SUC MatTest using isdiag passed (reldiff=%.4g)",diffmax);
    else wblog(FL, "ERR CMatTest failed (%d).", e);

    return e;
}

void RunBiCG(DATA &X, const DATA &B, OPTS &opts, IDAT &idat, double &tol) {

   unsigned iter, ibest=0, maxit=20, stag, rtflag=0;

   double dbl, dbl1=-1., dbl2=-1., normr, normrmin, nb, tolb, tolmin;

   wbcomplex z, alpha, beta, rho, rho1, ptq, gi;

   DATA AX, XX, XBEST, Q, QT, R, RT, P, PT;

   opts.getOpt("tol",   dbl1);
   opts.getOpt("rtol",  dbl2);
   opts.getOpt("maxit", maxit);

   if (dbl1>0) { tol=dbl1; rtflag=0; }
   if (dbl2>0) { tolmin=tol; tol=dbl2; rtflag=1; }

#ifdef __WBDEBUG__
   if (rtflag) {
   wblog(FL,"___ tolmin: %g", tolmin);
   wblog(FL,"___ rtol  : %g", tol   );
   } else
   wblog(FL,"___ tol   : %g", tol  );

   wblog(FL,"___ maxit : %d", maxit);
#endif

   opts.checkAnyLeft(); 

   if (((int)maxit)<0)
   wblog(FL,"ERR maxit=%d (%d)???", maxit);

   if (X.D.len) {
      if (!X.sameSize(B)) wblog(FL,"ERR Size mismatch of guess X.");
   }
   else X.initSize(B);
   XBEST=X;

   nb=sqrt(B.scalarProd(B).abs()); tolb=tol*nb;
   if (nb==0) {
       X.initSize(B); idat.init(1); idat.niter=0;
       return;
   }
   idat.init(maxit); idat.flag=1;

   B.minus( ATIMES(X,AX), R);

   normr=sqrt(R.scalarProd(R).abs());
   idat.resvec[0]=normr;

   if (rtflag) {
       dbl = MAX(tolmin, MIN(1.,tol*normr/nb));

#ifdef __WBDEBUG__
       wblog(FL,"___ rtol=%g -> tol=%.3g", tol, dbl);
#endif

       tol=dbl; tolb=tol*nb;
   }

   rho1=rho=1.; stag=0; normrmin=normr;

   if (normr<=tolb) {
       iter=ibest=0;
       idat.flag=0;
   }
   else {

       RT=R; P=R; PT=RT;

       for (iter=1; iter<maxit; iter++) {

#ifdef __WBDEBUG__
          printf("%s  %s:%d  %6d/%d",
          iter>2 ? "      \r" : "\n", FL, iter, maxit);
#endif
          rho1=rho; rho=R.scalarProd(RT);
          if (rho==0. || rho.isinf()) { idat.flag=4; break; }

          if (iter>1) {
             beta=rho/rho1;
             if (beta==0. || beta.isinf()) { idat.flag=4; break; }

             (P *=beta       ).plus(R );
             (PT*=beta.conj()).plus(RT);
          }

          ATIMES(P,Q); ATIMES(PT,QT,'C');

          ptq=Q.scalarProd(PT); if (ptq==0.) { idat.flag=4; break; }
          alpha=rho/ptq;        if (alpha.isinf()) { idat.flag=4; break; }

          if (alpha==0.) stag=1;
          if (stag==0) {
             dbl = nrmCheckStag(P,X);
             if (!isinf(dbl)) if (alpha.abs()*dbl<DBL_EPSILON) stag=1;
          }

          X.plus(P,alpha);

          B.minus( ATIMES(X,AX), XX);

          normr=sqrt(XX.scalarProd(XX).abs()); idat.resvec[iter]=normr;

          if (normr<=tolb) { idat.flag=0; break; }
          if (stag)        { idat.flag=3; break; }

          if (normr<normrmin) { normrmin=normr; XBEST=X; ibest=iter; }

          R .minus(Q, alpha);
          RT.minus(QT,alpha.conj());
       }

#ifdef __WBDEBUG__
       printf("\n\n");
#endif

   }

   if (iter<maxit)
   idat.resvec.Resize(iter+1);

   idat.niter=ibest;

   if (idat.flag==0) idat.relres=normr/nb;
   else {
      X=XBEST;
      idat.relres=normrmin/nb;
   }

   const char *fnames[] = {
     "flag", "niter", "relres", "resvec",
      "nb", "normr", "gi", "stag",
      "alpha", "beta", "rho", "rho1", "ptq", "tol" };

   gi = X.scalarProd(B);

   mxArray *S=mxCreateStructMatrix(1,1,14,fnames);

   mxSetFieldToNumber(S, 0,  0, (double)idat.flag);
   mxSetFieldToNumber(S, 0,  1, (double)idat.niter);
   mxSetFieldToNumber(S, 0,  2, idat.relres);
   idat.resvec.mat2mxs(S,    3, 1) ;

   mxSetFieldToNumber(S, 0,  4, nb);
   mxSetFieldToNumber(S, 0,  5, normr);
   mxSetFieldToNumber(S, 0,  6, gi);
   mxSetFieldToNumber(S, 0,  7, (double)stag);
   mxSetFieldToNumber(S, 0,  8, alpha);
   mxSetFieldToNumber(S, 0,  9, beta);
   mxSetFieldToNumber(S, 0, 10, rho);
   mxSetFieldToNumber(S, 0, 11, rho1);
   mxSetFieldToNumber(S, 0, 12, ptq);
   mxSetFieldToNumber(S, 0, 13, tol);

   mxPutAndDestroy(FL, S, "I", "caller");

#ifdef __WBDEBUG__

   printf("%s:%d %s %3d/%d iterations @ res=%10.4g (%d), |G|=%8.3g\n",
         __FILE__, __LINE__, idat.relres<tol ? "-. " : "ERR",
         iter, maxit, idat.relres, idat.flag, gi.abs());

   wblog(FL,"    normr = %.4g (%.4g; nb=%.4g)", normr, normrmin, nb);

   X.put("Xout"); B.put("B");
   ATIMES(X,XX);     XX.put("XX");
   ATIMES(X,XX,'C'); XX.put("XC");

   B.minus( ATIMES(X,AX), XX); XX.put("R","caller");

#else

   B.minus( ATIMES(X,AX), XX); XX.put("R","caller");

#endif

}

static void myCleanUp(void) {
    unsigned i;

    for (i=0; i<MALL.len; i++) {
        if (MALL[i]) {
           delete [] MALL[i];
           MALL[i]=NULL;
        }
    }

    MALL.init();
    MAP.init();
    sx.init();
}

static void CleanUpV(void) {  

    wblog(FL,"Clean up %s ...", mexFunctionName());
    if (MALL.len) wblog(FL,"   ~MALL (%d).", MALL.len);
    if (MAP .len) wblog(FL,"   ~MAP  (%d).", MAP .len);
    if (sx.D.len) wblog(FL,"   ~sx   (%d).", sx.D.len); fflush(0);

    myCleanUp();
}

void mexFunction(
   int nargout, mxArray *plhs[],
   int nargin, const mxArray *prhs[]
){
   unsigned i,d;
   double tol=1E-6;
   DATA A,B;
   IDAT idat;

   wbvector<XMAP_CMat> map1;
   OPTS opts;

   wbMatrix<double> D;

   if (nargin) if (isHelpIndicator(prhs[0])) { usage(); return; }

   mexAtExit(CleanUpV);

   if (nargin<4)       usage(FL,"Too few input arguments.");
   else if (nargout>2) usage(FL,"Too many output arguments.");

   A.init(prhs[0]); d=A.D.len;
   B.init(prhs[1]);

   if (!A.sameSize(B))
   wblog(FL,"ERR Dimension mismatch between A and B.");

   mxCellVec2Data(prhs[3],sx);

   if (sx.D.len!=3)
   wblog(FL,"ERR sx must be a 3D cell array of %dx%d matrices.",d,d);

   for (i=0; i<sx.D.len; i++)
   if (sx.D[i].R.dim1!=d || sx.D[i].R.dim2!=d || sx.D[i].I.data) {
       sprintf_str( "3d cell array of %s%dx%d matrices",
       sx.D[i].I.data ? "REAL " : "", d, d);
       if (sx.D[i].R.dim1!=d || sx.D[i].R.dim2!=d)
                           wblog(FL,"ERR sx must be %s.",   str);
       if (sx.D[i].I.data) wblog(FL,"WRN sx should be %s.", str);
       break;
   }

   getCtrPattern(prhs[2],MAP); 

   for (i=4; i<(unsigned)nargin; i++) {
       if (!mxIsCell(prhs[i])) break;

       getCtrPattern(prhs[i],map1);
       MAP.cat(map1);
   }

#ifdef __WBDEBUG__
   sx.put3D("sx_");
   outPut_CtrPattern();
#endif

   makeUniqueMap(MAP);

#ifdef __WBDEBUG__
   outPut_CtrPattern("Q2","M2");
#endif

   opts.init(prhs+i,nargin-i);

   RunBiCG(A,B,opts,idat,tol);

   if (nargout>0) plhs[0]=A.dat2mx3D();
   if (nargout>1) plhs[1]=numtoMx(tol);

   myCleanUp();
   return;
}

