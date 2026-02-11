/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : clebsch (for abelian and non-abelian symmetries)
 *
 * Copyright 2024 Andreas Weichselbaum
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * --------------------------------------------------------------------- */

#ifndef __WB_CLEBSCH_AUX_CC__
#define __WB_CLEBSCH_AUX_CC__

// See also header file clebsch.hh.
// Wb,Sep20,09 ; Wb,Dec17,14

template <class TD> inline
int CG::signFirstVal(const char *F, int L,
   const TD *d, SPIDX_T n, double eps1, double eps2
){
   SPIDX_T i; double a,x;

   for (i=0; i<n; ++i) { x=double(d[i]); a=fabs(x);
      if (a>eps1) return (x<0 ? -1 : +1);
      if (a>eps2 && F) wblog(F,L,
        "WRN %s() ignores small values (%.3g/%.3g)",FCT,x,eps2
      );
   }
   return +1;
};

template <class TD> inline
int CG::rangeSignConvention(const char *F, int L,
   TD *d, SPIDX_T n, double eps1, double eps2
){
   if (signFirstVal(F,L,d,n,eps1,eps2)<0) {
      for (SPIDX_T i=0; i<n; ++i) d[i]=-d[i];
      return -1;
   }
   return +1;
};

template <class TQ, class TD>
double CG::getSymmetryStates(const char *F, int L, const QType &q, 
   const wbvector< SPARRAY_TD > &Sp,
   const wbvector< SPARRAY_TD > &Sz,
   wbvector< SPARRAY_TD > &UK, 
   wbvector<unsigned> &dd, 
   wbvector<genRG_base<TQ,TD> > &RR,
   wbvector<double> &c2eps,  
   wbMatrix<unsigned> *iOM,  
   char vflag
){
   unsigned i,j,ip,i0=0, D,d0=0, nu=0, m,found,
        it=0, nt, r=q.sub, np=Sp.len, nz=Sz.len;
   double c2=0, r2=0; int largeD;

   wbvector< SPARRAY_TD* > uk;
   SPARRAY_TD U,x1,x2,v0,vi;
   wbMatrix<TQ> z2,JJ;
   wbvector<widx_t> dJ;
   wbvector<TD> sz;
   wbperm P;

   TD x,x12,vi2;

   TD eps=1e-8, eps2=1e-10;

#ifdef WB_CLK_SPARSE
   Wb::Clock clk("cgs:getSymStates",0); 
#endif

   if (!nz || np>nz) wblog(F_L,
      "ERR got invalid or empty Sp/Sz sets for %s (%d/%d)",
       q.toStr().data, np, nz);
   if (q.isNonAbelian(0,0,'l') && (np!=r || nz!=r)) wblog(FL,
      "ERR invalid number of Sz/Sp operators for %s (%d,%d/%d)",
       q.toStr().data, np,nz,r
   );

   D=Sz[0].dim();

   largeD=q.isLargeD(D);

   gStore.rclog(q, PFL, CG_VERBOSE>6 && largeD,
      "CGS %s() %s having D=%d",
      FCT, q.type==QT_UNKNOWN ? "general" : STR(q), D);

   U.init(D,0); 

   dd.init(D);  
   uk.init(D);

   for (it=0; it<D; ++it) {
      v0.initz(D,1,1); c2=0;
      for (; i0<D; ++i0) {
         v0.setRec(0, i0,0, 1.); if (!i0 && !it) break;

         Wb::MatProd(U,v0,x1,'C'); x=x1.norm();
         if (Wb::abs(x-1)<eps) {
            if (Wb::abs(x-1)>eps2) wblog(FL,"WRN %s() "
               "%g (%g,%g)",FCT,double(x),double(eps),double(eps2));
            continue;
         } else break;
      }

      if (i0==D) break; 
      if (it) {
         Wb::MatProd(U,x1,x2); 
         v0-=x2; v0.Normalize();
         Wb::MatProd(U,Wb::MatProd(U,v0,x1,'C'),x2);
         v0-=x2; v0.Normalize();
      }

      for (i=0; i<Sz.len; i++) {
          Wb::MatProd(Sz[i],v0,vi); if (vi.norm()<eps) continue;
          if (vi.sameUptoFac(v0)) wblog(FL,
             "ERR %s() failed to determine symmetry labels\n"
             "for starting vector (%d: %d,%d)",FCT,it+1,i0+1,i+1
          );
      }

      m=0; found=1;

      while (found) { found=0;
         for (ip=0; ip<Sp.len; ++ip) {
            Wb::MatProd(Sp[ip],v0,vi); x=vi.norm();
            if (x>eps) {
               vi*=(1/x); vi.save2(v0);
               ++found; ++m;
            }
         }
      }

      if (m && vflag && vflag!='v') wblog(FL,
      " *  applied %d Sp ops to get MW seed",m);

      SPARRAY_TD V(v0); 
      SPARRAY_TD Vi;    
      m=0; found=1;

      while (found) { found=0;
      for (ip=0; ip<Sp.len; ++ip) {

          Wb::MatProd(Sp[ip],v0,vi,'C'); 

          vi2=vi.norm2(); x=Wb::sqrt(vi2/vi.SIZE[1]);
          if (x<eps) { 
             if (x>eps2) wblog(FL,"WRN %s() got %.3g [%g %g]",
                 FCT,double(x),double(eps),double(eps2));
             continue;
          }

          Wb::MatProd(Vi,vi,x1,'C'); 
          x12=x1.norm2(); x=fabs(1-sqrt(x12/vi2));
          if (x<eps) { 
             if (x>eps2) wblog(FL,"WRN %s() got %.3g [%g %g]",
                 FCT,double(x),double(eps),double(eps2));
             continue;
          }

          if (x12) {
             vi-=Wb::MatProd(Vi,x1,x2);  
          }

          Wb::MatProd(V,vi,x1,'C'); 
          x=fabs(sqrt(x1.norm2()/vi2)); if (x>eps) {
            MXPut(FL,"q").add(U,"U").add(V,"V").add(Vi,"Vi")
               .add(vi,"vi").add(x1,"x1").add(double(x),"x");
            wblog(FL,"ERR %s() got overlap with V space (%.3g / %g)",
            FCT,double(x),double(eps));
          }
          vi-=Wb::MatProd(V,x1,x2); 

          Wb::MatProd(U,vi,x1,'C');
          x=x1.aMax(); if ((x*x)>eps) {
             MXPut(FL,"qU").add(U,"U").add(V,"V").add(vi,"vi")
               .add(Sp[ip],"Sp").add(x2,"x2");
             wblog(FL,"ERR %s() got overlap with U (%.3g / %g)",
             FCT,double(x),double(eps));
          }
          vi-=Wb::MatProd(U,x1,x2); 

          vi.SkipTinyCols(eps);

          vi.OrthoNormalizeColsQR(FL,CG_EPS1);

          if (CG_EPS1 < double(WbUtil<TD>().eps())) wblog(FL,
             "ERR %s() got eps=%.3g / %.3g (%s)",FCT,double(CG_EPS1),
             double(WbUtil<TD>().eps()),TSTR(TD)
          );

          vi-=Wb::MatProd(U, Wb::MatProd(U, vi,x1,'C'),x2); 

          vi-=Wb::MatProd(V, Wb::MatProd(V, vi,x1,'C'),x2); 
          vi-=Wb::MatProd(Vi,Wb::MatProd(Vi,vi,x1,'C'),x2); 
          vi.OrthoNormalizeColsQR(FL,CG_EPS1);          

          Vi.Cat(FL,vi,2); ++found;

          if (Vi.SIZE[1]>D) wblog(FL,
          "ERR %s() Vi space out of bounds (%s / %d)",FCT,SSTR(Vi),D);
      }
          if (found) {
             V.Cat(FL,Vi,2); Vi.save2(v0); ++m;
             if (V.SIZE[1]>D) wblog(FL,
             "ERR %s() V space out of bounds (%s / %d)",FCT,SSTR(V),D);
          }
      }

      U.Cat(FL,V,2); d0=V.SIZE[1]; dd[it]=d0;

      gStore.rclog(q, FL, CG_VERBOSE>6 && largeD,
         "  > d(%02d) =%6d @ %3.0f%% |%8d ... ",
         it+1, V.SIZE[1], 100.*U.SIZE[1]/double(U.SIZE[0]),
         U.SIZE[0]-U.SIZE[1]
      );

      WB_NEW_1(uk[it]);
      V.save2(*uk[it]);

      if (U.SIZE.len && U.SIZE[1]==D) {
         ++it; break;
      }
      else if (U.SIZE[1]>D) {
         MXPut(FL).add(U,"U").add(V,"V").add(vi,"vi").add(x1,"x1")
         .add(it+1,"it").add(ip+1,"ip").add(dd,"dd");
         wblog(FL,"ERR %s() D=%d/%d",FCT,U.SIZE[1],D);
      }
   }

   nt=it;
   if (U.SIZE[1]!=D || !nt) wblog(FL,
      "ERR %s() failed to obtain symmetry multiplets (%s/%d; %d)",
       FCT,U.sizeStr().data,D,nt);
   dd.len=nt;
   UK.init(nt);
   for (it=0; it<nt; ++it) {
      uk[it]->save2(UK[it]); WB_DELETE_1(uk[it]);
   }

#ifndef WB_SKIP_ASSERT
   Wb::MatProd(U,U,x1,'C'); 
   if (!x1.isIdentityMatrix(eps2)) {
      MXPut(FL,"q").add(Sp,"SP").add(Sz,"SZ").add(U,"U").add(UK,"UK");
      wblog(FL,"ERR %s() new space not orthogonal (%d)",FCT,nt+1);
   }
#endif

   RR.init(nt); JJ.init(nt,nz); c2eps.init(nt);

   for (i=0; i<nt; ++i) {
      genRG_base<TQ,TD> &R=RR[i];
      SPARRAY_TD &V=UK[i]; d0=dd[i];

      R.Sp.init(np);
      for (j=0; j<np; j++) { 
         Wb::MatProd(V,Wb::MatProd(Sp[j],V,x1),R.Sp[j],'C');
         c2eps[i]+=R.Sp[j].Compress(0,0,CG_SKIP_REPS); 
      }

      R.Sz.init(nz); R.Z.init(d0,nz);
      for (j=0; j<nz; j++) { 
         Wb::MatProd(V,Wb::MatProd(Sz[j],V,x1),R.Sz[j],'C');
         if (!R.Sz[j].isDiagMatrix(eps2)) wblog(FL,
            "ERR %s() got non-diagonal z-operator",FCT);
         c2eps[i]+=R.Sz[j].Compress(0,0,CG_SKIP_REPS); 

         R.Sz[j].getDiag(FL,sz);
         R.Z.setCol(j,sz);
      }

      CG::FixRational(FL,R.Z.data,R.Z.numel(),4);

      CG::findMaxWeight(q,R.Z,&R.J,&P); R.q=q;

      if (P.len && P.data[0]!=0 && vflag) {

          wblog(FL,"WRN state #1 in IREP-decomp is not MW "
            "(%s,%d; %d: %d/%d)", q==QT_UNKNOWN ? 
            "*":q.toStr('t').data, P.data[0]+1,i+1,d0,D);

          MXPut X(FL,"a"); X.add(Sp,"Sp").add(Sz,"Sz").add(UK,"UK")
           .add(dd,"dd"); if (iOM) X.add(*iOM,"M");
          X.add(D,"D").add(V,"V").add(R.Z,"Z").add(P,"P"); X.put();
      }

      if (!P.isIdentityPerm()) {
          R.Z.recPermute(P);
          V.ColPermute(P); 

          for (j=0; j<np; ++j) { R.Sp[j].MatPermute(P); }
          for (j=0; j<nz; ++j) { R.Sz[j].MatPermute(P); }

          CG::rangeSignConvention(FL,V.D.data, V.D.len);
      }

      c2=CG::FixRational(
          FL, V.D.data, V.D.len, 4, CG_SKIP_EPS1, CG_SKIP_EPS2
      );
      r2+=c2; c2eps[i]+=c2;

      nu+=V.D.len;

      r2+=V.Compress(FL,CG_SKIP_EPS2); 

      JJ.recSetP(i,R.J.data);

#ifdef CG_CHECK_MW_PERM
      R.P0=P;
#endif
   }

   r2/=nu; { double r=sqrt(r2);
      if (r<=CG_EPS2) gStore.rclog(q, FL, CG_VERBOSE>6 && largeD,
         "CGS ok @ e=%.3g",r);
      else gStore.rclog(q, FL, CG_VERBOSE>6 && largeD,
         "WRN e=%.3g/%.3g",r,CG_EPS2
      );
   }

   z2=JJ; z2.groupRecs(P,dJ);
   if (dJ.anyGT(1)) {
      if (iOM) {
         unsigned d, *iom; wperm_t *p=P.data;
         iOM->init(P.len,2);

         for (i=0; i<dJ.len; ++i) { d=dJ[i];
         for (j=0; j<d; ++j, ++p) { iom=iOM->rec(*p); iom[0]=j; iom[1]=d; }}
      }
      else { 
         wblog(FL,"NB! %s() got outer multiplicity (OM<=%d)",FCT,dJ.max());
      }
      return r2;
   }
   else if (iOM) iOM->init();

   return r2;
};

 namespace CG {

template<class TQ>
int FileLock::init(
   const QType &q_, const qset<TQ> &J1, const qset<TQ> &J2,
   const char *ext 
){
   int i=0; QSet<TQ> Q;

   if (!J1.len || J1.len!=q_.qlen() || J1.len!=J2.len) wblog(FL,
      "ERR %s() invalid qset for %s\ngot [J1;J2]=[%s; %s]",
      FCT,STR(q_),STR(J1),STR(J2));
   q=q_;

   if (J1<=J2)
        { Q.init2(q,J1,J2); }
   else { Q.init2(q,J2,J1); }

   i=gStore.get_file_name(FL,fname,Q,ext,RC_SAVE_S);

   i=GetLock(FL,'W');

   if (i<=0) { wblog(FL,"ERR %s() " 
      "failed to create lock (e=%d)\n%s",FCT,i,Wb::repHome(fname).data);
   }

   return i;
};

template<class T>
double FixRational(const char *F, int L,    
   T *d, SPIDX_T n,
   unsigned niter QS_UNUSED_VAR, 
   double eps1, double eps2
){

   double r2=0;
   if (!WbUtil<T>::isInt()) {
      double dx, dxmax=0; T q;

      for (SPIDX_T i=0; i<n; ++i) {
         q=round(d[i]); dx=double(fabs(d[i]-q));
         if (dx<eps1) {
            if (dx<eps2) { r2+=dx*dx; d[i]=q; }
            else if (dxmax<dx) { dxmax=dx; }
         }
      }

      if (dxmax) wblog(FL,
         "WRN %s() got small CGCs @ %.3g!",FCT,dxmax);

      if ((dx=std::sqrt(r2))>CG_EPS2) {
         wblog(F_L,"WRN skipped %.2g (%d)",dx,n);
      }
   }

   return r2;
};

template<>
double FixRational(const char *F, int L, 
   double *d, SPIDX_T n, unsigned niter,
   double eps1, double eps2
){
   double r2z,r2a,ra,r2;

   r2=Wb::FixRational(F_L,d,n,
      'r',   
      niter, 
      1024,  
      eps1>0 ? eps1 : CG_SKIP_DEPS1,
      eps2>0 ? eps2 : CG_SKIP_DEPS2,
      &r2z,  
      &r2a,  
      0      
   );

   if ((ra=std::sqrt(r2a))>CG_SKIP_DEPS1) { 
      wblog(F_L,"WRN %s() skipped %.2g [%.2g, %.2g; %d @ %3.g]",
      FCT, ra, std::sqrt(r2z/n),std::sqrt(r2/n),CG_EPS2,n);
   }

   return r2;
};

}; 

template<class TQ>
size_t DY::wdim_A(
   unsigned n,  
   const TQ *qs 
){
   unsigned i=0, j; long x,y, xa=1, ya=1;
   TQ rho=1; 

   for (; i<n; ++i) { x=y=0;
      for (j=i; j<n; ++j) {
         x+=(qs[j]+rho); y+=1;
         xa*=x; ya*=y;
      }
   }

   if (xa<=0 || ya<=0 || (x=xa/ya)*ya!=xa) wblog(FL,
      "ERR %s() got invalid d=%ld/%ld = %.3f",FCT,xa,ya,xa/double(ya));
   return size_t(x);
};

template<class TQ>
size_t DY::wdim_C(
   unsigned n,  
   const TQ *qs 
){

   unsigned i=1, j; long x,y, xa=1, ya=1;
   TQ rho=1; 

   if (!n) wblog(FL,"ERR %s() got n=%d",FCT,n);

   for (; i<=n; ++i) { x=y=0;
      for (j=i; j<=n; ++j) {
         x+=(j<n ? 1 : 2)*(qs[j-1]+rho); y+=(j<n ? 1 : 2);
         xa*=x; ya*=y;
      }
      for (j=n-1; j>=i; --j) {
         x+=(qs[j-1]+rho); y+=1;
         xa*=x; ya*=y;
      }
   }

   if (xa<=0 || ya<=0 || (x=xa/ya)*ya!=xa) wblog(FL,
      "ERR %s() got invalid d=%ld/%ld = %.3f",FCT,xa,ya,xa/double(ya));

   return size_t(x);
};

template<class TQ>
size_t DY::wdim_B(
   unsigned n,  
   const TQ *q  
){
   unsigned i,k=0,l=1; long P, Q, P_,Q_;
   if (!n) wblog(FL,"ERR %s() got n=%d",FCT,n);

   i=n-1; P_=P=(q[i]+1); Q_=Q=1; 
   for (--i; i<n; --i) {
      P+=(q[i]+1)*2; Q+=2; 
      P_*=P; Q_*=Q; l+=1;
   }

   for (k=1; k<n; ++k) {
      for (P=Q=0, i=k-1; i<n; ++i) {
          P+=(q[i]+1);   
          Q+=1;          
          P_*=P; Q_*=Q; l+=1;
      }
      for (i-=2; i>=k; --i) {
          P+=(q[i]+1);  
          Q+=1;         
          P_*=P; Q_*=Q; l+=1;
      }
   }

   if (l!=n*n) wblog(FL,
      "ERR %s() got l=%d/%d (having n=%d)",FCT,l,n*n,n);
   else {
      double x=double(P_)/double(Q_);
      if (x<=0 || fabs(x-round(x))>1e-6)
         wblog(FL,"ERR %s() got d = %d/%d = %g",FCT,P_,Q_,x);
   }

   return size_t(P_/Q_);
};

template<class TQ>
size_t DY::wdim_D(
   unsigned n,  
   const TQ *q  
){
   unsigned i,k=0,l=1; long P, Q, P_,Q_;
   if (n<2) wblog(FL,"ERR %s() got n=%d",FCT,n);

   i=n-2; P_=P=(q[i]+1); Q_=Q=1; 
   i=n-1; P=(q[i]+1); Q=1; P_*=P; Q_*=Q; l+=1;
   for (i=n-3; i<n; --i) {
      P+=(q[i]+1); Q+=1; P_*=P; Q_*=Q; l+=1;
   }

   for (k=1; k<n-1; ++k) {
      for (P=Q=0, i=k-1; i<n; ++i) {
          P+=(q[i]+1);   
          Q+=1;          
          P_*=P; Q_*=Q; l+=1;
      }
      if (n>3) { for (i-=3; i>=k; --i) {
          P+=(q[i]+1);  
          Q+=1;         
          P_*=P; Q_*=Q; l+=1;
      }}
   }

   if (l!=n*(n-1)) wblog(FL,
      "ERR %s() got l=%d/%d (having n=%d)",FCT,l,n*(n-1),n);
   else {
      double x=double(P_)/double(Q_); 
      if (x<=0 || fabs(x-round(x))>1e-6)
         wblog(FL,"ERR %s() got d = %d/%d = %g",FCT,P_,Q_,x);
   }

   return size_t(P_/Q_);
};

namespace DY {

Symmetry& Symmetry::init(const QType &q) {

   unsigned i=1;

   if (q.type==QT_SUN) { 
      if (q.sub<1) wblog(FL,"ERR %s() got symmetry %s",FCT,STR(q));
   } else
   if (q.type==QT_SpN) { 
      if (q.sub<2) wblog(FL,"ERR %s() got symmetry %s",FCT,STR(q));
   }
   else wblog(FL,
      "ERR %s() symmetry %s not implmented yet",FCT,STR(q));
   if (q.sub>32) wblog(FL,"WRN %s() got symmetry %s",FCT,STR(q));

   n=q.sub; n2.init2val(n,1);

   A.init(n,n); A(0,0)=2;
   for (; i<n; ++i) {
      A(i,i)=2; 
      A(i-1,i)=A(i,i-1)=-1;
   }

   R.init(n,n); R(0,0)=1; 
   for (i=1; i<n; ++i) {
      R(i-1,i)=-sqrt(double(i  )/double(2*(i+1)));
      R(i,  i)= sqrt(double(i+2)/double(2*(i+1)));
   }

   if (q.type==QT_SUN) {
      unsigned i2,j, k=0, l=0;
      M.init(n,(n*(n+1))/2); 
      for (k=1; k<=n; ++k) { 
         for (i2=n-k, i=0; i<=i2; ++i, ++l) {
            for (j=0; j<k; ++j) { M(i+j,l)=1; }
         }
      }
   }
   else if (q.type==QT_SpN) {
      unsigned i2,j,p, k=0, l=0;
      A(n-1,n-2)=-2; n2[n-1]=2;

      R(n-2,n-1)=-2*double(n-1)/sqrt(double(2*n*(n-1)));
      R(n-1,n-1)= sqrt(double(2)/double(n));

      M.init(n,n*n); 
      for (k=1; k<=n; ++k) { 
         for (i2=n-k, i=0; i<=i2; ++i) {
            for (j=0; j<k; ++j) { M(i+j,l)=1; }; ++l;
         }
         for (i=i2, p=2; p<=k; ++p) { 
            for (j=0; j<k-p; ++j) { M(i+j,l)=1; }
            for (; j<k; ++j) { M(i+j,l)=2; };
            M(i+j-1,l)=1; ++l;
         }
      }
   }

   return *this;
};

template <class TQ>
size_t Symmetry::getWeightsFT(const char *F, int L,
   const qset<TQ> &qm, 
   Weights<TQ> &R  
){
   size_t d=0; unsigned i,ip,np, iter=0, niter=99;
   map < qset<TQ>, const weight_info<TQ>* > x2;
   const TQ* md; int l;

   if (!qm.len || qm.len!=q.qlen()) wblog(F_L,
      "ERR %s() invalid qset %s (%s)",FCT,STR(q),STR2(qm,q));
   if (qm.len!=A.dim2 || A.dim1!=A.dim2) wblog(F_L,
      "ERR %s() invalid Cartan matrix A (%dx%d / %d; %s)",
      FCT,A.dim1,A.dim2,qm.len,STR(q));
   if (M.SIZE.len!=2 || M.SIZE[0]!=n || M.SIZE[1]<n) wblog(F_L,
      "ERR %s() invalid M data (%s /%d)",FCT,SSTR(M),n);

   qset<TQ> qs(n); np=M.SIZE[1];
   R.init(qm);

   while (++iter<niter) {
      for (auto Il=R.X.begin(); Il!=R.X.end(); ++Il) { if (Il->second) {
         weight_info<TQ> &w0=R.W[Il->first];

         if (w0.p.len!=n || !w0.m) wblog(F_L,
            "ERR %s() got missing weight (%s) @ m=%d",
            FCT, Il->first.wbvector<TQ>::toStr(2).data, w0.m);

         for (md=M.data, ip=0; ip<np; ++ip, md+=n) {

            qs.init(Il->first); 
            for (i=0; i<n; ++i) { if (md[i]) qs.Plus(A.rec(i),-md[i]); }

            if (!qs.anyLT(0)) { 
               weight_info<TQ> &w=R.W[qs];
               if (!w.m) {
                  w=w0; for (i=0; i<n; ++i) { w.p[i]+=md[i]; }
                  ++x2[qs];

                  w.m=-1;
               }
            }
         }
      }}
      R.X.clear(); x2.swap(R.X); 
      if (!R.X.size()) {
         break;
      }
   }
   if (iter>=niter) {
      wblog(F_L,"ERR %s() got iter=%d",FCT,iter);
   }
   else if (F) wblog(F_L," *  %s() got iter=%d",FCT,iter);

   map<int,   
      map<qset<TQ>,   
         weight_info<TQ>*  
      >
   > Ml;

   for (auto I=R.W.begin(); I!=R.W.end(); ++I) {
      l=I->second.p.sum();
      if (I->first.anyLT(0)) wblog(F_L,"ERR %s() got w=[%s] "
         "outside dominant Weyl chamber",FCT,STR(I->first));
      if (!I->second.p.len || !I->second.m) wblog(F_L,
         "ERR %s() got w=[%s] with p=[%s], m=%d", FCT,
         STR(I->first), STR(I->second.p), I->second.m);
      Ml[l][I->first]=&I->second;
   }

   for (auto Il=Ml.begin(); Il!=Ml.end(); ++Il) {
      auto &M2=Il->second;
      for (auto I2=M2.begin(); I2!=M2.end(); ++I2) {
         weight_info<TQ> &w=*I2->second;
         if (w.m<0)
         w.m=getInnerMultFT(I2->first,w,R);
      }
   }

   if (F) {
      wblog(F_L,"==> %s[%s]: q=(%s) ",FCT,STR(q),STR2(qm,q));
      R.print(F_L,q); 
   }
   wExpand(F_L,R); if (F) R.print(F_L,q);

   d=R.dim();
   i=q.wdim(qm.data); 
   if (d!=i) wblog(F_L,
      "WRN %s() got %s (%s) @ d=%d/%d",FCT,STR(q),STR2(qm,q),d,i);

   return d;
};

template <class TQ>
int Symmetry::getInnerMultFT(
   const qset<TQ> &qs, const weight_info<TQ> &w,
   Weights<TQ> &R, char vflag
){
   int m=1, mk, l=w.p.nnz();

   if (!n || qs.len!=n || w.p.len!=n || n2.len!=n) wblog(FL,
      "ERR %s() got qs=(%s), p=(%s), %d/%d",
      FCT,STR(qs),STR(w.p),n,n2.len);
   if (M.SIZE.len!=2 || M.SIZE[0]!=n || M.SIZE[1]<n) wblog(FL,
      "ERR %s() got M (%s), n=%d",FCT,SSTR(M),n);
   if (w.p.anyLT(0)) wblog(FL,"ERR %s() got p=(%s)",FCT,STR(w.p));
   if (!R.W.size()) wblog(FL,"ERR %s() got empty W",FCT);

   if (l<=1) { return m; }

   unsigned i, ip=0, np=M.SIZE[1];
   qset<TQ> ap,qk; const TQ *md=M.data;

   double P=0, Q=0;

   qk=qs; qk+=R.qm; qk+=2;
   for (i=0; i<n; ++i) {
      Q+=qk[i]*w.p[i]*n2[i]; 
   }

   for (ip=0; ip<np; ++ip, md+=n) {
      ap.init(n); qk=qs;
      for (i=0; i<n; ++i) { if (md[i]) { ap.Plus(A.rec(i),md[i]); }}
      if (vflag>1) printf("  ip=%d/%d, a=[%s]: ",ip+1,np,STR(ap));
      while (1) { 
         qk.Plus(ap.data);
         if (!(mk=Get(qk,R.W,vflag))) { break; }
         for (i=0; i<n; ++i) {
            P+=mk*qk[i]*md[i]*n2[i]; 
         }
      }
      if (vflag>1 && ip+1<np) printf("\n");
   }

   if (Q) { m=2*P/Q; if (m*Q!=2*P)
      wblog(FL,"ERR %s() got m=2*%g/%g = %g !?",FCT,P,Q,2*P/Q);
   } else wblog(FL,"ERR %s() got m=2*%g/%g (div/0)",FCT,P,Q);

   if (vflag>1) printf(": (%s) => m=2*%g/%g = %g\n",STR(qs),P,Q,2*P/Q);

   return m;
};

template <class TQ>
int Symmetry::Get(
   const qset<TQ> &qk, map < qset<TQ>, weight_info<TQ> > &W,
   char vflag
){
   auto it0 = W.find(qk);

   if (it0!=W.end()) {
      int m=it0->second.m;
      if (m>0) {
         if (vflag>1) printf(" => [%s]^%d",STR(qk),m);
         return m;
      }
   }

   if (!qk.anyLT(0)) { 
      if (it0==W.end()) {
         if (vflag>1) printf("  # [%s]",STR(qk));
         return 0;
      }
      const weight_info<TQ> &w0=it0->second;
      if (w0.m<=0) wblog(FL,"ERR %s() "
         "missing inner multiplicity value m=%d [%s]!?",FCT,w0.m,STR(qk));
      else {
         if (vflag>1) printf(" => [%s]^%d",STR(qk),w0.m);
         return w0.m;
      }
   }

   if (!n || qk.len!=n) wblog(FL,
      "ERR %s() got empty qk=[%s] /%d",FCT,STR(qk),n);
   if (n!=A.dim1 || n!=A.dim2) wblog(FL,
      "ERR %s() invalid Cartan matrix %s /%d",FCT,SSTR(A),n);

   unsigned i, iter=0;
   qset<TQ> qs(qk); TQ *qd=qs.data;

   while (1) {
      for (i=0; i<n; ++i) { if (qd[i]<0) break; }
      if (i<n) { qs.Plus(A.rec(i),-qd[i]); }
      else { break; }

      if (++iter>30) wblog(FL,"ERR %s() got iter=%d",FCT,iter);
   }

   auto it = W.find(qs);
   if (it!=W.end()) {
      const weight_info<TQ> &w = it->second;
      if (w.m<=0) wblog(FL,"WRN %s() "
         "missing inner multiplicity value m=%d [%s]!?",FCT,w.m,STR(qs));
      else {
         weight_info<TQ> &wk=W[qk];
         wk.m=w.m;
      }
      return w.m;
   }
   return 0; 
};

template <class TQ>
int Symmetry::wExpand(const char *F, int L, Weights<TQ> &R) const {

   unsigned i,iter=0, niter=999; int m0, nX=0;
   qset<TQ> qs; TQ *qd;

   map <
      qset<TQ>,  
      const weight_info<TQ>* 
   > W0,X,X_;

   if (!n || n!=A.dim1 || n!=A.dim2) wblog(F_L,
      "ERR %s() invalid Cartan matrix %s /%d",FCT,SSTR(A),n);

   for (auto Iw=R.W.begin(); Iw!=R.W.end(); ++Iw) { if (Iw->second.m>0) {
      W0[Iw->first]=&(Iw->second);
   }}

   for (auto Iw=W0.begin(); Iw!=W0.end(); ++Iw) {
      X.clear(); X[Iw->first]=Iw->second; m0=Iw->second->m; iter=0;
      while (++iter<niter) {
         for (auto Ik=X.begin(); Ik!=X.end(); ++Ik) { if (Ik->second) {
            if (Ik->second->m<=0) wblog(F_L,
               "ERR %s() [%s] got m=%d",FCT,STR(Ik->first),Ik->second->m);
            qs.init(Ik->first); qd=qs.data;

            for (i=0; i<n; ++i) { if (qd[i]>0) {
               qs.Plus(A.rec(i),-qd[i]);
               weight_info<TQ> &wi=R.W[qs]; if (wi.m<=0) {
                  wi.m=m0;
                  X_[qs]=&wi; ++nX; 
               }
               qs.Plus(A.rec(i),-qd[i]); 
            }}
         }}
         X.clear(); X.swap(X_); 
         if (!X.size()) break;
      }
      if (iter>=niter) wblog(F_L,
         "ERR %s() got iter=%d/%d",FCT,iter,niter);
   }
   return nX;
};

template <class TQ>
void Weights<TQ>::print(
   const char *F, int L, const QType &q, char vflag) const {

   const unsigned flen=64;
   size_t n=W.size(); int l, isn;
   char fmt[flen];
   unsigned mtot=0, merr=0, ndom=0;

   map<int,   
      map<qset<TQ>,   
         const weight_info<TQ>*  
      >
   > M;

   for (auto I=W.begin(); I!=W.end(); ++I) {
      if (I->second.m>0) { mtot+=I->second.m; } else { ++merr; }
      isn=I->first.anyLT(0); if (!isn) { ++ndom; }
      if (vflag || !isn) {
         l=I->second.p.sum(); if (!l && !I->second.p.len) l=-1;
         M[l][I->first]=&I->second;
      }
   }

   if (F) {
      char s[128];
      snprintf(s,128,"%s: total of %d + %d = %d weight%s (%d state%s)",
         STR(q), ndom, n-ndom, n, n!=1 ? "s":"", mtot, mtot!=1 ? "s":"",
         merr ? ", incomplete!":"");
      if (!merr) wblog(F_L," *  %s",s); else wblog(F_L,"WRN %s",s);
   }
   else PRINTF("\n");

   if (n) { n=3*W.begin()->first.len; 
      snprintf(fmt,flen,"\n  level  %%-%ds %%-%ds     m    # weights",n+3,n+3);
      PRINTF(fmt,"dcoeffs","pcoeffs"); 
      snprintf(fmt,flen,"\n  %%5d%%1s (%%%ds ) [%%%ds ] %%5d",n,n);
   }
   else fmt[0]=0;

   for (auto Il=M.begin(); Il!=M.end(); ++Il) {
      auto &M2=Il->second; n=M2.size();
      unsigned i=0;

      for (auto I2=M2.begin(); I2!=M2.end(); ++I2) {
         const weight_info<TQ> &w=*I2->second;

         PRINTF(fmt, Il->first,
            X.find(I2->first)==X.end() ? "" : "*",
            I2->first.wbvector<TQ>::toStr(2).data,
            w.p.wbvector<TQ>::toStr(2).data, w.m
         );

         if (++i==1) { PRINTF("   %4d",n); }
      }
   }
   PRINTF("\n\n");
};

}; 

template <class TD>
double cdata<TD>::SkipTiny(const char *F, int L){

   double r2=CG::FixRational(F_L, this->D.data, this->D.len, 12);
   r2+=this->Compress(F_L,CG_SKIP_EPS2);
   return r2;
};

template <class TQ, class TD>
double genRG_base<TQ,TD>::SkipTiny(const char *F, int L) {

    double e2=
        CG::FixRational(F_L,Z.data,Z.numel(),12);

    return e2;
};

template <class TQ, class TD>
double CData<TQ,TD>::SkipTiny(const char *F, int L){

   double e2=cgd.SkipTiny(F_L), e=std::sqrt(e2);
   if (e>CG_EPS1) wblog(F_L,"WRN %s() skipped %.3g",FCT,e);

   return e2;
};

template <class TQ>
SPIDX_T CG::findMaxWeight( 
   const QType &q, const wbMatrix<double> &Z, qset<TQ> *J, wbperm *P_
){
   unsigned i, r=(q.type ? q.qlen(): Z.dim2); SPIDX_T k;
   wbMatrix<double> z2(Z);
   wbperm P;

   if (Z.isEmpty()) wblog(FL,"ERR %s() "
      "got empty z-labels (%dx%d; %s)",FCT,Z.dim1,Z.dim2,STR(q));
   if (r && Z.dim2!=r) wblog(FL,
      "ERR %s()\ninconsistent z-labels for %s (%dx%d; %d)",
       FCT,STR(q),Z.dim1,Z.dim2,r);

   z2.FlipCols(); 
   z2.sortRecs_float(P,-1); 
   k=P[0];

   if (J) { 
      qset<double> qm(Z.dim2,Z.rec(k)); 

      if (q.type==QT_SUN) {

         if (r<1 || r>9) wblog(FL,
            "ERR %s() got symmetry %s",myname,STR(q)); 

         for (i=r-1; i>0; --i) {
            qm.data[i]=num2int(FL, (qm.data[i] - qm.data[i-1]) / double(i+1));
         }; qm.data[i]=num2int(FL, (qm.data[i] ));
      }
      else if (q.type==QT_SpN) {

         if (r<2 || r>9) { wblog(FL,
            "ERR %s() got symmmetry %s%s",myname,STR(q), 
             r==1? "\nhint: use SU(2) instead":"");
         }

         for (i=r-1; i>0; --i) {
            qm.data[i]=num2int(FL, (qm.data[i] - qm.data[i-1]) / double(i+1));
         }; qm.data[i]=num2int(FL, (qm.data[i] ));
      }
      else if (q.type==QT_SON) {

         if (r<2 || r>9) { wblog(FL,
            "ERR %s() got symmetry %s%s",myname,STR(q), 
            r==1? "\nhint: use SU(2) instead":"");      
         }

         unsigned j,l=(r-1)/2; TQ x=num2int(FL, 2*qm.data[0]);
         for (i=1; i<r; ++i) {
            qm.data[i-1]=num2int(FL, qm.data[i] - qm.data[i-1]);
         }; qm.data[i-1]=x;

         for (i=0; i<l; ++i) { j=r-i-2;
            x=qm.data[i]; qm.data[i]=qm.data[j]; qm.data[j]=x;
         }
      }
      else if (q.type==QT_SEN) {

         if (r<=2 || r>9) { wblog(FL,
            "ERR %s() got symmetry %s%s",myname,STR(q), 
            r==1? "\nhint: use U(1) instead" : (        
            r==2? "\nhint: use SU(2)xSU(2) instead":    
            ""));

         }

         unsigned j,l=(r-1)/2; TQ x=num2int(FL, qm.data[0] + qm.data[1]);
         for (i=1; i<r; ++i) {
            qm.data[i-1]=num2int(FL, qm.data[i] - qm.data[i-1]);
         }; qm.data[i-1]=x;

         for (i=0; i<l; ++i) { j=r-i-2;
            x=qm.data[i]; qm.data[i]=qm.data[j]; qm.data[j]=x;
         }
      }

      J->initT(FL,qm); 
   }

   for (i=1; i<z2.dim1; ++i) { if (z2.recDiff2(0,i)>1e-8) break; }
   if (i>1) { 
      MXPut(FL,"ans").add(i+1,"i").add(z2,"z2").add(Z,"Z").add(P,"P");
      wblog(FL,"ERR %s()\nmaximum weight state not unique (%d)",FCT,i);
   }

   if (P_) P.save2(*P_);

  return k;
};

template <class TQ, class TD>
void CG::get_SU2mat(
   TQ s2,
   wbvector<double> &sz, 
   wbsparray<TD> &Sp, wbsparray<TD> &Sz,
   wbsparray<TD> &S2, wbsparray<TD> &E
){
   QType q("SU2");

   { double q=s2;
     if (q!=double(int(q)) || q<0 || q>2E3) wblog(FL,
       "ERR invalid spin S=%g (@ %.3g)",q/2,q-round(q)
     );
   }

   SPIDX_T i, D=SPIDX_T(s2+1);
   TD z, s=s2/2.0; 
   wbsparray<TD> X1,X2,Sz2;

   E.initIdentity(D);
   Sz.initDiag(D); Sz2.initDiag(D); Sp.initz(D,D,D-1); sz.init(D);

   for (i=0; i<D; i++)
   for (z=+s, i=0; i<D; i++, z-=1) {
       sz[i]=double(Sz[i]=z); Sz2[i]=z*z; if (i) {
       Sp.setRec(i-1, i-1, i, sqrt(s*(s+1)-z*(z+1))); }
   }

   Wb::MatProd(Sp,Sp,X1,'N','C');
   Wb::MatProd(Sp,Sp,X2,'C','N'); X1.Plus(FL,X2)*=0.5; X1.Plus(FL,Sz2);

   X1.save2(S2);
};

gTQ get_qtot_abelian(
   const QType &t, const gTQ *q0, const widx_t n,
   const widx_t *idx, const widx_t stride
){
   gTQ q=0; 
   unsigned i=0;

   switch (t.type) {

     case QT_U1 : 
     case QT_ZN : 

       if (idx)
            { for (; i<n; ++i) q+=q0[idx[i]*stride]; }
       else { for (; i<n; ++i) q+=q0[    i *stride]; }
       if (t.type==QT_ZN) q%=t.sub;
       return q;

     case QT_P : q=+1; 

       if (idx)
            { for (; i<n; ++i) q*=q0[idx[i]*stride]; }
       else { for (; i<n; ++i) q*=q0[    i *stride]; }
       return q;

     default: 
       wblog(FL,"ERR %s() got non-abelian symmetry '%s'\n"
       "(n=%d, stride=%d)",FCT,STR(t),n,stride);
   }

   wblog(FL,"ERR %s() ",FCT); 
   return q;
};

template<class TD>
int BuildKrylovH(
   const QType &q,  
   wbarray<TD> &HK, unsigned nk, 
   SPARRAY_TD &X,        
   wbvector<TD> &EK,
   const SPARRAY_TD &HL,            
   const SPARRAY_TD &HR,            
   const wbvector< SPARRAY_TD > SL, 
   const wbvector< SPARRAY_TD > SR,
   double &xmin_
){
   unsigned rval=0, ik=0, d=HL.dim();
   SPARRAY_TD HX, U;
   wbarray<TD> UK;
   wbvector<TD> v;
   TD x=0; double xmin=1e99;

   if (nk<2) wblog(FL,"ERR %s() got invalid nk=%d",FCT,nk);

   if (nk>d) nk=d; 
   HK.init(nk,nk); X.Normalize();

   for (; ik<nk; ++ik) {
      GetHPsi(HX,X,HL,HR,SL,SR);
      HK(ik,ik)=X.dotProd(FL,HX);
      if (ik) {
         TD d=Wb::abs(HK(ik-1,ik-1))+Wb::abs(HK(ik,ik));
         if (d+x==d) { rval=2; 
            HK.Resize(ik,ik);  

            gStore.rclog(q, FL, CG_VERBOSE>6," *  %s() "
              "converged @ %.4g (ik=%d/%d)",FCT,double(x),ik,nk);
            break;
         }
      }

      if (!ik){
         X.save2(U); 
      }
      else { U.Cat(FL,X,3); }

      if (ik+1==nk) { break; }

      for (unsigned io=0; io<2; ++io) { 
         U.contract(FL,"12*",HX,"12",X); X.toFull(v);
         U.contract(FL,2,v,X);
         HX-=X;  
      }

      x=HX.norm(); if (xmin>double(x)) { xmin=double(x); }
      if (x>0) {
         HX*=(TD(1)/x); HK(ik+1,ik)=HK(ik,ik+1)=x;
      }
      else {
         HK.Resize(ik+1,ik+1);  
         gStore.rclog(q, FL, CG_VERBOSE>6," *  %s() "
            "converged @ %.4g (ik=%d/%d)",FCT,double(x),ik+1,nk);
         rval=3; break;
      }

      HX.save2(X);
   }

   HK.eigTriDiag(FL,EK,&UK,CG_EPS2);
   UK.getCol(0,v);
   U.contract(FL,2,v,X); 

   if (xmin<CG_EPS1 && xmin>xmin_) rval=1;
   xmin_=xmin;

   return rval;
};

template<class TD>
void GetHPsi(
   SPARRAY_TD &HX, const SPARRAY_TD &X,
   const SPARRAY_TD &HL, 
   const SPARRAY_TD &HR, 
   const wbvector< SPARRAY_TD > SL, 
   const wbvector< SPARRAY_TD > SR  
){
   SPARRAY_TD Q; HX.init();
   TD one(1);

   if (SL.len!=SR.len) wblog(FL,"ERR %s() "
      "got length mismatch (%d/%d)",FCT,SL.len,SR.len);

   Wb::MatProd(HL, X, HX, 'N','N',one,one); 
   Wb::MatProd(X, HR, HX, 'N','N',one,one); 

   for (unsigned i=0; i<SL.len; ++i) {
   Wb::MatProd(SL[i],X, Q,'N','N');
   Wb::MatProd(Q,SR[i],HX,'N','N',one,one); 
   }
};

template <class TQ, class TD>
int genRG_base<TQ,TD>::get1J_gen_aux(const char *F, int L,
   const genRG_base<TQ,TD> &G2, const CData<TQ,TD> &Z,
   wbsparray<TD> &C, double &dx, double &dn, char flag) const {

   unsigned j,n, nblks=0;
   const genRG_base<TQ,TD> &G1(*this);
   double eps=WbUtil<TD>().eps();

   SPIDX_T d=G1.dim(), rsym=G1.Sp.len;
   WBINDEX D; MVEC mark(rsym);
   wbvector< const wbarray<TD>*> AA(rsym);
   wbvector< wbarray<TD> > BB(rsym);
   wbarray<TD> A,B;

   wbMatrix< wbarray<TD> > SM1, SM2;
   wbvector< wbarray<TD> > X;
   wbMatrix<unsigned> Ib,Jb;

   if (!rsym || G2.Sp.len!=rsym || G1.Sz.len!=rsym || G2.Sz.len!=rsym) {
      wblog(FL,"ERR %s() %s [%s] / [%s]\n"
      "got Sp/Sz inconsistency (Sp[%d/%d], Sz[%d/%d])", FCT, STR(q),
      STR(G1.J),STR(G2.J),G1.Sp.len, G2.Sp.len, G1.Sz.len, G2.Sz.len);
   }

   if (G2.dim()!=d) { wblog(F_L,
      "ERR %s() invalid dual irep\n%s: [%s] <> [%s] (d=%d/%d)",
      FCT,STR(q),STR(G1.J),STR(G2.J),d,G2.dim());
   }

   G1.block_decompose_Sm(FL,2,SM1,Ib,&D); nblks=SM1.dim1;

   G2.block_decompose_Sm(FL,1,SM2,Jb);

   if (!nblks || SM2.dim1!=nblks || nblks>d) wblog(FL,
      "ERR %s() got block inconsistency (n=%d/%d / %d)",FCT,nblks,SM2.dim1,d);

   X.init(nblks); 
   RTD nrm=sqrt(RTD(1)/RTD(d));
   X[0].init(1,1,&nrm);

   if (flag & TP3_DBG) wblog(FL,"TST %s() "
      "%d blocks @ d<=%d (D=%d; %s)",FCT, D.len, D.max(), d, TSTR(TD));

   for (unsigned ib=1; ib<nblks; ++ib) {
      const wbarray<TD> *S1=SM1.rec((nblks-1)-ib);
      const wbarray<TD> *S2=SM2.rec(          ib);

      for (n=j=0; j<rsym; ++j) {
         if ((mark[j] =bool(S1[j]))) { ++n; }
         if ( mark[j]!=bool(S2[j]) ) {
            wblog(FL,"ERR %s() block set inconsistency "
            "(%d/%d: %d/%d)",FCT,ib,j,mark[j],bool(S2[j]));
         }

         AA[j]=NULL; BB[j].init();
         if (mark[j]) {
            AA[j]=&S1[j];
            X[Jb(ib,j)].contract(FL,2,S2[j],2,BB[j]);
         }
      }
      if (!n) wblog(FL,"ERR %s() empty block set for ib=%d/%d",FCT,ib,nblks);
      A.Cat(1,AA);
      B.Cat(1,BB); A.mldivide(FL,B,X[ib],-1);
   }

   dx=C.BlockDiag(X,2); 
   dn=fabs(double(C.norm()-TD(1)))/sqrt(double(d));

    { double dS;
      TD x2, x2max=0, one(1);
      wbsparray<TD> X;

      for (unsigned i=0; i<G1.Sp.len; ++i) {
         Wb::MatProd(G1.Sp[i],C, X,'N','N');
         Wb::MatProd(C,G2.Sp[i], X,'N','C',one,one);
         x2=X.norm2(); if (x2max<x2) { x2max=x2; }

         Wb::MatProd(G1.Sp[i],C, X,'C','N');
         Wb::MatProd(C,G2.Sp[i], X,'N','N',one,one);
         x2=X.norm2(); if (x2max<x2) { x2max=x2; }
      }
      dS=sqrt(double(x2max)/d);
      if (dS>sqrt(eps)) wblog(FL,"ERR %s() "
         "error with 1j-symbol @ dS=%.3g = %6.3g eps",FCT,dS,dS/eps); else
      if (flag & TP3_DBG) wblog(FL,
         "TST %s() dS = %.3g = %6.3g eps",FCT,dS,dS/eps);
      dn=MAX(dn,dS);
    }

   if ((flag & TP3_DBG) || dn>=1e-14) {
      wbvec<char> sout(128);
      if (dn<1e-3)
           { sout.catf(FL,"dC=[ sp=%.4g, ck=%.4g",dx/eps,dn/eps); }
      else { sout.catf(FL,"|C|=%.3g !? dx = %.4g",double(C.norm()),dx); eps=1; }

      if (!Z.isEmpty()) {
         wbsparray<TD> Q(Z.cgd); Q-=C; {
            double e=double(Q.norm()); if (e) {
            sout.catf(FL,", st=%.4g",e/eps); }
         }
      }
      sout.catf(FL," ] * %.3g",eps);

      MXPut(FL,"Ix").add(q,"q").add(G1.J,"J").add(D,"D")
         .add(dx,"dx").add(dn,"dn").add(X,"X").add(C,"C").add(Z,"Z");

      if (dn<sqrt(eps)) wblog(FL,
              " *  %s() accuracy [sparse | 1j-check | store]\n%s",FCT,sout.data);
      else wblog(FL,"ERR %s() failed to obtain 1J symbol\n%s",FCT,sout.data);

      return 1;
   }
   return 0;
};

cgdStatus& cgdStatus::init_time(char flag) {

   if ((*this)==CD_UNKNOWN || (!flag && (*this)==CD_ABELIAN)) {
      ctime=mtime=0; cID=0;
      return *this;
   }
   if (!flag && (*this)==CD_IMPLICIT) {
      ctime=mtime=0; cID=CID_RANK1_Q0;
      return *this;
   }

   double now=Wb::getTimeNow();

   if (!flag || flag=='f' || flag=='c') {
      if ((cID || ctime) && flag!='f')  wblog(FL,
         "ERR %s", STR2(*this,'V')); 
      mtime=ctime=now; setID();
   }
   else if (flag=='m') {
      if (!cID || !ctime || ctime>now) wblog(FL,"ERR %s",STR2(*this,'V'));
      mtime=now;
   }
   else wblog(FL,"ERR %s() invalid flag=%s",FCT,cSTR(flag));

   return *this;
};

int cgdStatus::cmp(const char *F, int L, const cgdStatus &b) const {

   if (mtime<ctime || b.mtime<b.ctime) { 
      if (F) wblog(F,L,"ERR %s() invalid mtime",FCT);
      return -12;
   }

   if (isComplete()!=3 || b.isComplete()!=3) {
      if (cID!=b.cID || ctime!=b.ctime) {
         if (F) wblog(F,L,"ERR %s() got CData status mismatch!"
            "%N   %s%N<> %s",FCT, toStr('V').data, b.toStr('V').data);
         return -11;
      }
      return (mtime<b.mtime? -1 : (mtime>b.mtime? +1 : 0));
   }
   else {
      if (cID!=b.cID) {
         if (F) wblog(F,L,"ERR %s() got CData status mismatch!"
            "%N   %s%N<> %s",FCT, toStr('V').data, b.toStr('V').data);
         return -10;
      }
      return 0; 
   }
};

int cgdStatus::sameAs(const cgdStatus &b, char lflag) const {

    int q=0; 

   if (cID!=b.cID || ctime!=b.ctime) {
      if (isComplete()==3 && b.isComplete()==3) {
         if (cID==b.cID) { q=2; } 
         else {
            wblog(FL,"WRN %s() %s with cID = %05X / %05X",FCT,cID,b.cID);
         }
      }
      return q;
   }

   if (lflag>2) {
      if (lflag=='l') lflag=1; else
      if (lflag=='L') lflag=2; else
      wblog(FL,"WRN %s() invalid lflag=%d (using 2='L')",FCT,lflag);
   }

   if (lflag<2 && mtime!=b.mtime) {
      char qa=0, qb=0;
      const cgdStatus &a=*this; int i;

         if ((i=a.isComplete(3))) { if (i==3) { qa|=1; } else { return q; }}
         if ((i=b.isComplete(3))) { if (i==3) { qb|=1; } else { return q; }}

         if (a==CD_STD3) { qa|=2; } else if (a==CD_REF_INIT) { qa|=4; }
         if (b==CD_STD3) { qb|=2; } else if (b==CD_REF_INIT) { qb|=4; }

      if ((qa<4 && qa!=3) || ((qa&2) && a.mtime<b.mtime)) { return q; }
      if ((qb<4 && qb!=3) || ((qb&2) && b.mtime<a.mtime)) { return q; }

      if ((qa|qb)!=7) { return q; }
   }

   if (lflag<1) {
      if (!(ctype & b.ctype & ((1<<CD_NUM_TYPES)-1))) { return q; }
   }

   return (q=1);
};

wbstring cgdStatus::to_vstr(double t) const {

   wbstring s; 

   if (!t) { 
      s="00/00/0000 -- -- --   "; 
   }
   else {
     #pragma omp critical (got_CPTR_TIME)
     { unsigned l=0; time_t tsec=t;
       struct tm *q=localtime(&tsec); s.init(40);

       l=snprintf(s.data,s.len,"%02d/%02d/%d %02d:%02d:%05.2f",
          q->tm_mon+1, q->tm_mday, q->tm_year+1900,
          q->tm_hour, q->tm_min, q->tm_sec+t-tsec  
       );

       if (l>=s.len) wblog(FL,
       "ERR %s() string out of bounds (%d/%d)",FCT,l,s.len);
     }
   }

   return s;
};

wbstring cgdStatus::tstr() const {

   wbvec<char> sout(32); 
   unsigned i=0, n=0, q=1, x=ctype;

   for (; i<CD_NUM_TYPES; ++i, q<<=1) {
      if (ctype&q) { x^=q; ++n; 
         sout.cat(FL,CD_TYPE_STR[i+1]+(sout.l? 1:0));
      }
   }
   if (!n ) { sout=CD_TYPE_STR[0]; }

   for (n=0; i<CD_NUM_TYPES_X; ++i, q<<=1) { 
      if (ctype&q) { x^=q; ++n;
         sout.cat(FL,CD_TYPE_STR[i+1]+(sout.l? (n==1?1:2):0));
      }
   }
   if (sout.l>4) { sout.tr(":?",",,",3); }

   if (x) { wblog(FL,
      "WRN %s() unexpected extraneous %s\n%s (≤%d; %d)",FCT,
      Wb::Bits(x).data, 
      BITS_(x,0),CD_NUM_TYPES,CD_COMPLETE);
   }
   return sout.data;
};

wbstring cgdStatus::toStr(char vflag) const {

   if (!vflag) { return tstr(); } 

   wbvec<char> sout;

   if (mtime<ctime) wblog(FL,
      "ERR %s() got invalid mtime=%.2f (ctime=%.2f)",FCT,mtime,ctime);

   sout.init(vflag? 96:64);
   sout.catf(0,0,"%s (%6s; #%05X", 
      to_vstr(ctime).data, Wb::sec2Str(mtime-ctime).data, cID);

   if (!vflag)
        { sout.cat(0,0,")"); }   
   else { sout.catf(0,0,"; %-5s)%s", tstr().data, u2Str(" u:").data); }

   return sout.data;
};

wbstring cgdStatus::u2Str(const char *istr) const {

   wbstring sout; 

   if (!user[0] && !user[1] && !user[3]) {
      if (user[2]) 
           { sout.init(8); snprintf(sout.data,sout.len," buf=%X",user[2]); }
      else { sout=""; }
   }
   else {
      unsigned i=(istr?strlen(istr):0), m=4;
      wbvec<char> sx(i<16 ? 32 : 16+i);

      sx.catf(FL," %.3s",istr?istr:"");

      for (i=0; i<m; ++i) { if (user[i]>15) break; }
      if (i<m) { 
         for (i=0; i<m; ++i) {
            if (user[i])
                 { sx.catf(FL,"%02X",user[i]); }
            else { sx.cat(FL,"--"); }
         }
      }
      else {
         for (i=0; i<m; ++i) {
            if (user[i])
                 { sx.catf(FL,"%X",user[i]); }
            else { sx.cat(FL,"-"); }
         }
      }
      sout=sx.data;
   }

   return sout;
};

template <class TQ>
QMap<TQ>& QMap<TQ>::init(const char *F, int L,
   const QVec &qv, const wbMatrix<TQ> &q1, const wbMatrix<TQ> &q2,
   const wbMatrix< wbMatrix<TQ> > &QQ,
   const wbMatrix< wbMatrix< const CRef<TQ>* > > *SM
){
   unsigned i,j,k, l=0, n=0;
   unsigned d1=q1.dim1, d2=q2.dim1, d=q1.dim2, mn=d1*d2;
   unsigned nsym=(qv.len ? qv.len : q1.dim2);
   wbindex I; wbperm P; TQ *qq;

   if (!q1.dim1 || !q2.dim1 || !nsym) wblog(F,L,
      "ERR got empty QIDX (%d/%d/%d)",q1.dim1, q2.dim1,nsym);
   if (q2.dim2!=d || QQ.dim1!=d1 || QQ.dim2!=d2) wblog(F,L,
      "ERR severe QDim mismatch (%d/%d, %d/%d, %d/%d)",
       q1.dim2, q2.dim2, QQ.dim1, d1, QQ.dim2, d2);
   if (SM && (SM->dim1!=QQ.dim1 || SM->dim2!=QQ.dim2)) wblog(FL,
      "ERR severe dimension mismatch in SS (%d/%d, %d/%d)",
       SM->dim1, QQ.dim1, SM->dim2, QQ.dim2);

   for (i=0; i<mn; ++i) { const wbMatrix<TQ> &q=QQ.data[i];
      n+=q.dim1; 
      if (q.dim2!=d) wblog(F,L, 
         "ERR inconsistency in qset record length (%d/%d)",q.dim2,d);
      if (SM) {
         const wbMatrix< const wbvector< CRef<TQ> >* > &S=SM->data[i];
         if (S.dim1!=q.dim1 || S.dim2!=nsym) wblog(FL,
            "ERR severe dimension mismatch in SM (%d/%d, %d/%d)",
             S.dim1, q.dim1, S.dim2, nsym
         );
      }
   }

   qvec=qv;
   if (!qvec.len) {
      qvec.init(nsym);
      for (i=0; i<nsym; ++i) qvec[i]=QT_U1;
   }

   Q1=q1; Q2=q2;
   Q.init(n,d); 
   I1.init(n); I2.init(n); II.init(n,2); D.init();

   if (SM)
        { cg3.init(n,nsym); }
   else { cg3.init(); }

   for (l=i=0; i<d1; ++i)
   for (  j=0; j<d2; ++j) {
      const wbMatrix<TQ> &q=QQ(i,j); qq=q.data;
      for (k=0; k<q.dim1; ++k, ++l, qq+=d) {
         I1[l]=i; I2[l]=j;
         Q.recSetP(l,qq);  if (SM) {
         cg3.recSetP(l, (*SM)(i,j).rec(k) ); }
      }
   }

   Q.groupRecs(P,D); I1.Select(P); I2.Select(P);
   if (SM) cg3.recPermute(P);

   for (l=i=0; i<D.len; ++i) { d=D[i];
   for (j=0; j<d; ++j,++l) {
      II(l,0)=i; 
      II(l,1)=j; 
   }}

   return *this;
};

template <class TQ>
void QMap<TQ>::Skip(const wbindex &Ix){

   if (!Ix.len) return;

   unsigned i;
   wbindex I,J;
   wbvector<char> m(D.len), mark(II.dim1);

   if (D.len!=Q.dim1) wblog(FL,
      "ERR %s() size inconsistencty (%d/%d)",FCT,D.len,Q.dim1);
   if (I1.len!=II.dim1 || I2.len!=II.dim1) wblog(FL,
      "ERR %s() size inconsistencty (%d,%d/%d)",FCT,I1.len,I2.len,II.dim1);
   if (!Ix.isUnique(D.len)) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,Ix.max()+1,D.len);

   Ix.invert(D.len,I);
   Q.Set2Recs(I); D.Select(I);

   m.set((const wbvector<unsigned>&)Ix,1);
   mark.set(1);
   for (i=0; i<II.dim1; i++) { if (m[II(i,0)]) mark[i]=0; }
   mark.find(J);

   if (!cg3.isEmpty()) {
      if (cg3.dim1!=II.dim1) wblog(FL,
         "ERR %s() index out of bounds (%d/%d)",FCT,cg3.dim1,II.dim1);
      cg3.Set2Recs(J); 
   }

   I1.Select(J); I2.Select(J); II.Set2Recs(J);
};

template <class TQ>
template <class TD>
int QMap<TQ>::getCGZlist(const char *F, int L,
   wbMatrix<TQ> &qc,       
   wbvector<TD> &dc,       
   wbvector<widx_t>* IC,   
   const wbvector<widx_t>* Ix, 
   const qset<TQ> *q3,     
   wbIndex *m3,            
   const wbvector<TD>* cc,
   double eps1,
   double eps2
) const {

   unsigned i,j,k,l,s,it,nt,ni,mQ,mq, N=0, n1=I1.len, n2=qvec.len;

   wbarray< const wbMatrix<double>* > Z(n1,n2,3); 
   wbMatrix< wbMatrix<SPIDX_T> > IJ(n1,n2);
   wbMatrix< wbvector<TD> > DD(n1,n2); 
   wbvector<unsigned> d0,dz;
   wbvector<char> mark(n1);
   wbindex I;
   double a;

   CDATA_TQ X;
   wbperm p3(3);

   if (I2.len!=n1 || II.dim1!=n1 || n1!=cg3.dim1 || n2!=cg3.dim2) wblog(F,L,
      "ERR %s() severe size mismatch (%d,%d,%d/%d; %d/%d)", 
       FCT, I1.len, I2.len, II.dim1, n1, n2, qvec.len
   );

   mark.set(1);
   if (Ix) {
      for (i=0; i<Ix->len; ++i) { k=Ix->data[i];
         if (k>=mark.len) wblog(FL,
            "ERR %s() index out of bounds (%d/%d)",FCT,k+1,mark.len);
         mark[k]=0;
      }
   }

   if (q3) {
      if (q3->len!=Q.dim2 || II.dim1!=n1) wblog(FL,
         "ERR %s() size mismatch (%d/%d; %d/%d)",
          FCT, q3->len,Q.dim2, II.dim1,n1);
      for (k=0; k<Q.dim1; ++k) { if (!Q.recCompareP(k,q3->data)) break; }
      if (k==Q.dim1) {
          MXPut(FL).add(*this,"M").add(*q3,"q3");
          wblog(FL,"ERR %s()\nno match of Q=[%s] within *this QMap (%d)",
          FCT,q3->toStr().data,Q.dim1);
      }
      for (i=0; i<II.dim1; ++i) { if (II(i,0)!=k) mark[i]=0; }

      if (m3) { 
         if (mark.sum()!=1) wblog(FL,
            "ERR %s() single record expected (%d)",FCT,mark.sum());
         i=mark.find1(); if (int(i)<0) wblog(FL,"ERR %s()",FCT);

         wbvector<unsigned> S(n2);
         for (j=0; j<S.len; ++j) {
            S[j]=cg3(i,j)->cgb->getOM(FL); 
         }

         if (m3->isEmpty()) { m3->init(S); }
         else {
            if (m3->SIZE!=S) wblog(FL,
               "ERR %s() multiplicity setting changed [%s] <> [%s]",
               FCT,STR(m3->SIZE), STR(S));
         }
         if (!(++(*m3))) return 0;
      }
   }
   else if (m3)
   wblog(FL,"ERR %s() m3 may only be used together with q3",FCT);

   if (cc) {
      if (!m3) wblog(FL,
         "ERR %s() cc may only be used together with q3 and m3",FCT);
      if (m3->SIZE.numGT(1)!=1) wblog(FL,
         "ERR %s() cc only applies with outer multiplicity\n"
         "in one symmetry only",FCT);
   }

   if (!mark.any()) {
      wblog(FL,"WRN %s() got empty list into CGR",FCT);
      qc.init(); dc.init(); if (IC) IC->init(N);
      return 0;
   }

   qvec.Qlen(d0,dz);

   for (i=0; i<n1; ++i) { if (mark[i]) { ni=1;
   for (j=0; j<n2; ++j) { const CRef<TQ> &R3=(*cg3(i,j));

       if (!R3) wblog(FL,"ERR %s() got empty cg3(%d,%d)",FCT,i+1,j+1);
       if (!R3.cgb) wblog(FL,"ERR %s() got null cg3(%d,%d).cgb",FCT,i+1,j+1);

       const CDATA_TQ &C3 = (*R3.cgb);

       R3.getP(p3,'i');

       if (p3.len!=3) wblog(FL,
          "ERR %s() got invalid cgp.len=%d",FCT,p3.len);
       if (C3.qs.len!=3*d0[j]) wblog(FL,"ERR CData "
          "qlabel length mismatch (%d / 3*%d)",C3.qs.len,d0[j]);

       if (dz[j]) { for (k=0; k<3; ++k) {
          Z(i,j,k) = &gRS.getR(FL,C3,p3[k]).Z;

          if (Z(i,j,k)->dim2!=dz[j]) wblog(FL,
             "ERR CData length inconsistency (%s) Z(%d,%d,%d): "
             "%s @ dz=%d", C3.qStr().data, i+1,j+1,k+1,
             SSTR_(Z(i,j,k)), dz[j]
          );
       }}

       if (!m3 || m3->SIZE[j]==1) { 
          if (C3.cgd.D.len) { ni*=C3.cgd.nnz();
             C3.cgd.IDX.getCols(p3,IJ(i,j));
             DD(i,j)=C3.cgd.D;
          }
          else if (C3.isAbelian()) { TD x=1;
             IJ(i,j).init(1,3); DD(i,j).init(1,&x);
          }
          else wblog(FL,"ERR %s() got invalid CData\n%s %s",
             FCT,STR(C3), C3.isRefInit() ? "(ref-init!)" : "!?");

          if (R3.getOM()!=1 || !R3.wscalar()) wblog(FL,
             "ERR %s() invalid CRef %s",FCT,SSTR(R3));

          DD(i,j)*=R3.wget0();
       }
       else if (cc) {
          unsigned m=R3.wdim2(); SPARR_RTD x;
          if (cc->len != m3->SIZE[j] || m!=cc->len) wblog(FL,
             "ERR %s() length mismatch of coefficients having OM\n"
             "(%d/%d/%d)",FCT,cc->len,m3->SIZE[j],m);

          for (l=0; l<m; ++l) {
             X.init(R3,l); 
wblog(FL,"ERR %s() check this",FCT);  

             X.cgd*=double(cc->data[l]);
             if (l)
                  { x+=X.cgd; }
             else { x =X.cgd; }
          }
          x.SkipTiny(1e-12);

          ni*=x.nnz();
          x.IDX.getCols(0,2,IJ(i,j));
          DD(i,j)=x.D;
       }
       else { 
          X.init(R3, m3->data[j]);
wblog(FL,"ERR %s() check this",FCT);  

          const SPARR_RTD &x=X.cgd;
          ni*=x.nnz();
          x.IDX.getCols(0,2,IJ(i,j));
          DD(i,j)=x.D;

       }

       if (eps1<=0) {
       k=C3.cgd.D.numZeros(1e-12); if (k) {
          wblog(FL,"WRN cg3(%dx%d): (%d,%d) got zero value%s (%d/%d)",
          n1, n2, i+1,j+1, k!=1 ? "s":"",k, C3.cgd.D.len);
       }}

   }; N+=ni; if (!ni) wblog(FL,"ERR got ni=%d",ni); }}

   mQ=d0.sum()+dz.sum(); l=n2-1;
   qc.init(N,3*mQ); dc.init(N); if (IC) IC->init(N);

   for (k=i=0; i<n1; ++i) { if (!mark[i]) continue;
      wbvec<unsigned> sz_(n2), I_(n2);
      unsigned *sz=sz_.data, *I=I_.data; nt=1;
      for (j=0; j<n2; ++j) { I[j]=0; sz[j]=DD(i,j).len; nt*=sz[j]; }

      for (it=0; it<nt; ++it, ++k) {
         TQ *q=qc.rec(k);
         TD &c=dc[k]; c=1; if (IC) IC->data[k]=i;

         for (j=0; j<n2; ++j) {
            if (!cg3(i,j) || !cg3(i,j) || !cg3(i,j)->cgb)
               wblog(FL,"ERR %s() unexpected cg3 data",FCT);

            const QSet<TQ> &Q=(*cg3(i,j));

            wbMatrix<SPIDX_T> &I1=IJ(i,j);
            mq=Q.t.qlen();

            c*=DD(i,j)[I[j]];

            if (d0[j]) { s=d0[j]*sizeof(TQ);
               memcpy(q     , Q.qs.data     , s);
               memcpy(q+  mQ, Q.qs.data+  mq, s);
               memcpy(q+2*mQ, Q.qs.data+2*mq, s); q+=d0[j];
            }
            if (dz[j]) {
               Wb::cpyRange(q,      Z(i,j,0)->ref(I1(I[j],0)), dz[j]);
               Wb::cpyRange(q+  mQ, Z(i,j,1)->ref(I1(I[j],1)), dz[j]);
               Wb::cpyRange(q+2*mQ, Z(i,j,2)->ref(I1(I[j],2)), dz[j]);
               q+=dz[j];
            }
         }

         a=Wb::abs(c); if (a<eps1) { --k; 
            if (a>eps2) wblog(FL,"WRN %s() got small CGC data "
              "(%.3g; %g %g)",FCT,c,eps1,eps2);
         }

         j=0; I[0]++; 
         while(I[j]>=sz[j] && j<l) { I[j]=0; ++I[++j]; }
      }
   }

   if (k!=N) {
      if (!k || k>N) wblog(FL,
         "ERR %s() got no CG data (%d/%d)",FCT,k,N);
      qc.dim1=k; dc.len=k; if (IC) IC->len=k;
   }

   return 1;
};

template <class TQ>
template <class TD>
void QMap<TQ>::getIdentityQ(const char *F, int L,
   const wbvector<widx_t> &Sa, 
   const wbvector<widx_t> &Sb, 
   QSpace<TQ,TD> &A,     
   char vflag
) const {

   unsigned i,j,k,l,M,i1,i2,i3,i4,S3;
   unsigned n=D.sum(), nq=qvec.len, QDIM=qvec.Qlen();
   char isa=qvec.allAbelian(), cgflag=0;

   double cfac=1, DS=0; 

   wbMatrix<widx_t> mm(n,nq); 
   wbvector<unsigned> d0,dz;
   wbvector< wbvector<widx_t> > DD(D.len), MM(D.len);
   wbvector<widx_t> S;

   qvec.Qlen(d0,dz);

   if (Q1.dim1!=Sa.len || Q1.dim2!=QDIM || Q.dim1!=D.len ||
       Q2.dim1!=Sb.len || Q2.dim2!=QDIM || Q.dim2!=QDIM) wblog(FL,
      "ERR QMap::Id block size mismatch\n"
      "%s / %dx%d ; %s / %dx%d ; %s / %dx%d", SSTR(Q1), Sa.len,QDIM,
       SSTR(Q2), Sb.len,QDIM, SSTR(Q), D.len,QDIM);

   if (I1.len!=n || I2.len!=n || II.dim1!=n || cg3.dim1!=n || cg3.dim2!=nq)
       wblog(F,L,"ERR QMap::Id size mismatch (%d,%d,%d / %s / %d)",
       I1.len, I2.len, II.dim1, SSTR(cg3), nq);

   for (i=0; i<D.len; ++i) { DD[i].init(D[i]); MM[i].init(D[i]); }
   for (i=0; i<n; ++i) { i1=I1[i]; i2=I2[i]; i3=II(i,0); i4=II(i,1);
      widx_t &Di=DD.at(i3).at(i4);
      if (Di) wblog(FL,"ERR QMap::Id " 
         "block index not unique (%d,%d,%d: %ld)",i+1,i3+1,i4+1,Di);
      Di = Sa.at(i1) * Sb.at(i2); 
   }

   for (i=0; i<n; ++i) {
      for (M=1, j=0; j<nq; ++j) { const CRef<TQ> *R3=cg3(i,j);
         if (!R3) wblog(FL,"ERR %s() got null cg3 data",FCT);
         if (!R3[0]) wblog(FL,"ERR %s() got empty cg3 data",FCT,STR_(R3));
         if (!R3->cgb) wblog(FL,"ERR %s() got null cg3.cgb",FCT);
         if (R3->cgb->qs.len!=3*d0[j]) wblog(FL,"ERR QSet length mismatch "
            "(%d,%d): %d / 3*%d)",i+1, j+1, R3->cgb->qs.len, d0[j]);
         M*=( mm(i,j) = R3->wdim() );
      }

      MM[II(i,0)][II(i,1)]=M; 

   }

   A.init(n,3,QDIM); 

   if (!isa ) { A.qtype=qvec; A.setupCGR(); cgflag=1; } else
   if (isa>1) { A.qtype=qvec; }

   for (j=0; j<cg3.dim2; ++j) { if (qvec[j].isAbelian()) {
   for (i=0; i<cg3.dim1; ++i) { l=0;
       if (cg3(i,j)) { if (!cg3(i,j)->wscalar1()) { l|=1; }}
       else if (cgflag) { l|=2; }
       if (l) wblog(FL,"ERR %s() invalid abelian cg3 data "
       "(%d,%d: %s, e=%d)",FCT,i+1,j+1, l&1 ? STR_(cg3(i,j)):"null", l);
   }}}

   for (i=0; i<n; ++i) {
      i1=I1[i]; i2=I2[i]; i3=II(i,0); i4=II(i,1);

      { size_t *d3=DD[i3].data, *m3=MM[i3].data, n4=D[i3];
        for (k=l=0; l<i4; ++l) { k  += d3[l]*m3[l]; } 
        for ( S3=k; l<n4; ++l) { S3 += d3[l]*m3[l]; }
      }

      A.QIDX.recSetB(i,0,QDIM, Q1.rec(i1));
      A.QIDX.recSetB(i,1,QDIM, Q2.rec(i2));
      A.QIDX.recSetB(i,2,QDIM, Q .rec(i3));

      if (cgflag) { M=MM[i3][i4];
         for (cfac=1, j=0; j<nq; ++j) {
            cfac *= A.CGR(i,j).init(*cg3(i,j)).NormSignW();
         }

         A.DATA[i]->initIdentityB3( Sa[i1], Sb[i2]*M, S3, k, TD(cfac));
         if (M>1) {
            size_t ss[]={ Sa[i1], Sb[i2], M, S3 };     
            A.DATA[i]->Reshape(S.init(4,ss)).Permute("1243");
         }
      }
      else {
         A.DATA[i]->initIdentityB3( Sa[i1], Sb[i2], S3, k, TD(1));
      }
      DS+=A.DATA[i]->numel();
   }

   A.otype=QS_AMATRIX;
   A.itags.init3();

   wbvector<widx_t> s;
   A.getDim(s,&S);

   if (S[2]!=S[0]*S[1]) { 
      MXPut(FL).add(*this,"M").add(A,"A").add(s,"sd").add(S,"S");
      wblog(FL,"ERR %s() overall dimension mismatch\n[ %s ]: [ %s ]",
      FCT,STR(s),STR(S));
   }

   if (vflag) {
      double Mtot=wbsys::getMemTot();
      char q=(DS>Mtot? 1 : 0); if (DS>2E8 && vflag=='V') { q|=2; }
      if (q) { 
         DS  /=double(1<<27); 
         Mtot/=double(1<<27);
         if (q&1)
              wblog(FL,"WRN %s() cID using %.3g / %.3g GB",myname,DS,Mtot);
         else wblog(FL,"WRN %s() cID using %.3g GB (%d blocks)",myname,DS,n);
      }
   }
};

template <class TQ>
int CStore<TQ>::checkInit(const char *F, int L, const QType &t) const {
    genRG_struct<TQ,RTD> &B=gRS.buf[t];
    if (B.q.isKnown() || t.isAbelian()) { return 0; }
    else { B.checkInit(F_L,t); return 1; }
};

template <class TQ>
CDATA_TQ* CStore<TQ>::BUF_find(const QSet<TQ> &Q) {

#ifdef QS_USING_OMP
   Wb::ompGuard gLK(CS_buf); 
#endif

   auto it = BUF.find(Q);
   return (it!=BUF.end() ? &(it->second) : nullptr);
};

template <class TQ>
int CStore<TQ>::Reduce2Ref(const char *F, int L, const QSet<TQ> &Q,
   char force, 
   unsigned Dmin 
){
   int rval=0; 

#ifdef QS_USING_OMP
   CG::Guard qLK(FL,Q,"buf");
#endif

   auto it=BUF.find(Q);

   if (force>3) {
      if (force=='!' || force=='f') { force=3; } else
      if (force=='l' || force=='t') { force=0; } 
      else wblog(FL,"ERR %s() invalid force=%d<%c>",FCT,force,force);
   }

   if (it==BUF.end()) {
      if ((force&1) && F) wblog(F,L,
         "ERR %s() missing QSet %s",FCT,STR(Q));
      rval=-99;
   }
   else if (it->second.cstat!=CD_REF_INIT) {
      if (force>1 || it->second.sizeable(Dmin)) {
         rval=it->second.Reduce2Ref(F,L, force>1 ? '!':0);
      }
      else { rval=-2; } 
   }

   return rval;
};

template <class TQ>
const CDATA_TQ& CStore<TQ>::getIdentityC(
   const char *F, int L,
   const QType &t, const TQ *qs, unsigned dim, unsigned loadRC
){
   unsigned dim_;

   if (t.isAbelian()) { dim_=1;
      if (int(dim)<0) { dim=1; } else 
      if (dim!=1) wblog(FL,"ERR %s() abelian CGC with dim=%d",FCT,dim);
   }
   else { dim_=t.qdim(qs); }

   QSet<TQ> Q; { Q.init2(t,qs); }
   CDATA_TQ &C = getBUF(0,0,Q,loadRC);

   if (int(dim)<0) { dim=dim_; } else
   if (dim!=dim_) wblog(F_L,
      "ERR %s() got dim=%d/%d having %s",FCT,dim,dim_,STR(Q));

   if (C.cgd.isEmpty() && C.QSet<TQ>::isEmpty()) {
      RTD nrm=Wb::sqrt(RTD(1)/RTD(dim)); 

      C.t=t; C.qs=Q.qs; C.qdir.init_iout(FL,2,2);
      C.cgd.initIdentity(t,dim,nrm);
      C.cstat_init_(CD_IDENTITY).setComplete(3);

      gStore.rclog(t, PFL, CG_VERBOSE>6 && F,
         "[+] CBUF[%03d|Id] #%05X %s",BUF.size(),C.cstat.cID,STR2(C,'v'));
      if (C!=Q) wblog(FL,
         "ERR %s() QSet inconsisteny\n%s",FCT,STR(Q));
      gStore.save_CData(FL,C); 
   }
   else {
      if (C.t!=t || C.qs!=Q.qs || C.qdir!=Q.qdir) wblog(F_L,
         "ERR %s() CGC Id symmetry mismatch (%s <> %s)",
         FCT,STR(C),STR(Q));
      if (!C.isAbelian() && !C.cgd.isSMatrix(F_L,dim)) wblog(F_L,
         "ERR %s() CGC Id size mismatch (%s; %d)",FCT,SSTR(C),dim);
      if (C.qdir!="+-") wblog(F_L,
         "ERR %s() invalid C.qdir=[%s]",FCT,STR(C.qdir)
      );
   }

   return C;
};

template <class TQ>
const CDATA_TQ& CStore<TQ>::getIdentity1J(
   const char *F, int L, CRef<TQ> &C,
   const QType &t, const TQ *qs, unsigned dim, unsigned loadRC
){
   char isa=t.isAbelian();
   if (isa) {
      if (dim!=1) wblog(FL,
         "ERR %s() abelian CGC with dim=%d",FCT,dim);
   }

   QSet<TQ> Q;
   Q.init1J(t,qs); Q.Sort(&C.cgp,'i');
   C.cgb = &getBUF(0,0,Q,loadRC);

   if (C.cgb->QSet<TQ>::isEmpty()) wblog(FL,
      "ERR %s() CData %s not yet defined\n(l=%d '%s')",
      FCT, STR(Q), loadRC, BITS(loadRC));

   const CDATA_TQ &c = *(C.cgb);

   if (int(dim)<0) dim=c.cgd.dim(); 
   if (dim!=t.qdim(qs)) wblog(F_L,
      "ERR %s() got dim=%d having %s",FCT,dim,STR(Q));

   C.cgw.init(1,1); C.cgw[0]=sqrt(double(dim));

#ifndef WB_SKIP_ASSERT
   if (c.t!=t || c.qs!=Q.qs || c.qdir!=Q.qdir) wblog(F_L,
      "ERR %s() CData symmetry mismatch (%s <> %s)",FCT,STR(c),STR(Q));
   if (!c.t.isAbelian()) { if (!c.cgd.isSMatrix(F_L,dim)) wblog(F_L,
      "ERR %s() CGC Id size mismatch (%s; %d)",FCT,SSTR(c),dim); }
   else if (dim!=1 || !c.cgd.isEmpty()) wblog(FL,
      "ERR %s() got abelian cgdata (d=%d)",FCT,dim);
   if (c.qdir!="++") wblog(F_L,
      "ERR %s() invalid c.qdir=[%s]",FCT,STR(c.qdir));
#endif

   return c;
};

template <class TQ>
int CStore<TQ>::valid_mp3_data(const QType &t,
   const qset<TQ> &J1, const qset<TQ> &J2, const qset<TQ> &J3,
   char rflag
) const {

   auto iJ12 = gCS.get_mp3_data(FL,t,J1,J2);
   const auto &M2 = iJ12->second;

   for (auto I2=M2.begin(); I2!=M2.end(); ++I2) {
      if (rflag)           
           { if (J3.isEqualR(I2->first)) return 1; }
      else { if (J3.isEqual (I2->first)) return 1; }
   }

   return 0;
};

template <class TQ>
void CStore<TQ>::getQfinal_1(
   const char *F, int L, const QType &t,
   const qset<TQ> &j1, const qset<TQ> &j2,
   wbMatrix<TQ> &jj, 
   wbMatrix< const CRef<TQ>* > *sm,
   char loadC 
){

   auto iJ12 = get_mp3_data(F,L,t,j1,j2);
   auto &M2 = iJ12->second;
   unsigned i=0, m, n=M2.size(), d=j1.len;

   if (d!=t.qlen()) wblog(FL, 
      "ERR invalid qset record length (%d/%d)",d,t.qlen());
   if (!n) wblog(FL,"ERR %s() got empty map3 data",FCT);

   jj.init(n,d); 

   if (sm) { sm->init(n,1); }

   for (auto I2=M2.begin(); I2!=M2.end(); ++I2, ++i) {
      const qset<TQ> &J=I2->first;

      if (J.len!=d) wblog(FL, 
         "ERR qset length mismatch (%d/%d)",J.len,d);

      jj.recSetP(i,J.data);
      if (sm) {
         CRef<TQ> &R3=(I2->second); 
         sm->data[i]=&R3; 

         if (!R3) wblog(FL,"ERR %s() got empty CG3s (%s)",FCT,STR(t));

         if (!R3.cgb) wblog(FL,
            "ERR %s() got empty CRef (%s; %s)",FCT,STR(t),STR(R3));

         if (loadC && R3.isRefInit()) { R3.LoadRef(FL); }

         const wbvector<unsigned> &S=R3.cgb->cgd.SIZE;
         if (S.len) { 
            if (S.len<3 || S.len>4) wblog(FL,
            "ERR %s() unexpected size %s for %s",FCT,SSTR(R3),STR(R3));
         }

         if ((m=R3.cgb->getOM())>1) {
            if (!R3.wdim_is(m) || R3.cgb->cgd.SIZE.len!=4 || m!=R3.Size(3))
            wblog(FL,"ERR %s() OM size mismatch (%s /%d)",FCT,SSTR(R3.cgw),m);
         }
      }
   }
};

template <class TQ>
void CStore<TQ>::getQfinal_v(
   const char *F, int L, const QVec &qvec,
   const qset<TQ> &q1, const qset<TQ> &q2, wbMatrix<TQ> &QQ,
   wbMatrix< const CRef<TQ>* > *SM, 
   char loadC 
){
   unsigned i,l,d=0,D;
   wbMatrix<TQ> X;
   qset<TQ> J1,J2;

   wbMatrix< const CRef<TQ>* > sm;
   wbvector<unsigned> dd;

   if (qvec.len) { D=qvec.Qlen(dd); }
   else {
      dd.init2val(q1.len,1); 
      D=q1.len; 
   }

   if (D!=q1.len || D!=q2.len) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d/%d)",q1.len,q2.len,D);

   QQ.init(); if (SM) { SM->init(); }

   for (l=i=0; i<dd.len; ++i, l+=d) { d=dd[i];
      J1.init(d,q1.data+l);
      J2.init(d,q2.data+l);

      getQfinal_1(F,L,
          qvec.len ? qvec[i] : QT_U1, 
          J1, J2, X, SM? &sm:nullptr, loadC
      );
      QQ.ColKron(X);   

      if (SM) { SM->ColKron(sm); } 
   }
};

template <class TQ>
void CStore<TQ>::getQfinal( 
  const char *F, int L, const QVec &qvec,
  const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, QMap<TQ> &M,
  unsigned char cgflag
){
  unsigned i,j;
  unsigned d1=Q1.dim1, d2=Q2.dim1, d=(qvec.len ? qvec.Qlen() : Q1.dim2);

  wbvector<unsigned> D;
  wbMatrix< wbMatrix<TQ> > QQ;
  wbMatrix< wbMatrix< const CRef<TQ>* > > SM;
  wbperm P;

  qset<TQ> q1,q2;

  if (Q1.isEmpty() || Q2.isEmpty()) wblog(F,L,
     "ERR got empty input spaces (%d,%d)",Q1.isEmpty(),Q2.isEmpty());
  if (!Q1.isUnique() || !Q2.isUnique()) wblog(F,L,
      "ERR expect input to be sorted and unique");

   if (!d1 || !d2) {
      wblog(F,L,"WRN %s() got empty Q-Space (%d/%d)",d1,d2);
      M.init(); return;
   }
   if (d!=Q1.dim2 || d!=Q2.dim2) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d/%d; %s)",
       Q1.dim2, Q2.dim2, d, STR(qvec));

   QQ.init(d1,d2); if (cgflag) {
   SM.init(d1,d2); }

   for (i=0; i<d1; ++i)
   for (j=0; j<d2; ++j) {
      q1.init2ref(d,Q1.rec(i));
      q2.init2ref(d,Q2.rec(j));
      getQfinal_v(F,L,qvec,q1,q2, QQ(i,j),
         cgflag ? &SM(i,j) : nullptr, 
         cgflag & QF_LOADC 
      );
   }

   M.init(FL,qvec,Q1,Q2,QQ, cgflag? &SM:nullptr);
};

template <class TQ>
int CStore<TQ>::getQfinal_zdim(
   const char *F, int L, const QVec &qvec,
   const wbMatrix<TQ> &Q1, const wbMatrix<TQ> &Q2, const wbMatrix<TQ> &Q,
   wbvector<widx_t> &s1, wbvector<widx_t> &s2, wbvector<widx_t> &s3,
   wbvector<unsigned> *M, wbMatrix< const CRef<TQ>* > *S3
){
   unsigned i,j,k,l;
   unsigned d1=Q1.dim1, d2=Q2.dim1, d3=Q.dim1,
       d=(qvec.len ? qvec.Qlen() : Q2.dim2),
       nsym=(qvec.len ? qvec.len : Q2.dim2); 

   wbvector<widx_t> D;
   wbvector<char> m3;
   wbindex i0,i2;
   wbMatrix< wbMatrix<TQ> > QQ;
   wbMatrix< wbMatrix< const CRef<TQ>* > > SS;
   wbperm P;

   qset<TQ> q1,q2;

   if (!d1 || !d2 || !nsym) {
      if (F) wblog(F,L,
         "WRN %s() got empty Q-Space (%d/%d/%d)",d1,d2,nsym);
      s1.init(); s2.init(); s3.init();
         if (M ) M ->init();
         if (S3) S3->init();
      return -1;
   }
   if (!Q1.isUnique() || !Q2.isUnique() || !Q.isUnique()) wblog(F,L,
      "ERR expect input to be sorted and unique [%d,%d,%d]",
       Q1.isUnique(), Q2.isUnique(), Q.isUnique());
   if (d!=Q1.dim2 || d!=Q2.dim2 || d!=Q.dim2) wblog(F,L,
      "ERR inconsistent QSpaces (%d,%d,%d/%d; %s)",
       Q1.dim2, Q2.dim2, Q.dim2, d, STR(qvec));

   QQ.init(d1,d2);
   SS.init(d1,d2);

   for (i=0; i<d1; ++i)
   for (j=0; j<d2; ++j) {
      q1.init2ref(d,Q1.rec(i));
      q2.init2ref(d,Q2.rec(j));
      getQfinal_v(F,L,qvec,q1,q2, QQ(i,j), &SS(i,j));
   }

   s1.init(d1).set(1);
   for (i=0; i<d1; ++i) {
       const wbMatrix< const CRef<TQ>* > &Si=SS(i,0);
       for (l=0; l<Si.dim2; ++l) { s1[i]*=Si(0,l)->Size(0); }
   }

   s2.init(d2).set(1);
   for (j=0; j<d2; ++j) {
       const wbMatrix< const CRef<TQ>* > &Sj=SS(0,j);
       for (l=0; l<Sj.dim2; ++l) { s2[j]*=Sj(0,l)->Size(1); }
   }

   if (Q.isEmpty()) { return 0; }

   s3.init(d3).set(1);
      if ( M) M ->init(d3).set(1);
      if (S3) S3->init(d3,SS(0,0).dim2); 
   m3.init(d3);

   for (i=0; i<d1; ++i)
   for (j=0; j<d2; ++j) {
       const wbMatrix< const CRef<TQ>* > &Sij=SS(i,j);
       matchIndex(QQ(i,j),Q,i0,i2);
       if (!i2.len) continue; 

       for (k=0; k<i0.len; ++k) { if (++m3[i2[k]]>1) continue;
          unsigned m_, id=i2[k];
          if (S3) { S3->recSetP(id,Sij.ref(i0[k])); }

          if (Sij.dim2!=nsym) wblog(FL,"ERR %s() inconsistent CG3 "
             "(%d; %s)",FCT,Sij.dim2,STR(qvec));

          for (l=0; l<Sij.dim2; ++l) {
             if (qvec.len && !qvec[l].isAbelian()) {
                const CRef<TQ> &r3=(*Sij(i0[k],l));
                s3[id] *= r3.Size(2);
                if ((m_=r3.getOM(FL))>1) {
                   if (M)
                        { M->data[id] *= m_; }
                   else { s3[id] *= m_; }
                }
             }
          }
       }
       if (!m3.contains(0)) break;
   }

   if (m3.contains(0)) {
       matchIndex(Q1,Q,i0,i2);
       for (k=0; k<i0.len; ++k) { unsigned id=i2[k];
           if (m3[id]) continue; else m3[id]+=100;
           s3[id]=s1[i0[k]];
       }

       if (m3.contains(0)) {
       matchIndex(Q2,Q,i0,i2);
       for (k=0; k<i0.len; ++k) { unsigned id=i2[k];
           if (m3[id]) continue; else m3[id]+=101;
           s3[id]=s2[i0[k]];
       }}
   }

   if (int(i=m3.find1(0))>=0) {
      MXPut X(FL,"a");
      X.add(Q1,"Q1").add(Q2,"Q2").add(Q,"Q").add(QQ,"QQ")
       .add(s1,"s1").add(s2,"s2").add(s3,"s3").add(m3,"m3").add(i,"i");
      if (M) { X.add(*M,"M"); }; X.put("caller");

      if (F) wblog(F_L,
         "ERR %s() invalid output symmetry [%s]",FCT,Q.rec2Str(i).data);
      else return -1;
   }

   return 0;
};

template <class TQ> inline
void CStore<TQ>::add_CData(
   const QType &q, const qset<TQ> &J1,const qset<TQ> &J2
){
   switch (q.type) {

     case QT_U1  :
     case QT_ZN  :
     case QT_P   :

        if (J1.len!=1 || J2.len!=1) wblog(FL,
           "ERR %s() got len=%d/%d for %s",FCT,J1.len,J2.len,STR(q));
        add_CData_abelian(q,J1[0],J2[0]);
        break;

     case QT_SUN :
     case QT_SpN :
     case QT_SON :
     case QT_SEN :
     case QT_A4  : {
        int i=gStore.load_mp3(FL,q,J1,J2);
        if (i<0) wblog(FL,"ERR %s() failed to load/generate\n"
           "CGC data for %s [%s; %s] (i=%d)",
           FCT,STR(q), STR(J1), STR(J2),i);
        }
        break;

     default:
        wblog(FL,"ERR symmetry %s not implemented yet (%d)",
        STR(q),q.type);
   }
};

template <class TQ>
void CStore<TQ>::add_CData_abelian(
   const QType &q, const TQ &J1, const TQ &J2
){
   unsigned i=0, m=(J1!=J2 ? 2 : 1);
   qset<TQ> J12(2), J21(2), J(1);

   J12[0]=J1; J12[1]=J2;
   J[0]=getTensorProdReps_abelian(q,J1,J2);

   if (J12[0]<=J12[1])
        { J21[0]=J12[1]; J21[1]=J12[0]; }
   else { J21[0]=J12[0]; J21[1]=J12[1]; J12[0]=J21[1]; J12[1]=J21[0]; }

   gStore.rclog(q, PFL, CG_VERBOSE>6,
      "[+] map3[%02d] %-5s %s [%2g %2g; %2g ]",
      map3[q].size(), STR(q), m>1 ? "X":"=",
      double(J12[0]), double(J12[1]), double(J[0])
   );

   CDATA_TQ S;
   S.init3_abelian(q,J12[0],J12[1],J[0]);

   CDATA_TQ &Sb=getBUF(0,0,(QSet<TQ>&)S, 0);
   S.save2(Sb);

   for (i=0; i<m; ++i) {
      CRef<TQ> &R=map3[q][i==0 ? J12 : J21][J];

      if (R.stat() ^ (R.cgb!=nullptr)) wblog(FL,
         "ERR %s() CRef status mismatch: %s", FCT, R.statStr().data);

      if (R.cgb) {
         if (!R.cgb->sameAs(Sb,'l')) wblog(FL,
            "ERR %s() inconsistent %s CGC\n   %s\n<> %s",
            FCT,STR(q),R.cgb->toStr().data, STR(Sb));
         if (R.cgw.numel()!=1 || R.cgw[0]!=1 || R.cgp.conj ||
            (!i &&  R.cgp.len) ||
            ( i && (R.cgp.len!=3 || R.cgp[0]!=1 || R.cgp[1]!=0))
          ) wblog(FL,"ERR %s() invalid scalar CRef (%d)\n%s",
            FCT, i+1, STR(R));
      }
      else {
         if (R.cgw || R.cgp) wblog(FL,
            "ERR %s() got unexpected empty CRef\n%s",FCT,STR(R));
         R.cgb=&Sb;
         R.cgw.init(1,1); R.cgw[0]=1; if (i) {
         R.cgp.initStr(FL,"213"); }
      }
   }
};

template <class TQ>
void CStore<TQ>::add2BUF(
   const char *F, int L, const CDATA_TQ &S,
   char nflag 
){
   unsigned r=S.rank(F_L);
   if (S.t.validType()<=0 || r!=3)
      wblog(FL,"ERR %s() got invalid rank-%d CData [%s]",
      FCT, r, S.qStr().data);
   if (S.cgd.SIZE.len<3 || S.cgd.SIZE.len>4)
      wblog(FL,"ERR %s() got invalid CData [%s: %s @ r=%d]",
      FCT, S.qStr().data, SSTR(S.cgd), r
   );

   CDATA_TQ &B = getBUF(FL,QSet<TQ>(S));
   char gotB=(B.isEmpty()? 0 : B.isRefInit()? -1 : 1);

   if (gotB<=0) {
      if (gotB) {
         if (S.cgd.SIZE!=B.cgd.SIZE || (QSet<TQ>&)S!=(QSet<TQ>&)B)
         wblog(FL,"ERR %s() inconsistent size-ref data\n"
         "%s <> %s",FCT,STR(S), STR(B));
      }
      gStore.rclog(S.t, PFL, CG_VERBOSE>6 && F,
         "[+] CBUF[%02d] %s %s", BUF.size(), STR(S),
         gotB==0 ? "": " [replaces Sref]"
      );
      B=S;

      gStore.save_CData(FL,S); 
   }
   else {
      char s[32]; s[0]=0;
      if (nflag)
           strcat(s,"already got existing entry"); else
      if (B.QSet<TQ>::operator!=(S))
           strcat(s,"inconsistent pre-existing data");
      else strcat(s,"[TST double check this]"); 

      if (s[0]) {
         MXPut(FL,"q").add(S,"A").add(B,"A_").add(nflag,"nflag")
           .add(wbstring(s),"str");
         wblog(FL,"ERR %s() %s\nB: %s\nS: %s",FCT,s,STR(B),STR(S));
      }
   }
};

template <class TQ>
CDATA_TQ& CStore<TQ>::getBUF(
   const char *F, int L, const QSet<TQ> &Q, unsigned loadRC) {

   int q, isnew=0; unsigned r=Q.rank(FL);
   CDATA_TQ *Cb=NULL; 

#ifndef WB_SKIP_ASSERT
   if (Q.isEmpty()) wblog(F_L,"ERR %s() got empty QSet",FCT);
   if (!Q.isSorted()) wblog(F_L,
      "ERR %s() CStore requires sorted QSet\n%s",FCT,STR(Q));
#endif

   {
#ifdef QS_USING_OMP
      Wb::ompGuard gLK(CS_buf);
#endif
      auto it=BUF.find(Q);
      if (it!=BUF.end()) { Cb=&(it->second);
         if (Cb->isEmpty(0)) {
            isnew|=2; 
         }
      }
      else { Cb=&(BUF[Q]);
         isnew|=1; 
      }

      if (Q.t.isAbelian()) { 
         if (isnew || Cb->isEmpty()) {
            Cb->initAbelian(Q); 
            Cb->setuser_BUF();
         }
         return *Cb;
      }
   }

#ifdef QS_USING_OMP

   CG::Guard qLK(FL,Q,"buf", 
      isnew || (loadRC && !(loadRC & LB_REF)) ? 0:"?"); 

   if (isnew && !Cb->isEmpty()) { isnew=0; }
#endif

   if (isnew) { Cb->setuser_BUF();
      if (!(loadRC & LB_ANY) && !F) { return *Cb; } 
   }
   else { 
      if (!(loadRC & LB_UPD) 
      && (!(loadRC & LB_GEN) || Cb->cstat!=CD_REF_INIT)) { return *Cb; }
   }

   Cb->setuser_BUF_active(); 

#if defined(DBG_QSX_BUF) && (DBG_QSX_BUF & 1)
   { wbvec<char> sx(128);
     sx.catf(FL,"%s ", isnew ? " + " :
        ((loadRC && !(loadRC & LB_REF)) ? " ->":"  r"));
     sx.catf(FL,"%-10s %3ld %p  %-20s%s",
        FCT, BUF.size(), Cb, STR(Q), Cb->cstat.u2Str(" ").data);

     if (!isnew && l<n) {
        sx.catf(FL," %p", Cb->cgd.SIZE.data);
     }
     BFF.blogf(PFL,sx.data);

     #if 0 && defined(DBSTOP)
        if (!isnew) { QSet<TQ> X;
           { TQ q[4]={0,1,0,1}; if (*Cb==X.init(Cb->t,q,"++--")) dbstop(FL); }
           { TQ q[3]={1,3,4  }; if (*Cb==X.init(Cb->t,q,"++-" )) dbstop(FL); }
        }
     #endif
   }
#endif

#if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 1 )
   if (Wb::omp_parallel()) { 
      LKF.blogf(FL,"%s() %-24s %p #%d",FCT,STR(Q),Cb,isnew);
   }
#endif

   unsigned mC=(isnew? 0 : Cb->numOM()); 

   if (isnew && Q.isZero()>0) {
      Cb->cstat.setComplete(5);    
      Cb->setuser_BUF_passive();
      return *Cb;
   }
   else if (isnew && (loadRC & LB_REF)) {
      q=gStore.load_CData(0,0,Q,*Cb,'r'); 
      if (q<=0 && F) {
         if (q==0) wblog(FL,"ERR %s() no CStore file for %s\n"
         "hint: check RC_STORE setting",FCT,STR(Q));
         else wblog(FL,"ERR %s() "
         "failed to read CRef from file (q=%d)\n%s",FCT,q,STR(Q));
      }
   }
   else if (isnew || Cb->cstat==CD_REF_INIT) { if (loadRC & LB_GEN) {
      q=gStore.load_CData(0,0,Q,*Cb); 

      if (q<=0 && r>2 && (loadRC & LB_CALC) && (r<=3 || Wb::envFullOM)) {
         if (CG::thread_fOM.me()) {
            loadRC &= ~LB_CALC; 
         }
         else if (Cb->gotuser_BUF_active()<=1) {
            if ((q=Cb->hasFullOM(FL))!=-1 || !isnew) wblog(FL,
               "WRN %s() got fOM=%d, isnew=%d",FCT,q,isnew);

            Cb->setuser_BUF_active(3);  
            q=Q.generateFullOM(FL);     
            Cb->setuser_BUF_passive();

            if (q<0) { 
               Cb->cstat.setComplete(5); 
               return *Cb;
            }
            else {
               q=gStore.load_CData(0,0,Q,*Cb); 
               if (q<=0) wblog(FL,"ERR %s() %s => q=%d",FCT,STR(Q),q);
            }
         }
      } 

      if (q>0) { if (CG_VERBOSE>8) { char tag[6]="";
         if (q==1) {
            if (CG_VERBOSE>4) { strcpy(tag,"(+)"); } 
         }
         else if (q==3) {
            if ((Cb->qdir.len>2 && CG_VERBOSE>1) || CG_VERBOSE>4) {
               strcpy(tag,"{+}");
            }
         }
         else if (q==4) {
            if ((Cb->qdir.len>2 && CG_VERBOSE>1) || CG_VERBOSE>4) {
               strcpy(tag,"{-}");
            }
         }
         else if (q==5) {
            if (CG_VERBOSE>4) { strcpy(tag,"(d)"); }
         }
         else if (q==6) {
            if (CG_VERBOSE>5) { strcpy(tag," ok"); }
         }
         else wblog(FL,"ERR %s() load_CData() returned q=%d",FCT,q);

         if (tag[0]) { 
            wblog(PFL,"%3s CBUF[%03ld] %s %s",tag,BUF.size(),STR_(Cb),
               CG_VERBOSE>8 && !Cb->isRefInit() ? 
               Wb::size2Str(Cb->memSize()).data : "");
            if (CG_VERBOSE>8) { fflush(0); }
         }
      }}
      else if (r==2) {
         unsigned flag=(loadRC & TP3_ITER);

         if (Q.isScalar()) { 
            qset<TQ> J(Q.t.qlen(),Q.qs.data,'r');
            if (gStore.load_RSet(0,1,Q.t,J)>0) {
               gCS.getIdentityC(FL,Q.t,Q.qs.data,-1,flag);
            }
         }
         else if (Q.is1J()) {
            if (Q.t.qrank()) {
               unsigned m=Q.t.qlen();
               qset<TQ> q(m,Q.qs.data,'r'), q_(m,Q.qs.data+m,'r');
               if (gStore.load_RSet(0,1,Q.t,q )>0 &&
                   gStore.load_RSet(0,1,Q.t,q_)>0
                  ){
                   genRG_struct<TQ,RTD> &B=gRS.Buf(Q.t);
                   if (B.q==QT_UNKNOWN) { B.SetupSym(FL,Q.t); }
                   B.get1J_gen(F_L,q,flag); 
               }
               else wblog(FL,
               "ERR %s() missing dual irep for 1J symbol\n%s",FCT,STR(Q));
            }
            else if (loadRC & LB_CALC) { 
               wblog(FL,"ERR %s() missing 1J %s\n%s",FCT,STR(Q),STR2_(Cb,3));
            }
         }
         else if (loadRC & LB_CALC__) { 
            if (*Cb || Cb->cgd || Cb->cstat.isComplete(7)
              || Q.isZero()<=0) { wblog(FL,"WRN %s() %s",FCT,STR2_(Cb,7)); }
            Cb->cstat.setComplete(5); 
         }
         else { wblog(FL,
            "ERR %s() got rank-%d QSet (loadRC=%s)\n%s",
            FCT, r, BITS(loadRC), STR(Q));
         }
      }
      else if (r<2) { 
         if (r!=1) wblog(F_L,"ERR %s() got QSet %s",FCT,STR(Q));
         Cb->initScalar(Q); if (Cb->cstat != CD_IMPLICIT) {
         gStore.save_CData(FL,*Cb); } 
      }
   }} 
   else if (loadRC & LB_UPD) { 
      char check_om=Cb->checkOM(); if (check_om) {
      wbstring file;

      q=gStore.get_file_name(F_L,file,Q,"cgd");
      if (q<=0) {
         if (F) wblog(F,L,
            "ERR missing CBUF %s (i=%d)\n%s",STR(Q),q,file.data);
         check_om=0;  
      }
      if (check_om) { 

      CDATA_TQ Cf;
      Cb->Load_CRef(FL,Cf,file.data);

      if (Cf.cmp(*Cb)<0 || Cf.cmpOM(*Cb)<0) { wblog(FL,
         "ERR %s() cstat mismatch with RCStore\n   %s\n=> %s"
         "%N%N    %s%N -> %s%N", FCT, STR_(Cb), STR(Cf),
         STR2(Cb->cstat,'V'), STR2(Cf.cstat,'V'));
      }

      if (Cb->isRefInit()) { Cf.save2(*Cb); } 
      else {
         q=gStore.load_CData(0,0,Q,*Cb);
         if (q<=0) { wblog(FL, 
            "ERR %s() failed to (re)load %s (q=%d)",FCT,STR_(Cb),q);
         }
      }
   }}}

   if (!Cb->isEmpty()) {
      Cb->setuser_BUF_passive();
   }
   else if (F || (CG_VERBOSE>8 && !isnew)) { char s[64];
      snprintf(s,64," CBUF[%03ld] %s (empty, %snew)",
         BUF.size(), STR(Q), isnew ? "":"!");
      if (F) wblog(PF_L,"ERR %s",s);
      wblog(PF_L,"TST %s",s);
   }
   else if (isnew) {
   }

   if (mC > Cb->numOM()) wblog(FL, 
      "WRN %s() OM decreased %d -> %d\n%s",FCT,mC,Cb->numOM(),STR_(Cb));

   if (!Cb->cgd.isEmpty() && !Cb->cstat.isComplete()) { 
      unsigned q=0;
      if (r<=2) { q=3; } else
      if (r==3) {
         if (!Q.t.permitsOM(r) || Cb->cstat==CD_STD3) { q=3; } else
         if (Cb->cstat==CD_STD3_X) { q=1; }
      }
      if (q) { Cb->cstat.setComplete(q); }
   }

   return *Cb;
};

template <class TQ>
size_t CStore<TQ>::reduceMemUsage(const char *F, int L, double mfac) {

   static long tlast=0; {
      struct timespec tp; int i;
      if ((i=clock_gettime(CLOCK_REALTIME,&tp))) { 
         wblog(FL,"WRN cloc_gettime() returned e=%d",i);
         return 0;
      }
      if (tlast+600 > tp.tv_sec) { return 0; } 
      else { tlast=tp.tv_sec; } 
   }

   size_t mfree=wbsys::getMemFree(), mtot=wbsys::getMemTot();

   if (mfac<0) {
      if (mtot>(long(1)<<34))  
           mfac=0.10; 
      else mfac=0.20; 
   }
   else if (mfac>1) { mfac=1; }

   if (mfree>((1.-mfac)*mtot)) {
      return 0;
   }

   wbvector<size_t> mc(BUF.size()); wbperm P;
   size_t m0=0, m1=0, mref=0; unsigned i, k=0, l=0, n=0;

   for (auto it=BUF.begin(); it!=BUF.end(); ++it, ++k) {
      m0+=(mc[k]=it->second.cgd.memSize());
   }

   if (m0<(mfac*mtot)) {
      return 0;
   }

#ifdef QS_USING_OMP
   Wb::ompGuard gLK(CS_buf); 
#endif

   mc.Sort(P); P.Invert(); mref=(0.80*mfac)*mtot; k=0;
   for (; l<mc.len; ++l) {
      if ((m1+=mc[l])>mref) break;
   }

   m1=m0;
   for (auto it=BUF.begin(); it!=BUF.end(); ++it, ++k) { i=P[k];
      if (i>=l) {
         CData<TQ,RTD> &Cb=it->second;
         Cb.Reduce2Ref(F_L); ++n;
         m1-=mc[i]; m1+=(mc[i]=Cb.cgd.memSize());
      }
   }

   if (CG_VERBOSE>5) { double gfac=1/double(size_t(1)<<30);
      wblog(F_L,"MEM Reduce2Ref() %d/%ld @ %.3g/%.3gG (%g; %.3g/%.3gG)",
      n,BUF.size(), (m0-m1)*gfac, m0*gfac, mfac, mfree*gfac, mtot*gfac);
   }

   return (m0-m1);
};

template <class TQ>
void CStore<TQ>::Info(const char *F, int L, char vflag) const {

   size_t mtot=0; unsigned nx=0, k=map3.size();

   if (k) {
      wblog(F_L," *  CStore::map3 got %d symmetr%s",k,k==1?"y":"ies");
      for (auto I0=map3.begin(); I0!=map3.end(); ++I0) {
         const auto &M1 = I0->second; unsigned m2=0;

         for (auto I1=M1.begin(); I1!=M1.end(); ++I1) {
            const auto &M2 = I1->second;
            m2+=M2.size();

         }
         wblog(F_L,"    > %-6s: %d->%d elements",
         I0->first.toStr('t').data, M1.size(),m2);
      }
   }
   else { wblog(F_L," *  CStore::map3 is empty",map3.size()); }

#ifdef QS_USING_OMP
#endif

   k=0;
   for (auto it=BUF.begin(); it!=BUF.end(); ++it, ++k) {
      const qset<TQ> &qs=it->first.qs;
      if (qs!=it->second.qs) {
         if (it->second.isEmpty()) { ++nx; continue; } 
         wblog(FL,"ERR %s() got inconsistent key [%s] <> [%s]",
         FCT,STR(qs), it->second.QStrS().data);
      }
      mtot+=(it->first.memSize() + it->second.memSize());
      if (vflag) printf("  BUF[%03d] %-24s %20ld :: %-34s 0x%lx\n",k,
         STR(it->first ), QHash<TQ>()(it->first),
         STR(it->second), (unsigned long)&(it->second)
      );
   }

   if (!BUF.size()) { wblog(F_L," *  CStore::BUF  is empty"); }
   else {
      wbvec<char> sx(8);
      if (nx)
           sx.catf(0,0,"%d",nx);
      else sx.catf(0,0,"no");
      wblog(F_L,"--> CStore::BUF total of %d entries (%s)",
         BUF.size(), Wb::size2Str(mtot).data);
      if (nx)
           wblog(FL," *  %d empty entr%s",nx,nx!=1?"ies":"y");
      else wblog(FL," *  with no empty entries",nx);
   }
};

template <class TQ, class TD>
void RStore<TQ,TD>::Info(const char *F, int L) const {

   size_t mtot=0; unsigned n=buf.size();

   if (n)
        wblog(F_L," *  RStore::buf got %d entr%s",n,n!=1?"ies":"y");
   else wblog(F_L," *  RStore::buf  is empty");

   for (auto it=buf.begin(); it!=buf.end(); ++it) {
      const QType &t = it->first;
      const genRG_struct<TQ,TD> &B = it->second;
      unsigned i=0; int l=0; n=B.RSet.size();

      PRINTF("\n  Symmetry %s:\n  > ",STR(t));

      for (auto I=B.RSet.begin(); I!=B.RSet.end(); ++I, ++i) {
         const qset<TQ> &qs=I->first;
         mtot+=I->second.memSize();
         if (l>=0) {
            if (l<50) {
               l+=printf("%s(%s)",i?", ":"", STR(qs));
            }
            else { l=-1;
               if (i+1<n) { PRINTF(" ..."); }
            }
         }
      }
      PRINTF(" (%d entries)\n",n);
      PRINTF("  > estimated memory usage: %s\n\n",Wb::size2Str(mtot).data);
   }
};

template <class TQ, class TD>
genRG_base<TQ,TD>& RStore<TQ,TD>::getR(
   const char *F, int L, const QType &t, const TQ *qs
 ) const {

   qset<TQ> J;
   J.wbvector<TQ>::init2ref(t.qlen(),qs);

   genRG_base<TQ,TD> &R=gRS.Buf(t).RSet[J];

   if (R.J.len && R.Z.dim2) { return R; }
   if (!R.isEmpty()) wblog(FL,
      "ERR %s() got partially empty RSet",FCT);

   int q=gStore.load_RSet(0,1,t,J);
   if (q>0) { if (CG_VERBOSE>8 && L)  {
      QSet<TQ> Q; Q.init1(t,J.data); wblog(PFL,
         "(+) RBUF[%03d] %s", gRS.buf[t].RSet.size(), STR(Q)
      );
   }}
   else { wblog(FL,
      "ERR %s() %s irep (%s) not yet in RCStore (e=%d)",
      FCT, STR(t), STR(J), q);
   }
   if (!R.J.len || !R.Z.dim2) { wblog(F_L, 
      "ERR %s() %s irep (%s) not yet generated (%dx%d)",
      FCT, STR(t), STR(J), R.Z.dim1, R.Z.dim2);
   }

   return R;
};

size_t RCStore::rclog(const QType &t,
   const char *F, int L, const char *fmt, ...
 ){
   size_t l=0; Wb::LogException e;
   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt);

   #ifdef QS_USING_OMP
      Wb::ompGuard myLK(rclog_lk,1); 
   #endif

   try { l=vrclog(t,F,L,fmt,args); }
   catch (Wb::LogException &e_) { e=e_; }
   catch (...) { l=-1; }

   if (e.type) { throw(e);  } else
   if (int(l)<0) { ExitMsg("");  }
   return l;
};

size_t RCStore::rclog(
   const QType &t, 
   const char *F, int L, char vflag, const char *fmt, ...
){
   size_t l=0; Wb::LogException e;
   va_list args; Wb::ARGV wd(&args); 
   va_start(args,fmt);

   #ifdef QS_USING_OMP
      Wb::ompGuard myLK(rclog_lk,1); 
   #endif

   try {
     if (vflag<0)
          { l=vrclog(t,0,0,fmt,args); }
     else { l=vrclog(t,F,L,fmt,args); }
     va_end(args); 

     if (vflag) {
        va_start(args,fmt); 

        if (vflag>=0 && F) { vwblogf(stdout,F,L,fmt,args); }
        else {
           wbstring s; unsigned i=0, l=0;
           for (; i<8; ++i) { s.init(s.len ? 2*s.len : 128);
              if (i) { va_end(args); va_start(args,fmt); } 
              l=vsnprintf(s.data,s.len,fmt,args);
              if (l<s.len) { break; }
           }
           if (l>=s.len) wblog(FL,
              "ERR %s() string out of bounds (%d/%d)%N%N'%s'",
              FCT,l,s.len,s.data
           );
           PRINTF("%s",s.data);
           va_end(args);
        }
     }
   }
   catch (Wb::LogException &e_) { e=e_; }
   catch (...) { l=-1; }

   if (e.type) { throw(e);  } else
   if (int(l)<0) { ExitMsg("");  }
   return l;
};

size_t RCStore::vrclog(
   const QType &t,  
   const char *F, int L,
   const char *fmt, va_list args
){
   unsigned l=0;
   const char *f=RC_LOG.data;

   if (t.isUnknown() || t.isAbelian() || (f && !f[0])) { return l; }

   wbstring file;
   int q=0, d1=0, d2=1; Wb::LogException ex;

 { try { char dstr[24]; strcpy(dstr,"(date-string)");
   if (!f) {
      wbvec<char> s1(16);
      file.init(256); char *fd=file.data; f=fd;

      s1.catf(FL,"%s.log",STR2(t,1));

      q=get_rcs_path(F_L,fd,file.len,t,0,0,s1.data,1); 
      if (abs(q)<20) wblog(FL,
         "ERR %s() failed to open RCStore (%d)\n'%s'",FCT,q,fd);
   }

  #pragma omp critical (got_CPTR_TIME)
   { struct stat fs;
     struct tm *T=NULL; time_t now; time(&now);
     if (stat(f,&fs)==0) {
        T=localtime(&(fs.st_mtime)); d1 = T->tm_yday + 366*T->tm_year;
        T=localtime(& now         ); d2 = T->tm_yday + 366*T->tm_year;
        if (d1!=d2) { strftime(dstr,23,"%d-%b-%Y %T", T); }
     }
     else {
        T=localtime(&now);
        strftime(dstr,23,"%d-%b-%Y %T", T);
     }
   }

   if (!fmt || !fmt[0]) { return (q==-20 ? l : 0); }

   if (F) { wblogf(stdout,F,L,""); }

   FILE *fid=fopen(f,"a"); l=0;
   if (fid) {
      if (d1!=d2) {
         l+=fprintf(fid,"\n>> TODAY %s on %s (%s)\n\n",
         dstr,Wb::hostname().data,myname);
      }
      if (F) { l+=vwblogf(fid,F,L,fmt,args); }
      else   { l+=vfprintf(fid,fmt,args); }
      fclose(fid);
   }
   else wblog(FL,
     "ERR %s() failed to write to log file\n'%s'\n"
     "hint: check environmental variable %s", FCT, f, f==RC_LOG.data ?
     "RC_LOG or --log options" : "RC_STORE"
   );

 } catch (Wb::LogException &e_) { ex+=e_; }
   catch (...) { ++ex; }} 

   ex.report();
   if (int(l)<0) { ExitMsg("");  }

   return l;
};

template <class TQ, class TD>
int x3map<TQ,TD>::contract_x3(const char *F, int L, 
   const CRef<TQ> &A_, const ctrIdx &ica_,
   const CRef<TQ> &B_, const ctrIdx &icb_,
   wbperm &cgp, char xCGR
){
   int rval=0;

   if (!A_.cgb || !B_.cgb) wblog(FL,
      "ERR %s() got empty CRef data (0x%lX, 0x%lX)",FCT,A_.cgb,B_.cgb);

   const CDATA_TQ &A=(*A_.cgb), &B=(*B_.cgb);
   ctrIdx ica(ica_), icb(icb_); { A_.adapt(ica); B_.adapt(icb); }
   QSet<TQ> &Qc=(QSet<TQ>&)c; 

   unsigned i,j, ka,kb, nq=A.t.qlen(), nc=ica.len,
      ra=A.rank(F_L), rb=B.rank(F_L), 
      rc, l=ra+rb;
   const wbperm &pa=A_.cgp, &pb=B_.cgp;

   char sa=(ica.conj!=0 ? -1 : +1),
        sb=(icb.conj!=0 ? -1 : +1), xd=-sa*sb;
   char calc=0, lflag=CG_FIX & cgfix_CID3; 

   wbvector<char> Ma(l); 
   char *ma=Ma.data, *mb=ma+ra; l=0;

#ifndef WB_SKIP_ASSERT
   if (A.t!=B.t || !nq) wblog(F_L,"ERR %s() got symmetry mismatch "
      "(%s, %s)", FCT, A.qStr().data, B.qStr().data);
   if (ra>32 || rb>32) wblog(F_L,
      "ERR %s() unexpected QSpace tensors of rank %d, %d",FCT,ra,rb);

   if (pa.len && pa.isValidPerm(ra)<=0) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pa.len,ra);
   if (pb.len && pb.isValidPerm(rb)<=0) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pb.len,rb);

   if (!nc || nc!=icb.len || nc>ra || nc>rb) wblog(FL,
      "ERR %s() invalid set of contraction indices (%d/%d; %d,%d)",
      FCT,ica.len,icb.len,ra,rb);
   if (A.qdir.len!=ra) wblog(FL,
      "ERR %s() invalid qdir set (A: %d/%d)",FCT,A.qdir.len,ra);
   if (B.qdir.len!=rb) wblog(FL,
      "ERR %s() invalid qdir set (B: %d/%d)",FCT,B.qdir.len,rb);
#endif

   for (i=0; i<nc; ++i) { j=ica[i];
      if (j>=ra) wblog(F_L, 
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,ica.len);
      if (!ma[j]) { ma[j]=i+1; } else wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,ica.len
      );
   }
   for (i=0; i<nc; ++i) { j=icb[i];
      if (j>=rb) wblog(F_L, 
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,icb.len);
      if (!mb[j]) { mb[j]=i+1; } else wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,icb.len);
      if (A.qdir[ica[i]] != xd*B.qdir[j]) wblog(FL,
         "ERR %s() matching in/out pairs of indices required\n"
         "A: %s @ [%s]\nB: %s @ [%s] having xd=%d",FCT,
         STR(A), STR(ica), STR(B), STR(icb),xd
      );
   }

   i=a.RefInit(0,0,A,lflag); 
   j=b.RefInit(0,0,B,lflag); 

   if (i || j) { 
      if ((xCGR&3)<=1) { 
         wblog(FL,"TST %s() unsetting test flag (%d)",FCT,xCGR);
         xCGR &= ~1; 
      }
   }
   calc=(xCGR&3? 0:1); 

   c.init(); cgp.init(); pab.init(); cgb=NULL;
   rtype=CR_DEFAULT;

   ka=ra-nc;
   kb=rb-nc; rc=ka+kb;

   Qc.t=A.t; Qc.qdir.init(rc); Qc.qs.init(nq*rc); 

   if (rc>32) wblog(FL, 
      "ERR %s() unexpected QSpace rank (rc=%d)",FCT,rc);

   if (ka) { for (i=0; i<ra; ++i) { j=pa.at_(i);
      if (!ma[j]) { ma[j]=-(++l); }
   }}
   if (kb) { for (i=0; i<rb; ++i) { j=pb.at_(i);
      if (!mb[j]) { mb[j]=-(++l); }
   }}

   if (l!=rc) wblog(FL,"ERR %s() %d != %d + %d",FCT,l,ka,kb);

   wbperm P0(rc);
   wperm_t *p0=P0.data;

   const char *da=A.qdir.data, *db=B.qdir.data; char *dc=Qc.qdir.data;
   const TQ *qa=A.qs.data, *qb=B.qs.data; TQ *qc=Qc.qs.data;
   int isz=-1; l=0;

   if (ka) { for (i=0; i<ra; ++i) { if (ma[i]<0) {
      dc[l]=sa*da[i]; p0[l++]=-ma[i]-1;
      MEM_CPY<TQ>(qc, nq, qa+i*nq); qc+=nq;
   }}}
   if (kb) { for (i=0; i<rb; ++i) { if (mb[i]<0) {
      dc[l]=sb*db[i]; p0[l++]=-mb[i]-1;
      MEM_CPY<TQ>(qc, nq, qb+i*nq); qc+=nq;
   }}}

   Qc.Sort(&pab);

   if (!A.cgd) { 
      if (!A.NP_zero()) wblog(FL,"ERR %s() got empty A.cgd\n%s",FCT,STR(A));
      isz=10*rc+5; 
   }
   else if (!B.cgd) {
      if (!B.NP_zero()) wblog(FL,"ERR %s() got empty B.cgd\n%s",FCT,STR(B));
      isz=10*rc+6; 
   }
   else {
      isz=Qc.isZero();
   }

   if (P0.extras() || pab.extras(2) || (pab.len && pab.len!=rc)) wblog(FL,
      "ERR %s() unexpected pab=%s, P0=%s @ rc=%d",FCT,STR(pab),STR(P0),rc);

   cgp=pab; 
      if (!pab.len) { cgp.WBPERM::init(rc);
           for (i=0; i<rc; ++i) { cgp[p0[    i ]]=i; }}
      else for (i=0; i<rc; ++i) { cgp[p0[pab[i]]]=i; }
   cgp.compact(); 

   char largeD=0,  
        fOM_C=0,   
        done=0;    

   size_t D1=(1<<14), DX=(1<<20); 
   double cmax=0;
   RTD cfac=1;

   wbsparray<RTD> c3;

   if (A.cgd.D.len>DX) { largeD|=1; } 
   if (B.cgd.D.len>DX) { largeD|=2; } 
   if (CG_VERBOSE>5  ) { largeD|=8; } 

   if (!rc) {
#ifdef QS_USING_OMP
      CG::Guard LK2(FL,A,B);
#endif
      if (A_.cgb!=B_.cgb || ra!=nc || rb!=nc) wblog(FL,"ERR %s() "
         "invalid contraction to scalar (%d/%d/%d)\n%s  %p\n%s  %p",
         FCT,ra,rb,nc, STR_(A_.cgb), A_.cgb, STR_(B_.cgb), B_.cgb);

      if (calc) { 
         try {
            if (A.isRefInit()) A_.LoadRef(FL); 
            if (B.isRefInit()) B_.LoadRef(FL);
         }
         catch (...) {
            throw Wb::LogException(ERR,RC_LOAD_ERR_AB_MSG);
         }

         A.cgd.SPARR_RTD::contract(F,L,ica,B.cgd,icb,c3);
         cmax=double(c3.aMax());
         if (cmax<CG_EPS1) {
            rval=-2; 
         }
         c3.toFull(X3);

         X3.Reshape(A.getOM(),B.getOM(),size_t(1)); 
         X3.SkipTiny(CG_SKIP_REPS);
         x3.initT(X3);

         rtype=CR_CTR_SCALAR; cgb=NULL;

         if (largeD&7 || a.t.sub>4) { 
            if (largeD>8) wblog(FL,
               "FINISHED %s() contracted to scalar (@ %.3g)",
               FCT, x3.data ? double(x3[0]) : -1);

            if ((largeD&1) && A.cgd.D.len>D1) { gCS.Reduce2Ref(FL,a); }
            if ((largeD&2) && B.cgd.D.len>D1) { gCS.Reduce2Ref(FL,b); }
            doflush();
         }
      }
      return rval; 
   }

   if (xCGR&2) {  
      return isz; 
   }

 # ifdef QS_USING_OMP
   CG::Guard LK3;
 # endif

   cfac=(isz>0? 0:1);

   if (calc) {

      CDATA_TQ *Cb=&gCS.getBUF(0,0,Qc,
          !isz || (isz<0 && rc<=3) ? LB_GEN 
       : ( isz<0? LB_LOAD : 0 ));           

      if (Cb->NP_zero()) {
         if (!isz || !Cb->cgd.isEmpty()) wblog(FL, 
            "ERR %s() %s with isz=%d",FCT,STR_(Cb),isz);
         isz=10*rc+7; 
         cfac=0;
      }

   if (isz<=0) { 
      if (Wb::envFullOM) {
         if (A.isComplete()!=3 || B.isComplete()!=3) { wblog(FL,
            "ERR %s() got incomplete OM with input\n"
            "  A: %-32s @ %-4s (%d)\n  B: %-32s -> %-4s (%d)",FCT,
            STR(A), STR(ica), A.isComplete(7),
            STR(B), STR(icb), B.isComplete(7));
         }
      }

    # ifdef QS_USING_OMP
      if (Wb::envFullOM) { const QSet<TQ>* Q1[1]={&Qc};
         LK3.acquire(FL,Q1,1); 
      }
      else {
         LK3.acquire(FL,A,B,Qc); 
      }
    # endif

      if (A  .isRefInit()) { A_ .LoadRef(FL); } 
      if (B  .isRefInit()) { B_ .LoadRef(FL); }
      if (Cb->isRefInit()) { Cb->LoadRef(FL); } 

      if (CG::thread_fOM.me()) {  
         int iA=A.isComplete(), iB=B.isComplete();
         if (!iA || !iB) { wblog(FL,
            "ERR %s() incomplete OM in input\nA %d : %s\nB %d : %s",
            FCT,iA,STR(A),iB,STR(B));
         }
      }

      fOM_C=Cb->hasFullOM(FL);

      if (!fOM_C) wblog(FL,"ERR %s() got %s",FCT,STR_(Cb));

      if (Cb->isEmpty() && !CG::thread_fOM.me()) {
      if (Cb->gotuser_BUF_active()<=1 && fOM_C<32 && Wb::envFullOM) {

         if (fOM_C!=-1) wblog(FL,"WRN %s() got fOM_C=%d",FCT,fOM_C);

         Cb->setuser_BUF_active(3);    
         char q=Qc.generateFullOM(FL); 
         Cb->setuser_BUF_passive();

         if (A.isRefInit()) { A_.LoadRef(FL); } 
         if (B.isRefInit()) { B_.LoadRef(FL); } 

         if (q>=0) { isz=0; }
         else {
            isz=10*rc+4;
            Cb->init(Qc);
            Cb->cstat.setComplete(5); 
         }
         fOM_C=Cb->hasFullOM(FL); 
         if (q<0 && (!Cb->NP_zero() || fOM_C<32)) { wblog(FL,
            "ERR %s() %s",FCT,STR2_(Cb,3));
         }
      }}
   }

   if (Cb->cgd.D.len>DX) { largeD|=4; } 

   if (isz<=0 && fOM_C>0 && (nc<=1 || rc<=3 || !ka || !kb)) {
      unsigned ma=A.numOM(), mb=B.numOM(), mc=Cb->numOM();

      if ( A.reportRefInit(FLF,"A") || 
           B.reportRefInit(FLF,"B") ) wblog(FL,"WRN %s() isz=%d",FCT,isz);

      if (A.cgd.D.len<=B.cgd.D.len && kb) { 
         wbsparray<RTD> BC;
         wbperm Pbc; 
         ctrIdx kb,kcb,iA;

         ica.extend2Perm( 
            ra,ma, icb,rb,mb, pab,mc,
            Pbc,kb,kcb);
         iA.Index(ra);

         B.cgd.SPARR_RTD::contract(F_L,kb,Cb->cgd,kcb, BC,Pbc);
         A.cgd.SPARR_RTD::contract(F_L,iA,BC,iA,c3);

         c3.Reshape( ma?ma:1, mb?mb:1, mc?mc:1 );
      }
      else {
         wbsparray<RTD> AC;
         wbperm Pac; 
         ctrIdx ka,kca,iB;

         ica.extend2Perm(
            ra,ma,icb,rb,mb,pab,mc,
            Pac,ka,kca, 0 
         );
         iB.Index(rb);

         A.cgd.SPARR_RTD::contract(F_L,ka,Cb->cgd,kca, AC,Pac); 
         B.cgd.SPARR_RTD::contract(F_L,iB,AC,iB,c3);

         c3.Reshape(mb?mb:1,ma?ma:1,mc?mc:1).Permute("213");
      }

      cmax=double(c3.aMax());
      c3.toFull(X3); 

      double e2=double(X3.SkipTiny(CG_SKIP_REPS));

      if (e2) { gStore.rclog(A.t, FL, largeD&7 && CG_VERBOSE>6 && F,
         " *  %s() x3 @ %.3g",FCT,std::sqrt(e2)); }

      x3.initT(X3); 

      c.RefInit(0,0,*Cb);
      cgb=Cb; done|=1;
   }
   else if (isz<=0 && Cb->isEmpty() && (ra<=2 || rb<=2)) {
      unsigned ma=A.Mdims(), mb=B.Mdims();
      double e;

      cdata__ x2; wbperm Pab;
      ctrIdx icc; icc.Index(rc);

    # ifdef QS_USING_OMP
      if (!LK3) { 
         if (Wb::envFullOM) { const QSet<TQ>* Q1[1]={&Qc};
                LK3.acquire(FL,Q1,1);   } 
         else { LK3.acquire(FL,A,B,Qc); } 
      }
    # endif

      Cb->init(Qc);
      if (Wb::envFullOM)
           { Cb->cstat_init_(CD_FROM_CTR); }
      else { Cb->cstat.init (CD_FROM_CTR); }

      if (ma && nc==1) {
         if (ra<3 || mb) wblog(FL,
            "ERR %s() got ma=%d with %s, mb=%d",FCT,ma,SSTR(A),mb);
         Pab.initMove(ka,rc,rc+ma+mb);
      }
      Pab.Permute(pab);

      cfac=A.cgd.contract(F_L,ica,B.cgd,icb, Cb->cgd, nullptr, &Pab);

      Cb->cgd.contract(F_L,icc,Cb->cgd,icc,x2);

      e=std::sqrt(double(x2.SkipTiny()));
      cmax=fabs(double(cfac));

      if (e>CG_SKIP_EPS2 || !x2) wblog(FL,
         "ERR %s() invalid x2 (e=%g, %s)",FCT,SSTR(x2));
      if (!x2.isProptoId()) wblog(FL,
         "ERR %s() invalid CGT (nrm=%g; %s)",FCT,double(x2.D[0]),SSTR(x2));
      if (cmax<CG_EPS2) wblog(FL,
         "ERR %s() got small cfac=%g",FCT,cmax);
      if (e) gStore.rclog(A.t, FL, largeD&7 && CG_VERBOSE>6 && F,
         " *  %s() x3 @ %.3g",FCT,e);

      if (ra==2) 
           { x2.wbsparray::Permute("312"); }  
      else { x2.wbsparray::Permute("132"); }  

      if (x2.SIZE.len!=3) wblog(FL,"ERR %s",SSTR(x2));
      (x2*=cfac).toFull(X3);
      x3.initT(X3);

      i=A.cstat.isComplete(7);
      j=B.cstat.isComplete(7);

      if (i && j) {
         Cb->cstat.setComplete(MIN(i,j));
      }
      else if (Wb::envFullOM) { 
         if (!i) wblog(FL,"WRN %s() incomplete OM for A: %s",FCT,STR(A));
         if (!j) wblog(FL,"WRN %s() incomplete OM for B: %s",FCT,STR(B));
      }

      if (Cb->gotuser_BUF()) {
         gStore.rclog(Cb->t, PFL,
            (CG_VERBOSE>8 || (CG_VERBOSE>6 && Cb->rank(FL)<4)) && F,
            "[+] CBUF[%d] new %s",gCS.BUF.size(),STR_(Cb));
         gStore.save_CData(FL,*Cb); 
      }

      c.RefInit(0,0,*Cb);
      cgb=Cb; done|=2;

   }
   else if (isz<=0 && (ra<=2 || rb<=2)) {
      if (Wb::envFullOM) 
      wblog(FL,"WRN %s() isz=%d, ra=%d, rb=%d",FCT,isz,ra,rb);
   }

   if (isz<=0 && largeD>8) { 
      Wb::MemStat(FL);
      wblog(FL,"START %s()\n"
      "    %-35s @ %s (%s)\n    %-35s @ %s (%s)\n    %s",FCT,
           STR(A), STR(ica), SSTR(A.cgd),
           STR(B), STR(icb), SSTR(B.cgd), STR(Qc));
      doflush();
   }

   if (isz<=0 && !done) { 

      if (isz<=0 && ( 
          A.reportRefInit(FLF,"A") ||
          B.reportRefInit(FLF,"B") )) wblog(FL,"WRN %s() isz=%d",FCT,isz);

      unsigned ma=A.numOM(), mb=B.numOM(), mc=0;

      if (!ma) { ma=1; if (!A.cgd.SIZE.len && !A_.isDiagCSC())
         wblog(FL,"ERR %s() got empty A: %s",FCT, STR(A)); }
      if (!mb) { mb=1; if (!B.cgd.SIZE.len && !B_.isDiagCSC())
         wblog(FL,"ERR %s() got empty B: %s",FCT, STR(B)); }

      if (isz<=0) { unsigned irep=0, nrep=1;
         wbarray< wbvector<RTD> > ww(ma,mb);

         double afac;

      for (; irep<nrep; ++irep) { 
         wbIndex IA,IB; unsigned ia=-1, jb=-1;
         cdata__ ai,bj,cx;

         l=0; done|=(1<<(irep+1));

         if (irep) {
            ma=A.numOM(); if (!ma) { ma=1; }
            mb=B.numOM(); if (!mb) { mb=1; }
            ww.Resize(ma,mb);
         }

         a.RefInit(0,0,A); 
         b.RefInit(0,0,B); 

         while (B.getCG_set(FL,IB,bj)) { ++jb; ia=-1; IA.init();

           #ifdef LD_CLEBSCH_QS
            if (largeD>8) {
               Wb::MemStat(FL,'m');
               wblog(FL,"om: jb=%2d | m=(%d,%d,%d)",jb+1,ma,mb,mc);
            }
           #endif

         while (A.getCG_set(FL,IA,ai)) { ++ia; 
            cx.init(); 
            cfac=ai.contract(F,L,ica,bj,icb,cx, nullptr, &pab);

            afac=fabs(double(cfac));

            if (cmax<afac) { cmax=afac; }

            if (afac<CG_EPS2) { ++l; continue; }
            else if (*Cb) {
               if ((*Cb)!=Qc) wblog(FL,
                 "ERR %s() got QSet mismatch\n%s <> %s",FCT,STR_(Cb),STR(Qc));
               mc = Cb->getOM();
            }
            else {
               Cb = &gCS.getBUF(0,0,Qc,LB_GEN);
               mc = Cb->getOM();

               if (!mc || Cb->isEmpty()) {
                  (*Cb)=Qc;
                  if (CG::thread_fOM.me()) {
                     Cb->cstat_init_(CD_FROM_CTR); if (xCGR>3) {
                        Cb->cstat.setComplete(3);  
                     }
                  }
                  else {
                     if (Wb::envFullOM) { wblog(FL,
                        "WRN %s() got thread_fOM=0 with QS_FULL_OM (tid=%d/%d)"
                        "-> %s [init DEC]",FCT, omp_get_thread_num(),
                        omp_get_num_threads(), STR_(Cb));
                     }
                     Cb->cstat.init(CD_FROM_DEC);
                  }

                  gStore.rclog(A_.cgb->t, PFL, CG_VERBOSE>6 && F,
                     "[+] CBUF[%03d] new/dec: %s #%05X",
                     gCS.BUF.size(), STR(Qc), Cb->cstat.cID
                  );
               #ifndef WB_SKIP_ASSERT
                  if (Qc.isStd3() && CG_VERBOSE>6) wblog(FL,
                  " *  %s() new %s [%s]",FCT,STR(Qc),BITS(LB_GEN));
               #endif
               }
            }

            if (ia>=ma || jb>=mb) {
               unsigned ma_=A.getOM(); { if (!ma_) ma_=1; }
               unsigned mb_=B.getOM(); { if (!mb_) mb_=1; }
               char vflag=(CG_VERBOSE>2 ? ma_>ww.dim(1) || mb_>ww.dim(2) : 0);

               if (vflag && Wb::envVRB&8) {
                  wblog(FL,"NB! %s() got extended OM space!",FCT);
                  if (ma_!=ma) wblog(FL,
                   " A: %s #%05X (%d->%d)",STR(A), A.cstat.cID, ma, ma_);
                  if (mb_!=mb) wblog(FL,
                   " B: %s #%05X (%d->%d)",STR(B), B.cstat.cID, mb, mb_);
               }

               if (ia>=ma_ || jb>=mb_ || ma_<ma || mb_<mb) wblog(FL,
                  "ERR %s() OM out of bounds (%d/%d<-%d; %d/%d<-%d)",
                  FCT,int(ia),ma_,ma,int(jb),mb_,mb);

               ww.Resize(ma_,mb_);

               if (ma==ma_ || ia==0) { ma=ma_; mb=mb_; }
               else if (nrep<irep+2) {
                  if (vflag) { wblog(FL,"TST i_ab=(%d,%d): "
                     "OM_ab = (%d,%d) -> (%d,%d), repeat n=%d -> %d",
                     ia,jb,ma,mb,ma_,mb_,nrep,irep+2);
                  }
                  nrep=irep+2;
               }
            }
            if (cx) {
               if (CG::thread_fOM.me() && xCGR>3 && !Cb->cstat.isComplete()) {
                  Cb->cstat.setComplete(3); 
               }

               Cb->Project(FL,cx,ww(ia,jb)); 

               if (mc<ww(ia,jb).len) { mc=ww(ia,jb).len; }
               ww(ia,jb)*=cfac;
            }
            ++l;

         }}

         if (mc) { cgb=Cb; X3.init(ma,mb,mc);
            for (i=0; i<ma; ++i) {
            for (j=0; j<mb; ++j) { wbvector<RTD> &w=ww(i,j);
            for (l=0; l<w.len; ++l) { X3(i,j,l)=w[l]; }}}

            double e=std::sqrt(double(X3.SkipTiny(CG_SKIP_REPS)));
            if (e>CG_EPS) {
               char s[32]; snprintf(s,32,"%s() x3 @ %.3g",FCT,e);
               gStore.rclog(A.t, FL, CG_VERBOSE>6 && F,
               e<CG_EPS1 ? " *  %s" : "WRN %s",s);
            }
            x3.initT(X3);
         }
         if (nc==1 && fOM_C==-1) {
            if (x3.nnz(1e-12)!=x3.ndiag()) { wblog(FL,
               "WRN %s() %s\ngot non-diagonal x3 (%s: %d/%d) t=%d:%d",
               FCT,STR_(Cb),SSTR(x3),x3.nnz(1e-12),x3.ndiag(),
               omp_get_thread_num(), CG::thread_fOM.me());
            }
         }
      }
      }

      if (cgb) {
         c.RefInit(0,0,*cgb); 
      }
      else if (mc || cmax>=CG_EPS2) { 
         wblog(FL,"ERR %s() %d 0x%lx %.3g",FCT,mc,cgb,cmax);
      }

      #ifndef WB_SKIP_ASSERT
      { int e=0;
        if (ma!=a.getOM(FL)) { e|=1; wblog(FL,
           "WRN %s() OM mismatch (%d/%d)\nA: %s\na: %s",
           FCT, ma, a.getOM(FL),STR(A),STR(a)); }

        if (mb!=b.getOM(FL)) { e|=2; wblog(FL,
           "WRN %s() OM mismatch (%d/%d)\nB: %s\nb: %s",
           FCT, mb, b.getOM(FL),STR(B),STR(b)); }

        if (mc && mc!=c.getOM(FL)) { e|=4; wblog(FL,
           "WRN %s() OM mismatch (%d/%d)\nC: %s\nc: %s",
           FCT, mc, c.getOM(FL),STR_(cgb),STR(c)); }

        if (e) wblog(FL,"ERR %s() e=%d (isz=%d)",FCT,e,isz);
      }
      #endif

   } 
   } 

   if (calc || isz>0) { 
      if (cmax<CG_EPS2) {
         if (calc) {
         if (largeD&7 || a.t.sub>4) { 
            if (largeD>8) {
               if (isz>0) wblog(FL,
               "FINISHED %s() yields zero (by QSet)\n%s",FCT,STR(Qc));
               else wblog(FL,
               "FINISHED %s() contracted to zero (@ %.4g)",FCT,cmax);
            }
            if ((largeD&1) && A.cgd.IDX.dim1>4096) { gCS.Reduce2Ref(FL,a); }
            if ((largeD&2) && B.cgd.IDX.dim1>4096) { gCS.Reduce2Ref(FL,b); }
            doflush();
         }
         else  {
            if (isz>0) {
               if (CG_VERBOSE>8) wblog(FL,  
               " *   %s() %s to zero [z=%d,%d]",FCT,STR(Qc),isz,xCGR);
            }
            else if (CG_VERBOSE>8) { wblog(FL, 
             " *  %s() %s zero @ %.3g (%d/%d)",FCT,STR(Qc),cmax,isz,xCGR);
            }
         }}

         if (!isz && nc<2) { wblog(FL, 
            "WRN %s() got cmax=%g with isz=%d, rc=(%d+%d-2*%d)=%d legs",
            FCT,cmax,isz,ra,rb,nc,rc);
         }
         if (done && isz>0) { wblog(FL, 
            "WRN %s() got done=%d with isz=%d (cmax=%g)",
            FCT,done,isz,cmax);
         }

         c.init();   
         pab.init(); 
         cgb=NULL; rtype=CR_CTR_ZERO;

         x3.init(); 
         #ifdef QS_USING_MPFR
            X3.init();
         #endif

         if (isz>0) { zflag=isz; } 

         return rval;
      }

      if (cmax<1e-5) wblog(FL,"WRN %s() got small cmax=%.4g",FCT,cmax);
   }

  #ifndef WB_SKIP_ASSERT
   if (cgb) { int e=0;
      if ((e=cgb->cstat.inValid())) wblog(FL,"ERR %s() %s (e=%d)\n%s",
         FCT,STR_(cgb),e,STR2(cgb->cstat,'V'));
      if (!cgb->gotuser_BUF()) { wblog(FL,
         "WRN %s() missing user_BUF flag",FCT); 
      }
   }
  #endif

   if (largeD&7 || a.t.sub>4) { 
      if (largeD>8) wblog(FL,
         "FINISHED %s() %s [done=%d]",FCT,x3Str().data,done);
      if ((largeD&1) && A.cgd.IDX.dim1>4096) { gCS.Reduce2Ref(FL,a); }
      if ((largeD&2) && B.cgd.IDX.dim1>4096) { gCS.Reduce2Ref(FL,b); }

      if (cgb) {
         if (cgb->cgd.IDX.dim1>4096) { gCS.Reduce2Ref(FL,c); }
         else if (largeD&7 && CG_VERBOSE>3) {
            wblog(FL,"TST %s() %s [%ld/%ld]",FCT,
            STR_(cgb),cgb->cgd.IDX.dim1,cgb->cgd.D.len);
         }
      }
      doflush();
   }

   if (largeD>8) Wb::MemStat(FL);

   return rval;
};

template <class TQ, class TD>
void x3map<TQ,TD>::print(
   const char *F, int L, const char *istr, char vflag,
   const ctrIdx *ica, const ctrIdx *icb 
 ) const {

   wbvector<size_t> S;
   wbperm P;

   wblog(F_L,"I/O x3map() data %s",istr? istr:"");

   PRINTF("\nx3map data"
     "\n> QSets\n     a: %-36s%s%s\n   * b: %-36s%s%s\n   = c: %s\n",
      STR2(a,vflag), ica? " @ ":"", ica? STR_(ica):"",
      STR2(b,vflag), icb? " @ ":"", icb? STR_(icb):"",
      STR2(*this,vflag));

   PRINTF("\n> cstat\n     a: %s\n   * b: %s\n   = c: %s\n\n",
      STR2(a.cstat,vflag), STR2(b.cstat,vflag), STR2(c.cstat,vflag));

   Wb::hpsort(S=x3.SIZE,P,-1); 

   if (P.len>=2 && P[0]>P[1] && x3.SIZE[P[0]]>x3.SIZE[P[1]]) {
      SWAP(P[0],P[1]); 
   }

   PRINTF("> |x3|=%s %s\n",RATS(double(x3.norm())),SSTR(x3));
   if (x3.numel()<64) { PRINTF("\n");
      if (!P.relevant()) { x3.printdata("x3"); }
      else {
         char sx[16]; snprintf(sx,16,"x3_%s",STR(P));
         x3.permute(P).printdata(sx);
      }
      PRINTF("\n");
   }
};

template <class TQ, class TD>
wbstring x3map<TQ,TD>::toStr(char vflag) const {

   wbvec<char> s;

   if (!pab.len && !cgb && c.isEmpty('c')) { 
      s.init(128);
      if (x3.isEmpty()) {
         if (rtype==CR_CTR_ZERO)
              s.catf(FL,"contracts to zero"); else
         if (rtype==CR_DEFAULT)
              s.catf(FL,"empty x3map");
         else s.catf(FL,"empty x3map (having rtype=%s)",STR(rtype));
         return s.data;
      }

      if (rtype==CR_CTR_SCALAR) {
         s.catf(FL,"contracts to scalar");
         if (x3.numel()==1) {
            s.catf(FL," (x3=%.4g)",double(x3.data[0]));
         }
         else {
            double q=x3.sparsity();
            s.catf(FL," (x3: %s / %.4g",SSTR(x3), double(x3.norm2()));

            if (q<1)
                 { s.catf(0,0," @ %.0f%%)",100*q); }
            else { s.append(0,0,')'); }
         }
         return s.data;
      }
   }

   s.init(1024);

   s.catf(FL,"%s",STR((QSet<TQ>&)c));  if (pab) {
   s.catf(FL,", pab=%s",STR(pab));     }

   if (vflag) {
      if (vflag>2) {
         if (vflag=='v') vflag=1; else
         if (vflag=='V') vflag=2;
         else wblog(FL,"ERR %s() invalid vflag=%d",FCT,vflag);
      }
      s.catf(FL," => %s",x3Str(vflag<4 ? vflag-1 : 3).data);
   }

   return s.data;
};

template <class TQ, class TD>
wbstring x3map<TQ,TD>::x3Str(char vflag) const {

   if (x3.isEmpty() && !rtype) { return "x3=[]"; }

   if (vflag>8) {
      if (vflag=='v') { vflag=1; } else
      if (vflag=='V') { vflag=2; }
      else {
         wblog(FL,"WRN %s() got vflag=%s",FCT,cSTR(vflag));
         vflag=3;
      }
   }

   size_t n=x3.numel();
   wbvec<char> sx;

   if (!n) { sx.init(16);
      sx.catf(FL,"%s (%s)", rtype==CR_CTR_ZERO? "zero":"x3=[]",
      rtype.tostr());
   }
   else if (n==1 || vflag<=1) { sx.init(64);
      if (n==1)
           { sx.catf(0,1,"x3=%s",RATS(x3.data[0])); }
      else { sx.catf(0,1,"|x3|=%s",RATS(double(x3.norm()))); }

      if (rtype) {
         sx.catf(0,0," %s",
         vflag && rtype==CR_CTR_SCALAR? "scalar" : rtype.tostr());
      }

      if (x3.SIZE.len!=3 || n!=1) {
      sx.catf(0,1," (%s)", SSTR(x3)); }
   }
   else {
      if (n<=4 || vflag>3) {
         wbstring dstr=wbvector<TD>(n,x3.data,'r').toStr();
         sx.init(32+dstr.len);
         sx.catf(0,0,"x3=[%s]",dstr.data);
      }
      else {
         wbstring dstr=wbvector<TD>(3,x3.data,'r').toStr();
         sx.init(32+dstr.len);
         sx.catf(0,0,"x3=[%s, ...]",dstr.data);
      }

      if (rtype) {
      sx.catf(0,0," %s",rtype.tostr()); }
      sx.catf(FL," (%s)",SSTR(x3));
   }

   return sx.data;
};

template <class TM> 
template <class TQ, class TD> inline
cgc_contract_id<TM>::cgc_contract_id(
   const CData<TQ,TD> &A, const ctrIdx &ica,
   const CData<TQ,TD> &B, const ctrIdx &icb
){
   unsigned n=256, ra=A.rank(FL), rb=B.rank(FL), nc=ica.len;
   unsigned l, m=sizeof(TM); wbstring s1;
   wbvector<unsigned> S;

   wbvec<char> s_(n);
   char *s=s_.data;

   if (A.t!=B.t) wblog(FL,"ERR %s() qtype mismatch: "
      "%s <> %s",FCT, A.qStr().data, B.qStr().data);
   if (!A.qs.len || !B.qs.len || !nc || nc!=icb.len)
       wblog(FL,"ERR %s() got empty data set",FCT);
   if (A.qdir.len!=ra || B.qdir.len!=rb) wblog(FL,
      "ERR %s() qdir length mismatch (%d/%d; %d/%d",
      FCT, A.qdir.len, ra, B.qdir.len, rb);
   if (nc>ra || nc>rb) wblog(FL,
      "ERR %s() invalid contraction (%d/%d; %d/%d)",
      FCT, ica.len,ra, icb.len,rb);
   if (!ica.isSorted()) wblog(FL, 
      "ERR %s() got non-sorted ica=%s",FCT,STR(ica));
   if (n<(l=4*(A.qdir.len+B.qdir.len)+16)) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,n,l); 

   ++nc; 

   set_CData(FL,s,n-=nc,A);
   set_vec(FL,s,ica);
   set_val(FL,s,Wb::conj2bool(ica.conj)); 

   l=s-s_.data;
   if (int(l)<2 || l+4*m>s_.len) wblog(FL, 
      "ERR %s() unexpected string length (l=%d+%d/%d)",FCT,l,4*m,s_.len);
   n-=l;

   set_CData(FL,s,n-=nc,B,0); 
   set_vec(FL,s,icb);
   set_val(FL,s,Wb::conj2bool(icb.conj)); 
   l=s-s_.data;

   if (int(l)<2 || l+m>s_.len) wblog(FL, 
      "ERR %s() unexpected string length (l=%d+%d/%d)",FCT,l,m,s_.len);

   if ((l%=m)) { for (; l<m; ++l, ++s) { *s=0; } }
   if ((l=s-s_.data)%m) wblog(FL,"ERR %s() l=%d / m=%d",FCT,l,m);

   this->init(l/m,(TM*)s_.data);
};

template <class TM> 
template <class TQ, class TD> inline
void cgc_contract_id<TM>::extract(
   CData<TQ,TD> &A, ctrIdx &ica,
   CData<TQ,TD> &B, ctrIdx &icb,
   char xflag 
 ) const {

   if (!this->len) wblog(FL,"ERR %s() object not yet initialized",FCT);
   const char *sx = (const char*) this->data;

   get_CData(FL,sx,A); A.rank(FL); 
   get_vec(FL,sx,ica);
   get_val(FL,sx,ica.conj);

   get_CData(FL,sx,B,0); B.t=A.t; B.rank(FL);
   get_vec(FL,sx,icb);
   get_val(FL,sx,icb.conj);

   if (ica.anyGT(20) || icb.anyGT(20)) { wblog(FL, 
      "ERR %s() unexpected contraction index (%s; %s)",
      FCT, STR(ica), STR(icb));
   }
   if (!xflag) { return; } 

 # ifdef QS_USING_OMP
   CG::Guard xLK(FL,*this);
 # endif
   const x3map<TQ,double> &M=gXS.getXBUF(0,0,*this);

   if (!M.isEmpty()) { 
      const wbvector<size_t> &Sa=M.a.cgd.SIZE, &Sb=M.b.cgd.SIZE;
      unsigned i, ra=A.cgd.SIZE.len, rb=B.cgd.SIZE.len;

      A.cstat=M.a.cstat;
      B.cstat=M.b.cstat;

      if (ra<Sa.len) { A.cgd.SIZE=Sa; } else
      if (ra>(i=Sa.len) || A.cgd.D.len) { 
         if (i || xflag!='l') wblog(FL,   
         "ERR %s() unexpected OM setting (%d/%d/%d; %d)\nA: %s %s <> %s %s",
         FCT,ra,i,Sa.len, A.cgd.D.len, STR(M.a), SSTR(M.a), STR(A), SSTR(A));
      }

      if (rb<Sb.len) { B.cgd.SIZE=Sb; } else
      if (rb>(i=Sb.len) || B.cgd.D.len) { 
         if (i || xflag!='l') wblog(FL,   
         "ERR %s() unexpected OM setting (%d/%d/%d; %d)\nB: %s %s <> %s %s",
         FCT,rb,i,Sb.len, B.cgd.D.len, STR(M.b), SSTR(M.b), STR(B), SSTR(B));
      }
   }
};

template <class TQ, class TD>
x3map<TQ,TD>& X3Map<TQ,TD>::getXBUF( 
   const char *F, int L, 
   const cgc_contract_id<MTI> &idc
){

 # ifdef QS_USING_OMP
   Wb::ompGuard xLK(XS_buf);
 # endif

   auto it=XBUF.find(idc);
   if (it==XBUF.end()) { if (!F) {
         auto &x=XBUF[idc];
         #if defined(DBG_QSX_BUF) && (DBG_QSX_BUF & 2)
            BFF.blogf(FL," +  %-10s %3ld %p  %s",
            FCT,XBUF.size(),&x,STR2(idc,0));
         #endif
         return x;
      }
      wblog(F,L,"ERR %s() got non-existing entry",FCT);
   }

   #if defined(DBG_QSX_BUF) && (DBG_QSX_BUF & 2)
      BFF.blogf(FL," *  %-10s %3ld %p  %s",
      FCT,XBUF.size(),&(it->second),STR2(idc,0));
   #endif

   if (F && it->second.isEmpty()) {
      wblog(F,L,"ERR %s() got empty XBUF entry",FCT); }
   return it->second;
};

template <class TQ, class TD>
int X3Map<TQ,TD>::contains(const char *F, int L,
   const x3map<TQ,TD> &x, const cgc_contract_id<MTI> &idc) {

   auto it=XBUF.find(idc);
   if (it!=XBUF.end()) {
      if (&x == &it->second) return 1; 
      if (F) wblog(F,L,
         "ERR %s() got x3map object outside XBUF\nthis: %s\nXBUF: %s",
         FCT, STR(x), STR(it->second));
      return -2;
   }
   else {
      if (F) wblog(F,L,"ERR %s() non-existing XBUF entry",FCT);
      return -1;
   }
};

template <class TQ, class TD> 
int X3Map<TQ,TD>::contractCGR( const char *F, int L,
   const CRef<TQ> &A, const ctrIdx &ica_,
   const CRef<TQ> &B, const ctrIdx &icb_, CRef<TQ> &C,
   char xCGR 
) {

   char isa=0; 
   if (!A.cgb) { isa|=1; } else if (A.isAbelian()) { isa|=4; }
   if (!B.cgb) { isa|=2; } else if (B.isAbelian()) { isa|=8; }

   if (!(isa&3)) { 
      if (A.cgb->t!=B.cgb->t) wblog(F_L,
      "ERR %s() got mixed symmetry (isa: %d)\n   %s\n<> %s",
      FCT,BITS(isa), STR_(A.cgb), STR_(B.cgb));
   }

   if (isa) {
      if (!(isa&5) || !(isa&10)) wblog(FL, 
         "ERR %s() inconsistent abelian CRefs (isa: %s)",FCT,BITS(isa));
      if (xCGR&1) { return 0; } 

      if (isa&3) { 
         C.initAbelian();
         if (A.cgw) {
            if (!A.wscalar1()) wblog(FL,
               "ERR %s() invalid abelian CRef %s",FCT,STR(A));
            C.cgw=A.cgw;
         }
         if (B.cgw) {
            if (!B.wscalar1()) wblog(FL,
               "ERR %s() invalid abelian CRef %s",FCT,STR(B));
            if (!C.cgw.data) { C.cgw=B.cgw; }
         }
         return 0;
      }

      if (!A.wscalar() || !B.wscalar()) wblog(F_L,
         "ERR %s() unexpected scalar CRef data (%s / %s)",
         FCT, STR(A.cgw), SSTR(B.cgw));

      C.initAbelian(A.cgw[0]*B.cgw[0]);

      QSet<TQ> Q(F_L, A, ica_, B, icb_);
      Q.checkQ_abelian(F_L);

      Q.Sort(&C.cgp,'i');

      if (!Q.qs && !Q.qdir) {
         C.rtype=CR_CTR_SCALAR;
         C.cgb=NULL; 
      }
      else {
         CDATA_TQ &S=gCS.getBUF(FL,Q); 

         if (S.isEmpty()) {
            S=Q; S.cstat.init(CD_ABELIAN); S.cgd.init(); 
            if (CG_VERBOSE>2)
            wblog(PFL,"[+] add2BUF() new/ctr: %s",S.toStr('v').data);
         }
         else if (S!=Q) wblog(F_L,"ERR %s() "
            "CData inconsistency\n%s <> %s",FCT, STR(S), STR(Q));
         else if (S.getScalar(F_L)!=1) wblog(F_L,
            "ERR %s() invalid scalar cgd\n%s",FCT,STR(S));

         C.cgb=&S; 
      }

      return 0; 
   }

   if (!A.cgb || !B.cgb) wblog(F_L,
      "ERR %s() got empty cref data (0x%lX, 0x%lX)",FCT,A.cgb,B.cgb);

   if (!A.isw2() || !B.isw2()) wblog(FL,
      "ERR %s() invalid cgw (A: %s; B: %s)",FCT,SSTR(A.cgw),SSTR(B.cgw));

   ctrIdx ica(ica_), icb(icb_); wbperm cgp;
   A.adapt(ica);
   B.adapt(icb);

   ica.Sort_(icb); 

   cgc_contract_id<MTI> idc(*A.cgb,ica, *B.cgb,icb); 

#ifdef QS_USING_OMP
   CG::Guard xLK(FL,idc);

#endif

   x3map<TQ,TD> &Cm=getXBUF(0,0,idc);

   int loaded=0;       
   unsigned calc=0;    
   char sq[3]={0,0,0}; 

   int ma=A.wdim1(), Ma=A.getOM(),
       mb=B.wdim1(), Mb=B.getOM();

   int mc=0, Mc=1; 

   if (!ma || ma>Ma || !mb || mb>Mb) wblog(FL,
      "ERR %s() invalid OM data (%dx%d; %dx%d)",FCT,ma,mb,Ma,Mb);

   if (Cm.isEmpty()) {
      if (CG_FIX & cgfix_CTR) wblog(FL,
         "WRN %s() ignore XStore data, recompute in any case",FCT);
      else if (gStore.load_XMap(0,0,idc)>0) { loaded=1; }
   };

   if (xCGR&2) { static size_t count_new=0; 
      if (Cm.isEmpty()) { x3map<TQ,TD> Cx;
         int isz=Cx.contract_x3(FL,A,ica_,B,icb_, cgp, xCGR);
         if (!isz) { PRINTF(           
            "  CTR_PREVIEW[%3ld]  %s_%s %s_%s\n",++count_new,
               STR((QSet<TQ>)(*A.cgb)), STR(ica),
               STR((QSet<TQ>)(*B.cgb)), STR(icb));
            return 1;
         }
      }
      return 0;
   }

   if (Cm.isEmpty()) { calc|=1; } else
   if (!Cm.cgb && !Cm.x3 && Cm.rtype==CR_DEFAULT) { calc|=2; } 

   if (calc) {
      Cm.contract_x3(FL,A,ica_,B,icb_, cgp, xCGR&~3); 
      if (Cm.cgb) {
         mc=Cm.x3.SIZE[2];
         Mc=Cm.cgb->getOM(); 
      }

      if (Cm.rtype!=CR_CTR_ZERO) { 
      char fmt[128]; snprintf(fmt,128,
         "\b[+] XBUF[%%03d] %%s { %s %s %s }\n"
         "    %%s @ %%s->%%s\n    %%s @ %%s->%%s\n--> %%s",
         Cm.a.cstat.cID ? "#%05X":"#%X?",
         Cm.b.cstat.cID ? "#%05X":"#%X?",
         Cm.c.cstat.cID ? "#%05X":"#%X");

      gStore.rclog(A.cgb->t, PFL,
         CG_VERBOSE>8 || (CG_VERBOSE>6 && (
            Cm.x3.numel()>1 || Cm.rtype==CR_CTR_SCALAR ||
            (A.cgb->qdir.len + B.cgb->qdir.len - ica.len - icb.len) > 3
         )),
         fmt, XBUF.size(),
         Cm.x3Str().data, Cm.a.cstat.cID, Cm.b.cstat.cID, Cm.c.cstat.cID,
         STR(A), STR(ica_), STR(ica),
         STR(B), STR(icb_), STR(icb), STR(Cm));
      }

   }
   else {

      x3map<TQ,TD> Cx;
      Cx.contract_x3(FL,A,ica_,B,icb_,cgp,xCGR|1); 

#ifndef WB_SKIP_ASSERT
      if (loaded) { 
         char e1=!Cx.a.sameAs(Cm.a,'L'),  
              e2=!Cx.b.sameAs(Cm.b,'L'),
              e3=(Cx.cgb ? !Cx.c.cstat.sameAs(Cm.c.cstat,'L') : 0);

         if (e1 || e2 || e3) { PRINTF("\n");
            wblog(FL,"WRN x3map::%s() inconsistency [%d%d%d]",FCT,e1,e2,e3);
            wblog(FL,"--- "
               "> A %s @ %s -> %s\n"
               "  B %s @ %s -> %s\n> c %s '%s'",
               STR(A),STR(ica_),STR(ica),
               STR(B),STR(icb_),STR(icb), STR(Cx.c), STR(cgp));

            if (e1) PRINTF("\n a: %-30s %s\n != %-30s %s\n",
               STR(Cx.a),STR2(Cx.a.cstat,'V'),
               STR(Cm.a),STR2(Cm.a.cstat,'V'));
            if (e2) PRINTF("\n b: %-30s %s\n != %-30s %s\n",
               STR(Cx.b),STR2(Cx.b.cstat,'V'),
               STR(Cm.b),STR2(Cm.b.cstat,'V'));
            if (e3) PRINTF("\n c: %-30s %s\n != %-30s %s\n",
               STR(Cx.c),STR2(Cx.c.cstat,'V'),
               STR(Cm.c),STR2(Cm.c.cstat,'V'));

            if (CG_FIX)
                 { calc|=16; PRINTF("\n"); }
            else { wblog(FL,"ERR %s()",FCT); }
         }
      }
#endif

      if (Cm.rtype!=CR_CTR_ZERO && (
         (QSet<TQ>&)(Cm.c)!=(QSet<TQ>&)(Cx.c) || Cm.pab!=Cx.pab)
      ){ wblog(FL,
          "ERR %NX3Map::%s() -> C QSet inconsistency%N" 
          "   %-50s @ %-3s -> %s%N"
          "   %-50s @ %-3s -> %s%N=> %s <> %s%N",FCT, 
          STR(A), STR(ica_), STR(ica),
          STR(B), STR(icb_), STR(icb), STR(Cm), STR(Cx));
      }

      if (Cm.cgb && !Cm.x3.isEmpty()) {
         if (Cm.x3.SIZE.len!=3) wblog(FL,
            "ERR %s() invalid x3 (%s)",FCT,SSTR(Cm.x3));
         mc=Cm.x3.SIZE[2];   
         Mc=Cm.cgb->getOM(); 

         sq[0]=NUMCMP(Ma,int(Cm.x3.SIZE.data[0]));
         sq[1]=NUMCMP(Mb,int(Cm.x3.SIZE.data[1]));
         sq[2]=NUMCMP(Mc,mc); 

         if (sq[2]<0) wblog(FL,
            "ERR %s() got x3 with larger OM !? %s (%dx%dx%d)%N%N   "
            "A: %-25s| a: %-25s @ %-3s%N   "
            "B: %-25s| b: %-25s @ %-3s%N   "
            " %-25s => c: %-25s%N   %s%N", FCT,SSTR(Cm.x3),Ma,Mb,Mc,
            STR_(A.cgb), STR(Cm.a), STR(ica),
            STR_(B.cgb), STR(Cm.b), STR(icb),
            "", STR(Cm.c), STR2(Cm.cgb->cstat,'V')
         );

         if (sq[0]>0 || sq[1]>0 || sq[2]>0) { calc|=4; }
      }
      else {
         int qa=A.cgb->cstat.cmp(0,0,Cm.a.cstat),
             qb=B.cgb->cstat.cmp(0,0,Cm.b.cstat);

         if (qa<-1 || qb<-1) {
            if (CG_FIX) {
               wblog(FL,"WRN %s() using CG_FIX=%d",FCT,CG_FIX);
               calc|=16; 
            }
            else { char istr[8];
               if (qa<-1 && qb<-1) strcpy(istr,"A and B"); else
               if (qa<-1) strcpy(istr,"A"); else strcpy(istr,"B");

               PRINTF("\n A: %s @ %s\n a: %s\n\n    %s\n <> %s\n",
                  STR_(A.cgb),STR(ica),STR(Cm.a),
                  STR2(A.cgb->cstat,'V'), STR2(Cm.a.cstat,'V'));
               PRINTF("\n B: %s @ %s\n b: %s\n\n    %s\n <> %s\n",
                  STR_(B.cgb),STR(icb),STR(Cm.b),
                  STR2(B.cgb->cstat,'V'), STR2(Cm.b.cstat,'V'));
               PRINTF("\n c: %s ; %s\n", STR(Cm.c),
                  Cm.cgb ? STR2(Cm.cgb->cstat,'V') : "(null)");

               wblog(FL,"ERR %s() "
              "got cstat mismatch in %s (loaded=%d)",FCT,istr,loaded);
            }
         }
         else if (qa || qb) {
            unsigned om[4]={
               Cm.a.getOM(), A.cgb->getOM(),
               Cm.b.getOM(), B.cgb->getOM() };

            gStore.rclog(A.cgb->t, FL,
              (CG_VERBOSE>6 || ((qa<0 || qb<0) && CG_VERBOSE>2)) && F,
               " *  %s() recalc (A%s om=%d->%d; B%s om=%d->%d)",FCT,
               qa>0 ? " newer" : (qa ? " older!?" : ""), om[0], om[1],
               qb>0 ? " newer" : (qb ? " older!?" : ""), om[2], om[3]
            );
            if (qa<0) gStore.rclog(A.cgb->t, FL, CG_VERBOSE>2 && F,
               "WRN A: %s\n   %s\n-> %s",
               STR(A), STR2(Cm.a.cstat,'V'), STR2(A.cgb->cstat,'V')
            );
            if (qb<0) gStore.rclog(B.cgb->t, FL, CG_VERBOSE>2 && F,
               "WRN B: %s\n   %s\n-> %s",
               STR(B), STR2(Cm.b.cstat,'V'), STR2(B.cgb->cstat,'V')
            );
            calc|=8;
         }

         if (Cm.x3.SIZE.len) { mc=1; 
            if (Cm.x3.SIZE[2]!=1) wblog(FL,"ERR %s() "
               "got unexpected x3 (%s)",FCT,SSTR(Cm.x3));
            sq[0]=NUMCMP(Ma,int(Cm.x3.SIZE[0]));
            sq[1]=NUMCMP(Mb,int(Cm.x3.SIZE[1]));

            if (!calc && (sq[0] || sq[1])) wblog(FL, 
               "ERR %s() got cstat inconsistency",FCT);

            if (sq[0]<=0 && sq[1]<=0 && !(calc&16)) { 
               if (calc>7 && CG_VERBOSE>1) wblog(FL,
                  " *  %s() x3: %s (%dx%dx%d) -> unset calc\nA: %s\nB: %s",
                  FCT,SSTR(Cm.x3),Ma,Mb,mc,STR(A),STR(B));
               calc=0;
            }
         }
      }

      if (calc) {

         wbarray<TD> x3_; Cm.x3.save2(x3_);
         if (x3_.isEmpty()) {
            if (Cm.rtype!=CR_CTR_ZERO) wblog(FL,
               "ERR %s() got rtype=%s",FCT,Cm.rtype.tostr());
            x3_.init(size_t(0),size_t(0),size_t(0)); 
         }

         if (long(x3_.SIZE[2])>Mc) wblog(FL,
            "ERR %s() x3 got larger OM (%s -> %dx%dx%d)\n"
            "    %s @ %s -> %s\n    %s @ %s -> %s\n"
            "--> %s", FCT, SSTR(x3_), Ma,Mb,Mc,
            A.toStr().data, STR(ica_), STR(ica),
            B.toStr().data, STR(icb_), STR(icb),
            STR(Cm) 
         );

         Cm.contract_x3(FL,A,ica_,B,icb_, cgp,xCGR&~3); 

         gStore.rclog(A.cgb->t, PFL, CG_VERBOSE>6 && F,
            "\b[u] XBUF[%03d] @ %s { #%05X #%05X #%05X }\n"
            "    %s @ %s -> %s\n    %s @ %s -> %s\n"
            "--> %s", XBUF.size(), Cm.x3Str().data, 
            Cm.a.cstat.cID, Cm.b.cstat.cID, Cm.c.cstat.cID,
            STR(A), STR(ica_), STR(ica),
            STR(B), STR(icb_), STR(icb), STR(Cm)
         );

        #ifndef WB_SKIP_ASSERT
         if (!CG_FIX) {
            unsigned i,j,k, d1=x3_.SIZE[0], d2=x3_.SIZE[1], d3=x3_.SIZE[2];
            double xd;

            for (k=0; k<d3; ++k)
            for (j=0; j<d2; ++j)
            for (i=0; i<d1; ++i) x3_(i,j,k) -= Cm.x3(i,j,k);

            if ((xd=double(x3_.norm()))>1e-12) wblog(FL,
            "ERR %s() got x3 data mismatch @ %.3g",FCT,xd);
         }
        #endif
      }
   }

   if (calc && !Cm.zflag) {
      if (Cm.c.cstat!=CD_IMPLICIT) { 
      gStore.save_XMap(FL,idc); }
   }

   if (Cm.rtype==CR_CTR_ZERO) {
      if (Cm.cgb) wblog(FL,"ERR %s() got %s",FCT,STR_(Cm.cgb));
      C.init(); C.rtype=Cm.rtype;
      return 0;
   }

   if (!mc) {
      if (Cm.x3.SIZE.len!=3) wblog(FL,
         "ERR %s() got compromised x3 %s",FCT,SSTR(Cm.x3)); 
      mc=Cm.x3.SIZE[2]; 
   }

   wbarray<double> x3, a_x3;
   const wbvector<size_t> &S=Cm.x3.SIZE;

   if (S.len!=3 ||
      ma!=int(A.wdim1()) || ma>int(S[0]) ||
      mb!=int(B.wdim1()) || mb>int(S[1]) || mc>int(S[2])) {

      MXPut(FL,"Idbg").add(Cm,"Cm")
        .add(ma,"ma").add(mb,"mb").add(mc,"mc")
        .add(Ma,"Ma").add(Mb,"Mb").add(Mc,"Mc");
      wblog(FL,
        "ERR %s() unexpected OM setting (%dx%dx%d / %dx%d @ %s)",
         FCT, ma,mb,mc, Ma,Mb, SSTR(Cm.x3)
      );
   }

   if (ma<int(S[0]) || mb<int(S[1])) {
      wbvector<size_t> S_(Cm.x3.SIZE); S_[0]=ma; S_[1]=mb;
      Cm.x3.resize(S_,x3);
   }
   else { x3.init(Cm.x3,'r'); }

   if (Cm.rtype==CR_CTR_SCALAR) {
      int e=(S[2]!=1 ? 1 : 0);
      if (C.cgw && (C.cgw.SIZE.len!=3 || C.cgw.SIZE[2]!=1)) { e|=2; }

      if (Cm.cgb) wblog(FL,"ERR %s() got %s",FCT,STR_(Cm.cgb));
      if (e) { wblog(FL,"ERR %s() "
         "unexpected OM setting (x3: %s / %s; e=%d)",SSTR(x3),SSTR(C.cgw),e);
      }
      C.init(); C.rtype=Cm.rtype;
   }

   if (ma==1 && mb==1) { 
      if (!mc || mc>long(x3.SIZE[2])) { wblog(FL,
         "ERR %s() got scalar x3=%s / OM=%dx%dx%d /%d\n"
         "  A: %-20s %-8s @ %s\n"
         "  B: %-20s %-8s @ %s",FCT,SSTR(Cm.x3),ma,mb,mc,Mc,
         STR_(A.cgb), SSTR(A.cgb->cgd), STR(ica),   
         STR_(B.cgb), SSTR(A.cgb->cgd), STR(icb));  
      }
      (C.cgw=x3)*=(A.wget1()*B.wget1()); 
   }
   else {
      wbperm p132("132");
      C.cgw.init(); 
      A.cgw.contract(FL,1,x3,1,a_x3)              
           .contract(FL,2,B.cgw,1,C.cgw,p132);    

      C.cgw.SkipTrailingZeroSpace(2,1,1e-14);  
   }

   if (Cm.rtype==CR_CTR_SCALAR) {
      return 0;
   }

   if (!Cm.cgb) wblog(FL,"ERR %s() got null cgb %s",FCT,STR(Cm));
   if (mc>(Mc=Cm.cgb->getOM()) || mc<=0) wblog(FL,
      "ERR %s() invalid mc=%d/%d (%s)",FCT,mc,Mc,SSTR(Cm.x3));

   if (C.cgb) { 
      if (C.cgb!=Cm.cgb) wblog(FL,
      "WRN %s() cgb changed\n   %s\n-> %s",FCT,STR_(C.cgb),STR_(Cm.cgb));
      if (C.cgp!=cgp) wblog(FL,
      "WRN %s() cgp changed (%s / %s)",FCT,STR(C.cgp),STR(cgp));
   }

   C.cgb=Cm.cgb;
   C.cgp=cgp;

   if (calc) { loaded=-loaded; }
   if (loaded>0) {
      if (Cm.rtype!=CR_DEFAULT) wblog(FL,
         "ERR %s() unexpected rtype=%s",FCT,Cm.rtype.tostr());
      C.rtype=Cm.rtype;
   }
   else if (!calc) {
      C.rtype=Cm.rtype;
   }
   else { C.rtype=CR_DEFAULT; }

   A.Reduce2Ref(FL); 
   B.Reduce2Ref(FL); 

   if (calc) { 
   gCS.reduceMemUsage(FL); } 

   return 1;
};

template <class TQ, class TD>
CRef<TQ>& X3Map<TQ,TD>::contractDegQ(const char *F, int L,
   const CRef<TQ> &A, 
   CRef<TQ> &B        
){
   if (B.cgb!=A.cgb || !B.cgw.sameSize(A.cgw)) wblog(FL, 
      "ERR %s() invalid CRefs (%p / %p)\ncgw: %s / %s",
      FCT, A.cgb, B.cgb, SSTR(A.cgw), SSTR(B.cgw));

   if (!A.cgb) { B.cgp.init(); return B; } 
   if (A.isAbelian()) {
      if (!A.cgb->isScalar() || A.cgw.numel()!=1 ||
         (A.cgb->cgd.D.len && A.cgb->cgd.D[0]!=1)) wblog(F_L,
         "ERR %s() unexpected scalar CRef %s",FCT,STR(A));
      return B;
   }

   unsigned r=A.rank(FL);

   if (r==2 && A.cgb->qdir.prod()<0) { 
      if (!A.isIdentityCG() || A.cgb->qdir.len!=2) wblog(FL,
         "ERR %s() unexpected rank-%d CGT %s",
         FCT,r,SSTR_(A.cgb));
      return B;
   }
   else if (r<2) {
      wblog(FL,"ERR %s() got rank-%d CRef\n%s",FCT,r,SSTR_(A.cgb));
   }

   unsigned ma=A.wdim2(), m=A.wdim1(), M=A.getOM();
   char loaded=0, calc=0; int isz=0, xOM=0;

   wbperm pa(A.cgp);
   pa.Permute(wbperm(B.cgp,'i')); 

   if (pa.isIdentityPerm()) wblog(FL,
      "ERR %s() %s @ p=%s",FCT,STR_(A.cgb),STR(pa));

   ctrIdx ica(pa,pa.conj+1);  
   ctrIdx icb, ica_, icb_;
   wbperm cgp; 

   icb.Index(r);
   ica.Sort_(icb); 

   A.adapt(ica_=ica,'i');
   B.adapt(icb_=icb,'i');

   cgc_contract_id<MTI> idc(*A.cgb,ica, *B.cgb,icb);

 # ifdef QS_USING_OMP
   CG::Guard xLK(FL,idc);
 # endif
   x3map<TQ,TD> &Cm=getXBUF(0,0,idc);

   cdata__ a; 

   if (!m || m>M) wblog(FL,
      "ERR %s() invalid OM data (%d/%d)",FCT,m,M);

   if (Cm.isEmpty()) {
      if (gStore.load_XMap(0,0,idc)>0) { loaded=1; }
   };

   if (Cm.isEmpty()) { calc|=1;
      xOM=((CDATA_TQ*)(B.cgb))->completeOM_DegQ(FL,pa); 
      isz=Cm.contract_x3(FL,A,ica_,B,icb_, cgp);

      if (Cm.cgb) wblog(FL,
         "ERR %s() got cgb data\n%s",FCT,STR_(Cm.cgb));
      if (Cm.c.qs.len || Cm.c.qdir.len) wblog(FL,
         "ERR %s() got non-empty c\n%s",FCT,STR(Cm.c));
      if (Cm.x3.SIZE.len!=3 || Cm.x3.SIZE[2]!=1) wblog(FL,"ERR %s() "
         "got unexpected x3 data (%s)",FCT,SSTR(Cm.x3));

      gStore.rclog(A.cgb->t, PFL,
         CG_VERBOSE>8 || (CG_VERBOSE>6 && Cm.x3.numel()>1),
         "\b[+] XBUF[%03d] @ %s { #%05X #%05X }\n"
         "    %s\n  @ %s * %s -> %s\n--> %s", XBUF.size(), SSTR(Cm.x3),
         Cm.a.cstat.cID, Cm.b.cstat.cID,
         STR(A), STR(ica), STR(icb), STR(pa), STR(Cm)
      );
   }
   else {
      if (Cm.cgb || Cm.x3.SIZE.len!=3 ||
          Cm.x3.SIZE[0]!=Cm.x3.SIZE[1] || Cm.x3.SIZE[2]!=1) {
         wblog(FL,"ERR %s() invalid x3 (%s; l=%d)",
         FCT,SSTR(Cm.x3),loaded);
      }

      if (M!=Cm.x3.SIZE.data[0]) { calc|=2;
         wbarray<TD> x3_; Cm.x3.save2(x3_);
         if (!A.isRefInit() && M < x3_.SIZE.data[0]) wblog(FL,
            "ERR %s() x3map() size inconsistency\n(%dx%dx1 <> %s)\n"
            "having %s",FCT,M,M,SSTR(x3_),STR_(A.cgb)
         );

         xOM=((CDATA_TQ*)(A.cgb))->completeOM_DegQ(FL,pa); 
         isz=Cm.contract_x3(FL,A,ica_,B,icb_, cgp);

         gStore.rclog(A.cgb->t, PFL, CG_VERBOSE>6 && F,
            "\b[u] XBUF[%03d] @ %s { #%05X #%05X }\n"
            "    %s\n  @ %s * %s -> %s\n--> %s", XBUF.size(), SSTR(x3_),
            Cm.a.cstat.cID, Cm.b.cstat.cID,
            STR(A), STR(ica), STR(icb), STR(pa), STR(Cm)
         );

      #ifndef WB_SKIP_ASSERT
         if (!x3_.fitsSize(Cm.x3)) {
            wblog(FL,"ERR %s() %s -> %s",FCT,SSTR(x3_),SSTR(Cm.x3)); }

         unsigned i,j,k, d1=x3_.SIZE[0], d2=x3_.SIZE[1], d3=x3_.SIZE[2];
         double e3=0;

         for (k=0; k<d3; ++k)
         for (j=0; j<d2; ++j)
         for (i=0; i<d1; ++i) e3+=Wb::norm2(double(x3_(i,j,k) - Cm.x3(i,j,k)));

         if ((e3=sqrt(e3))>1e-12) {
            wblog(FL,"ERR %s() got x3 data mismatch @ %.3g",FCT,e3);
         }
      #endif
      }
   }

   if (calc) { 
      if (Cm.zflag) gStore.rclog(A.cgb->t,PFL,F, 
         "WRN %s() got zero-contraction (%d/%d)",FCT,xOM,calc);
      gStore.save_XMap(FL,idc);
   }

   B.cgw.init(M,ma); 

   if (Cm.rtype==CR_CTR_SCALAR) { 
      const wbvector<SPIDX_T> &S=Cm.x3.SIZE;

      if (Cm.cgb) wblog(FL,"ERR %s() got %s",FCT,STR_(Cm.cgb));
      if (!m || m>M) wblog(FL,"ERR %s() invalid OM data (%d/%d)",FCT,m,M);
      if (S.len!=3 || S[0]!=S[1] || S[2]!=1 || M>S[0]) {
         wblog(FL,"ERR %s() "
           "OM out of bounds (%dx%dx1; %s; c=%d)\nA: %s\nB: %s",
           FCT,M,M,SSTR(Cm.x3),calc,STR(A),STR(B));
      }
      if (M<S[0]) { M=S[0]; B.cgw.init(M,ma); }

      if (S.prod()==1) { 
         B.cgw[0] = Cm.x3[0] * A.cgw[0];  
      }
      else { long M_=M, ma_=ma, m_=m; double afac=1, cfac=0;
         dgemm("T","N",&M_,&ma_,&m_,
            &afac, Cm.x3.data, &M_, 
                   A.cgw.data, &m_, 
            &cfac, B.cgw.data, &M_  
         );
      }

      double e=-1, x[2]; x[0]=x[1]=0;
      if (!A.cgw.isOrthoCols(x,0,CG_SKIP_DEPS1)) wblog(FL,
         "ERR %s() got non-orthogonal A.cgw",FCT);
      if (x[0]<CG_SKIP_DEPS1) wblog(FL,
         "ERR %s() got small A.cgw (%.3g)",FCT,x[0]);

      if (B.cgw.isOrthoCols(x+1,0,CG_SKIP_DEPS1)) {
         e=Wb::abs(x[1]-x[0])/MAX3(double(1),x[0],x[1]);

         if (e<CG_SKIP_DEPS1) {
            if (e>CG_SKIP_DEPS2) wblog(FL,"WRN %s() "
               "got cgw difference %.3g (%.3g / %.3g)",e,FCT,x[0],x[1]);
            return B;
         }
      }

      MXPut(FL,"Ix")
        .add(QSet<TQ>(A),"A").add(A.cgw,"a")
        .add(QSet<TQ>(B),"B").add(B.cgw,"b")
        .add(Cm,"Cm").add(wbvector<double>(2,x),"x");
      wblog(FL,"ERR %s() check new OM component (e=%g)",FCT,e);
   }
   else {
      wblog(FL,"ERR %s() check new OM component (zero-ctr)",FCT);
   }

   gStore.rclog(B.cgb->t, PFL,
      (CG_VERBOSE>3 && isz) 
   || (CG_VERBOSE>6 && F),"[+] CBUF[%d] new/dec #%05X %s %+d",
      gCS.BUF.size(), B.cgb->cstat.cID, STR_(B.cgb),xOM);

   if (B.cgb->qdir.len==3) wblog(FL,
      "WRN %s() decomposing rank-3 CData\n%s",FCT,STR_(B.cgb));

   return B;
};

template <class TQ, class TD>
void X3Map<TQ,TD>::Info(const char *F, int L) const {

 # ifdef QS_USING_OMP
   Wb::ompGuard xLK(XS_buf); 
 # endif

   if (XBUF.size()) {
      size_t mtot=0;
      for (auto it=XBUF.begin(); it!=XBUF.end(); ++it) {
         mtot+=(it->second.memSize() + it->second.memSize());
      }

      wblog(F_L," *  XStore::XBUF got %d entries @ %s",
      XBUF.size(), Wb::size2Str(mtot).data);
   }
   else { wblog(F_L," *  XStore::XBUF is empty"); }
};

#endif

