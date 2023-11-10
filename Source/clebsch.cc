/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : clebsch (for abelian and non-abelian symmetries)
 *
 * Copyright 2022 Andreas Weichselbaum
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

#ifndef __WB_CLEBSCH_GORDAN_CC__
#define __WB_CLEBSCH_GORDAN_CC__

/* major changes in cgw->mat update
   >> isIdentityCG() nrm: double -> wbvector<double> // Wb,Jun08,20
   >> getTensorProdReps_gen() storing cgw as matrix in map3 // Wb,Jun11,20
*/

int get_CG_VERBOSE(const char *F, int L) {

   int k=0,i;

   if ((i=Wb::GetEnv(0,0,"CG_VERBOSE",k))==0) {
      if (k>9) wblog(F_L,
         "WRN %s() invalid CG_VERBOSE=%d (0..9)",FCT,k);
      else {
         if (F && k && k!=CG_VERBOSE) {
            wblog(F,L,"ENV %s() CG_VERBOSE = (%d -> %d)",
               myname, CG_VERBOSE, k);
            CG_VERBOSE=k;
         }
      }
   }
   else if (int(i)!=-1) wblog(F_L,
      "WRN %s() invalid CG_VERBOSE (e=%d)",FCT,i);

   Wb::GetEnv(0,0,"RC_LOG",RC_LOG); 
   return k;
};

unsigned get_CG_FIXIT(const char *F, int L) {

   unsigned k=0,i;

   if ((i=Wb::GetEnv(0,0,"CG_FIXIT",k))==0) {
       if (F && k && (k!=CG_FIXIT || CG_VERBOSE>6)) {
          wblog(F,L,"ENV %s() CG_FIXIT = (%d -> %d)",
             myname, CG_FIXIT, k);
          CG_FIXIT=k;
       }
       return k;
   }
   else if (int(i)!=-1) wblog(F_L,
      "WRN %s() invalid CG_FIXIT (e=%d)",FCT,i);
   return k;
};

unsigned get_CG_PREVIEW(const char *F, int L) {

   unsigned k=0,i;

   if ((i=Wb::GetEnv(0,0,"CG_PREVIEW",k))==0) {
       if (F && k && (k!=CG_PREVIEW || CG_VERBOSE>6)) {
          wblog(F,L,"ENV %s() CG_PREVIEW = (%d -> %d)",
             myname, CG_PREVIEW,k);
          CG_PREVIEW=k;
       }
       return k;
   }
   else if (int(i)!=-1) wblog(F_L,
      "WRN %s() invalid CG_PREVIEW (e=%d)",FCT,i);
   return k;
};

inline bool QType::validType(const char *F, int L) const {

   if (type<=0 || type==QTYPE_ASEP) {
      if (F) wblog(F,L,
         "ERR %s() unspecified symmetry (%d/%d; %s)",
         FCT, type, QTYPE_NUM_TYPES, STR_(this));
      return 0;
   }
   if (type>=QTYPE_NUM_TYPES) {
      if (F) wblog(F,L,
         "ERR %s() type out of bounds (%d/%d)",FCT,type,QTYPE_NUM_TYPES);
      return 0;
   }

   if (type<QTYPE_ASEP) {
      if (type==QTYPE_ZN) { if (sub>=2) return 1;
         if (F) wblog(F,L,
            "ERR %s() invalid sub=%d for %s [>=2]",FCT,sub,QTYPE_STR[type]);
         return 0;
      }
      if (!sub) return 1;
      if (F) wblog(F,L,
         "ERR %s() invalid sub=%d for %s [0]",FCT,sub,QTYPE_STR[type]);
      return 0;
   }

   if (type>QTYPE_ASEP && sub>=1) return 1;

   if (F) wblog(F,L, 
      "ERR %s() invalid sub=%d for %s [>=2]",FCT,sub,QTYPE_STR[type]);

   return 0;
};

int QType::init_s(const char *s) {

   int k=1, e=0;

   if (!s || !s[0]) { type=QTYPE_UNKNOWN; sub=0; return 1; }

   for (; k<QTYPE_NUM_TYPES; ++k) {
      if (!strcmp(s,QTYPE_STR[k])) {
         type=QTYPE_SET(k); sub=0; 
         validType(FL); return 0;
      }
   }

   if (s[0]=='Z') {
      if ((e=this->atoi(s+1,k,3))>0 && k>1) {
         type=QTYPE_ZN; sub=k; 
      }
      else return (-60+e);
   }
   else if (s[0]=='S') {
      if (s[1]=='U') { 
         if ((e=this->atoi(s+2,k,2))>0 && k>1) {
            type=QTYPE_SUN; sub=k-1;
         }
         else return (-70+e);
      }
      else if (s[1]=='p') { 
         if ((e=this->atoi(s+2,k,2))>0 && k>=4 && !(k%2)) {
            type=QTYPE_SpN; sub=k/2;
         }
         else return (-80+e);
      }
      else if (s[1]=='O') { 
         if ((e=this->atoi(s+2,k,2))>0 && k>=3) {
            if (k%2)
                 { type=QTYPE_SON; sub=(k-1)/2; }
            else { type=QTYPE_SEN; sub=(k  )/2; }
         }
         else return (-80+e);
      }
      else return -90;
   }
   else return -90;

   validType(FL);
   return 0;
};

template <class TQ> inline
unsigned QType::qdim(const TQ* q) const {

   if (isAbelian()) { return 1; }

   unsigned d=0;
   QSet<TQ> Q; { Q.init2(*this,q); }

   CDATA_TQ *C = gCS.BUF_find(Q);
   if (C && !C->isEmpty()) { d=C->cgd.cgsparray::dim(); }

   if (!d) { 
      genRG_base<TQ,RTD> *R=gRS.find_RSet(*this,q);
      if (R && R->Z.dim2) d=R->Z.dim1;

   if (!d) { 
      C=&gCS.getBUF(0,0,Q,LB_LOAD); 
      if (C->cgd.SIZE.len) {
         if (C->cgd.SIZE.len!=2 || C->cgd.SIZE[0]!=C->cgd.SIZE[1]) wblog(FL,
            "ERR %s() got unexpected scalar CData\n%s\nSIZE=[%s]/%ld",
            FCT,C->toStr().data, C->cgd.sizeStr().data, C->cgd.D.len
         );
      }
      d=C->cgd.D.len; 

      if (!d) { 
         const genRG_base<TQ,RTD> &R=gRS.getR(FL,*this,q);
         if (R.Z.dim2) d=R.Z.dim1;
      }

   if (!d) {
      if (isSU2()) { 
         int n=q[0]; 
         if (n<0 || double(n)!=q[0]) wblog(FL, 
            "ERR invalid %s symmetry q=%g",STR_(this),double(q[0]));
         d=(n+1);
      }
      else wblog(FL,
        "ERR %s() failed to determine multiplet dimension for\n%s",
         FCT,STR(Q)
      );
   }}}

#  ifndef WB_SKIP_ASSERT
   if (isSU2()) { 
      if (d!=unsigned(q[0]+1)) wblog(FL,
         "ERR %s() unexpected multiplet dimension (%d/%d)\nhaving %s",
         FCT,d,q[0]+1,STR(Q)
      );
   }
#  endif

   return d;
};

template <class TQ> inline
wbvector<unsigned>& QType::QDim(
   const TQ* q, unsigned r, unsigned stride, wbvector<unsigned> &S
 ) const {

   S.init(r);

   if (isAbelian()) { S.set(1); } else
   for (unsigned i=0; i<r; ++i, q+=stride) { S[i]=this->qdim(q); }

   return S;
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
      if (x<=0 || fabs(x-round(x))>1E-6)
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
      if (x<=0 || fabs(x-round(x))>1E-6)
         wblog(FL,"ERR %s() got d = %d/%d = %g",FCT,P_,Q_,x);
   }

   return size_t(P_/Q_);
};

namespace DY {

Symmetry& Symmetry::init(const QType &q) {

   unsigned i=1;

   if (q.type==QTYPE_SUN) { 
      if (q.sub<1) wblog(FL,"ERR %s() got symmetry %s",FCT,STR(q));
   } else
   if (q.type==QTYPE_SpN) { 
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

   if (q.type==QTYPE_SUN) {
      unsigned i2,j, k=0, l=0;
      M.init(n,(n*(n+1))/2); 
      for (k=1; k<=n; ++k) { 
         for (i2=n-k, i=0; i<=i2; ++i, ++l) {
            for (j=0; j<k; ++j) { M(i+j,l)=1; }
         }
      }
   }
   else if (q.type==QTYPE_SpN) {
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

   if (!qm.len || qm.len!=q.qlen()) wblog(FL,
      "ERR %s() invalid qset %s (%s)",FCT,STR(q),STR2(qm,q));
   if (qm.len!=A.dim2 || A.dim1!=A.dim2) wblog(FL,
      "ERR %s() invalid Cartan matrix A (%dx%d / %d; %s)",
      FCT,A.dim1,A.dim2,qm.len,STR(q));
   if (M.SIZE.len!=2 || M.SIZE[0]!=n || M.SIZE[1]<n) wblog(FL,
      "ERR %s() invalid M data (%s /%d)",FCT,SSTR(M),n);

   qset<TQ> qs(n); np=M.SIZE[1];
   R.init(qm);

   while (++iter<niter) {
      for (auto Il=R.X.begin(); Il!=R.X.end(); ++Il) { if (Il->second) {
         weight_info<TQ> &w0=R.W[Il->first];

         if (w0.p.len!=n || !w0.m) wblog(FL,
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
      wblog(FL,"ERR %s() got iter=%d",FCT,iter);
   }
   else if (F) wblog(FL," *  %s() got iter=%d",FCT,iter);

   map<int,   
      map<qset<TQ>,   
         weight_info<TQ>*  
      >
   > Ml;

   for (auto I=R.W.begin(); I!=R.W.end(); ++I) {
      l=I->second.p.sum();
      if (I->first.anyLT(0)) wblog(FL,"ERR %s() got w=[%s] "
         "outside dominant Weyl chamber",FCT,STR(I->first));
      if (!I->second.p.len || !I->second.m) wblog(FL,
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

   if (F) wblog(FL,"==> %s[%s]: q=(%s) ",FCT,STR(q),STR2(qm,q));
                  if (F) R.print(FL,q); 
   wExpand(FL,R); if (F) R.print(FL,q);

   d=R.dim();
   i=q.wdim(qm.data); 
   if (d!=i) wblog(FL,
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

   unsigned i,k, ip=0, np=M.SIZE[1];
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
      for (k=0;; ++k) {
         qk.Plus(ap.data);
         if (!(mk=Get(qk,R.W,vflag))) break;
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

   if (!n || n!=A.dim1 || n!=A.dim2) wblog(FL,
      "ERR %s() invalid Cartan matrix %s /%d",FCT,SSTR(A),n);

   for (auto Iw=R.W.begin(); Iw!=R.W.end(); ++Iw) { if (Iw->second.m>0) {
      W0[Iw->first]=&(Iw->second);
   }}

   for (auto Iw=W0.begin(); Iw!=W0.end(); ++Iw) {
      X.clear(); X[Iw->first]=Iw->second; m0=Iw->second->m; iter=0;
      while (++iter<niter) {
         for (auto Ik=X.begin(); Ik!=X.end(); ++Ik) { if (Ik->second) {
            if (Ik->second->m<=0) wblog(FL,
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
      if (iter>=niter) wblog(FL,
         "ERR %s() got iter=%d/%d",FCT,iter,niter);
   }
   return nX;
};

template <class TQ>
void Weights<TQ>::print(
   const char *F, int L, const QType &q, char vflag) const {

   size_t n=W.size(); int l, isn;
   char fmt[64];
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
      snprintf(s,128,"total of %d + %d = %d weight%s (%d state%s)",
         ndom, n-ndom, n, n!=1 ? "s":"", mtot, mtot!=1 ? "s":"",
         merr ? ", incomplete!":"");
      if (!merr) wblog(FL," *  %s",s); else wblog(FL,"WRN %s",s);
   }
   else printf("\n");

   if (n) { n=3*W.begin()->first.len; 
      sprintf(fmt,"\n  level  %%-%ds %%-%ds     m    # weights",n+3,n+3);
      printf(fmt,"dcoeffs","pcoeffs"); 
      sprintf(fmt,"\n  %%5d%%1s (%%%ds ) [%%%ds ] %%5d",n,n);
   }
   else fmt[0]=0;

   for (auto Il=M.begin(); Il!=M.end(); ++Il) {
      auto &M2=Il->second; n=M2.size();
      unsigned i=0;

      for (auto I2=M2.begin(); I2!=M2.end(); ++I2) {
         const weight_info<TQ> &w=*I2->second;

         printf(fmt, Il->first,
            X.find(I2->first)==X.end() ? "" : "*",
            I2->first.wbvector<TQ>::toStr(2).data,
            w.p.wbvector<TQ>::toStr(2).data, w.m
         );

         if (++i==1) printf("   %4d",n);
      }
   }
   printf("\n\n");
};

}; 

QVec& QVec::init(
   const char *F, int L, const char *s,
   unsigned D0 
){

   if (!s || !s[0] || !strcmp(s,CGC_ALL_ABELIAN)) {
      init(); return *this;
   }

   unsigned l,i, n=strlen(s), i1=0, i2=0;
   unsigned D=(D0 ? 2*D0 : 8); 
   char sep=0, sx[n+1];
   QType q;

   init(D);

   for (l=i=0; i<=n; i++) { 
      if (s[i]==' ' || s[i]=='\t' || s[i]=='\n' ||
          s[i]==',' || s[i]==';') { sep++; continue; }

      if (sep || s[i]==0) {
         if (l>=D) { l++; break; } 

         sx[i2]=0; q.init(F,L,sx+i1);
         data[l++]=q; i1=i2; sep=0;
      }
      sx[i2++]=s[i];
   }

   if (l>D) {
      wblog(F,L,"ERR invalid qtype specs: %s\n"
     "(string contains too many entries; %d/%d)",s,l,D);
   }

   if (l) len=l; else init();

   if ((D0 && D0!=l) || l>Qlen()) wblog(F,L, 
   "ERR qtype length inconsistency (%d/%d; %d)",D0,l,Qlen());

   return *this;
};

unsigned QVec::Qlen(unsigned n) const {

   unsigned l=0, i=0;
   if (int(n)<0) { n=len; } else
   if (n>len) wblog(FL,"ERR %s() size out of bounds (%d/%d)",FCT,n,len);

   for (; i<n; ++i) { l+=data[i].qlen(); }
   return l;
};

unsigned QVec::Qlen(wbvector<unsigned> &dd) const {
   unsigned i=0, n=0;
   dd.init(len); for (; i<len; ++i) { n+=(dd[i]=data[i].qlen()); }
   return n;
};

unsigned QVec::Qlen(wbvector<unsigned> &dd, wbvector<unsigned> &dz) const {
   unsigned i=0, n=0;
   dd.init(len); for (   ; i<len; ++i) { n+=(dd[i]=data[i].qlen()); }
   dz.init(len); for (i=0; i<len; ++i) { dz[i]=data[i].qrank(); }
   return n;
};

unsigned QVec::Qpos(wbvector<unsigned> &dc) const {

   dc.init(len); if (len) { unsigned i=1;
      for (; i<len; ++i) {
         dc[i] = dc[i-1] + data[i-1].qlen();
      }
      return (dc[i-1] + data[i-1].qlen());
   }
   return 0;
};

unsigned QVec::Qrank() const {
   unsigned n=0, i=0;
   for (; i<len; ++i) { n+=data[i].qrank(); }
   return n;
};

template <class TQ>
wbvector<unsigned>& QVec::QDim(
   const TQ *q, wbvector<unsigned> &S
 ) const {

   unsigned d, i=0; S.init(len);

   for (; i<len; ++i, q+=d) {
      d=data[i].qlen();
      S[i]=data[i].qdim(q);
   }

   return S;
};

template <class TQ>
wbvector<unsigned>& QVec::QDim(
   const TQ *q, unsigned k, unsigned r, wbvector<unsigned> &S
 ) const {

   unsigned i, d0, n=0;

   if (k>=len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,len);

   if (!r) { S.init(); return S; }
   if (data[k].isAbelian()) { S.init(1).set(1); return S; }

   for (i=0; i<len; ++i, n+=d0) {
      d0=data[i].qlen();
      if (i==k) q+=n;
   }

   if (S.len!=r) S.init(r);
   for (i=0; i<r; i++, q+=n) { S[i]=data[k].QDim(q); }

   return S;
};

template <class TQ>
unsigned QVec::QDim(const TQ *q) const {

   unsigned i,d0, n=(len ? 1 : 0);

   for (i=0; i<len; ++i) {
      d0=data[i].qlen();
      n*=data[i].qdim(q); q+=d0;
   }

   return n;
};

template <class TQ>
unsigned QVec::QDim(const TQ *q, unsigned k) const {

   for (unsigned d0=0, i=0; i<len; ++i, q+=d0) {
      d0=data[i].qlen();
      if (i==k) return data[i].QDim(q);
   }

   wblog(FL,"ERR index out of bounds (%d/%d)",k,len);
   return 0;
};

template <class TQ>
wbMatrix<TQ>& QVec::getQsub(
   const wbMatrix<TQ>& Q0, const wbindex &I,
   wbMatrix<TQ>& QI,
   wbMatrix<TQ> *Qx  
 ) const {

   if (Q0.isEmpty()) { QI.init(); if (Qx) Qx->init(); return QI; }

   unsigned i,j,k,l,n;
   wbvector<unsigned> s0,s,ic0,ic;
   const TQ *d0=Q0.data; TQ *d;

   Qlen(s0); n=s0.sum();

   if (Q0.dim2!=n) wblog(FL,
      "ERR %s() Q-set size mismatch (%dx%d/%d)",FCT,Q0.dim1,Q0.dim2,n);

   s0.cumsum_(ic0); ic0.select(I,ic);
   s0.get(I,s); QI.init(Q0.dim1,s.sum()); d=QI.data;

   for (i=0; i<Q0.dim1; ++i, d0+=Q0.dim2)
   for (j=0; j<I.len; ++j, d+=n) { 
      for (n=s[j], k=ic[j], l=0; l<n; ++l) d[l]=d0[k+l];
   }

   if (Qx) {
      wbindex Ix; I.invert(len,Ix);
      s0.get(Ix,s); Qx->init(Q0.dim1,s.sum()); d=Qx->data;
      ic0.select(Ix,ic); d0=Q0.data;

      for (i=0; i<Q0.dim1; ++i, d0+=Q0.dim2)
      for (j=0; j<Ix.len; ++j, d+=n) { 
         for (n=s[j], k=ic[j], l=0; l<n; ++l) d[l]=d0[k+l];
      }
   }

   return QI;
};

wbstring QVec::toStr(const char vflag) const {

   wbstring s; 

   if (!len) {
      if (vflag!='V')
           s=(vflag ? CGC_ALL_ABELIAN : "");
      else s="(all abelian U1)";
      return s;
   }

   char tflag=(vflag ? 0:1);
   s.init(8*len+1); s.data[0]=0;

   for (unsigned k=0, i=0; i<len; i++) {
      if (i) { s.cat(k,tflag? ",":"*"); }
      s.cat(k, data[i].toStr(tflag).data);
   }
   return s;
};

int QVec::print_qset(
   const char *F, int L, const gTQ *qs, char *s, unsigned n) const {

   unsigned l=0; 
   unsigned k=0, m=0, i;

   if (!qs) { wblog(F_L,"ERR %s() got null qs",FCT); }
   if (! s) { wblog(F_L,"ERR %s() got null string buffer (n=%d)",FCT,n); }
   if (int(n)<=0 || n<2*len) { wblog(F_L,
      "ERR %s() buffer size too small (n=%d @ nsym=%d)",FCT,n,len); }
   if (!len) { s[l]=0; return l; }

   if (!WbUtil<gTQ>::isInt()) wblog(F_L,
      "ERR %s() got non-integer type for qset  '%s'",FCT, TSTR(gTQ));
   if (!len && qs) wblog(FL, 
      "WRN %s() got empty qtype (q=%d,..)",FCT,*qs);

   for (; k<len && l<n; ++k, qs+=m) { if (k) { s[l++]=' '; }
      if (data[k].isAbelian()) { m=1;
         if (l<n) { l+=snprintf(s+l,n-l,"%2d",qs[0]); }
      }
      else {
         m=data[k].qlen();
         if (l+m<n) {
            if (CG::gotQAlpha(qs,m)) {
               if (CG::qset2cstr(qs,s+l,m)) wblog(F_L, 
                  "ERR %s() failed to obtain compact qset string",FCT);
               l+=m;
            }
            else {
               for (i=0; i<m && l<n; ++i) {
               l+=snprintf(s+l,n-l," %d",qs[i]); }
            }
         }
         else {
            if (l<n) { snprintf(s+l,n-l,".."); }
            l+=m;
         }
      }
   }

   if (l<n) { s[l]=0; }
   else if (F) { s[n-1]=0; wblog(F,L,
      "ERR %s() string out of bounds (%s, l=%d/%d",FCT,s,l,n); }

   return l;
};

QDir& QDir::init(const iTags& it) {

   init2val(it.len,+1);
   for (unsigned i=0; i<len; ++i) {
      if (it[i].isConj()) data[i]=-1;
   }
   return *this;
};

QDir& QDir::init(const char *F, int L,const char *s) {

   if (!s || !s[0]) { init(); }
   else {
      unsigned i=0, l=(s? strlen(s): 0);
      if (l>32) wblog(F_L,
         "ERR %s() unexpected QDir string '%s' (len=%d)",FCT,s?s:"",l);
      for (init(l); i<l; ++i) {
         if (s[i]=='+') { data[i]=+1; } else
         if (s[i]=='-') { data[i]=-1; } else
         if (s[i]=='0') { data[i]= 0; } else
         wblog(F_L,"ERR QDir::%s() invalid string `%s'",FCT,s);
      }
   }
   return *this;
};

QDir& QDir::init_iout(const char *F, int L,
   unsigned r, unsigned iout 
){
   if (!r) {
      if (iout) wblog(F_L,
         "ERR %s() invalid iout=%d (having r=%d)",FCT,iout,r);
      return init();
   }
   if (iout) {
      if ((--iout)>=r) wblog(F_L,
         "ERR %s() index out of bounds (%d/%d)",FCT,iout+1,r);
      init2val(r,+1); data[iout]=-1;
   }
   else { init(r); }

   return *this;
};

bool QDir::operator==(const char *s) const {

   if (!s) {
      if (len) wblog(FL,"ERR %s() got null string",FCT);
      return 1;
   }
   if (len!=strlen(s)) { return 0; }

   for (unsigned i=0; i<len; ++i) {
      if (s[i]=='+') { if (data[i]<=0) { return 0; }} else
      if (s[i]=='-') { if (data[i]>=0) { return 0; }}
      else wblog(FL,"ERR got qdir=='%s' (%d/%d)",s,i+1,len);
   }

   return 1;
};

wbstring QDir::toStr() const {
   wbstring s(16); unsigned i=0, l=0, n=s.len;
   for (; i<len && l<len; ++i) {
      if (data[i]==+1) { s[i]='+'; ++l; } else
      if (data[i]==-1) { s[i]='-'; ++l; } else
      if (data[i]== 0) { s[i]='0'; ++l; } else
      l+=snprintf(s.data+l,n-l,"<%+d>",data[i]);
   }
   if (l>=n) { s.data[n-1]=0; wblog(FL,"ERR %s() "
      "string out of bounds (%d/%d)\n`%s...'",FCT,l,n,s.data);
   }
   return s;
};

wbstring QDir::toTag() const { wbstring s(len+1); 
   for (unsigned i=0; i<len; ++i) {
      if (data[i]==+1) s[i]='p'; else
      if (data[i]==-1) s[i]='m'; else
      wblog(FL,"ERR %s() invalid qdir (%d: %d)",FCT,i,data[i]);
   }; s[len]=0;

   return s;
};

template <class TQ>
qset<TQ>::qset(const qset<TQ> &q1, char op, const qset<TQ> &q2) {

   if (q1.len!=q2.len) wblog(FL,"ERR size mismatch (%d/%d)",q1.len,q2.len);
   this->data=NULL;

   switch (op) {

      case '+': 
         INIT(q1.len,q1.data);
         for (unsigned i=0; i<q1.len; i++) q1.data[i]+=q2.data[i];
         break;

      case '-': 
         INIT(q1.len,q1.data);
         for (unsigned i=0; i<q1.len; i++) q1.data[i]-=q2.data[i];
         break;

      case '.': 
         INIT(q1.len,q1.data,q2.len,q2.data);
         break;

      default: wblog(FL,"ERR invalid operator=%c<%d> ???",op,op);
   }
};

template <class TQ>
qset<TQ>& qset<TQ>::init(
   unsigned l1, const TQ *d1,
   unsigned l2, const TQ *d2
){
   unsigned i=0;
   wbvector<TQ>::RENEW(l1+l2);
   TQ *d=(this->data);

   if (l1!=l2) wblog(FL,
      "WRN cat() qlabel length mismatch (%d+%d)",l1,l2);
   if (l1 && (!d1 || !d2)) wblog(FL,
      "ERR cat(%d+%d) got NULL data (0x%Xl, 0x%Xl)",l1,l2,d1,d2);

   for (   ; i<l1; ++i) { d[i]=d1[i]; }; d+=l1;
   for (i=0; i<l2; ++i) { d[i]=d2[i]; }

   return *this;
};

template <class TQ>
qset<TQ>& qset<TQ>::init(
   unsigned l1, const TQ* d1,
   unsigned l2, const TQ *d2,
   unsigned l3, const TQ *d3
){
   unsigned i=0;
   wbvector<TQ>::RENEW(l1+l2+l3);
   TQ *d=(this->data);

   if (l1!=l2 || l1!=l3) wblog(FL,
      "WRN cat() qlabel length mismatch (%d+%d+%d)",l1,l2,l3);
   if (l1 && (!d1 || !d2 || !d3)) wblog(FL,"ERR cat(%d+%d+%d) "
      "got NULL data (0x%Xl, 0x%Xl, 0x%Xl)",l1,l2,l3,d1,d2,d3);

   for (   ; i<l1; ++i) { d[i]=d1[i]; }; d+=l1;
   for (i=0; i<l2; ++i) { d[i]=d2[i]; }; d+=l2;
   for (i=0; i<l3; ++i) { d[i]=d3[i]; }

   return *this;
};

template <class TQ>
qset<TQ>& qset<TQ>::init(const char *F, int L, const char *s, QType *t) {

   unsigned n=(s ? strlen(s) : 0);
   if (n>32) wblog(F_L, 
      "ERR %s() qset out of bounds (len=%d)\n%s",FCT,n,s?s:"");

   unsigned i=0, j=0, l=0, bflag=0;
   wbvector<TQ>::init(n);

   for (; i<n; ++i) {
      if (s[i]!=' ') {
      if (s[i]=='(') { if (!bflag) { ++bflag; } else break; }
      else break; }
   }

   l=CG::cstr2qset(s+i,this->data+j); 
   i+=l; j=+l;

   for (; i<n; ++i) {
      if (s[i]!=' ') { 
         if (s[i]==')' && bflag==1) { ++bflag; }
         else break;
      }
   }

   l=(t? t->qrank() : 0);
   if (j && i==n && (!l || l==j)) {
      this->len=j; return *this;
   }

   for (i=0; i<n; ++i) { if (isalpha(s[i])) wblog(FL,
      "ERR %s() got unexpected string `%s'",FCT,s);
   }

   Str2Idx(F_L,s,(wbvector<TQ>&)(*this),0,0); 

   if (!this->len || (l && this->len!=l)) wblog(FL,"ERR %s() "
      "got unexpected qset `%s' (len=%d/%d)",FCT,s,this->len,l);
   return *this;
};

template <class TQ>
int qset<TQ>::checkNQs(const QVec &qq, const char *F, int L) {

   unsigned i,d0, l=0;

   for (i=0; i<qq.len && l<this->len; i++) {
      d0=qq[i].qlen(); l+=d0;
   }

   if (i!=qq.len || l!=this->len) {
      if (F) wblog(F,L,"ERR qset/QType inconsistency (%d/%d, %d/%d)",
      i,qq.len,l,this->len); return -1;
   }

   return i;
};

template <class TQ>
bool qset<TQ>::isConsistent(
   const QVec &qq, int rank, const char *F, int L
) const {

   unsigned i,d0, l=0, m=(rank>0 ? rank : 1)*qq.len;
   if (qq.isEmpty()) wblog(F_L,"ERR got empty QVec");

   for (i=0; l<(this->len) && i<m; i++) {
      d0=qq[i<qq.len ? i : i%qq.len].qlen();
      l+=d0;
   }

   if (i>qq.len && (rank<=0 || i%qq.len)) {
      if (F) wblog(F,L,
      "ERR incompatible QVec (%d/%d; %d)",i,qq.len,rank);
      return 0;
   }

   if (l!=this->len || i!=(rank>0 ? rank : 1)*qq.len) {
      if (F) wblog(F,L,"ERR qset out of bounds "
         "(%d/%d, %d/%d; %d)",l,this->len,i,qq.len,rank);
      return 0;
   }

   return 1;
};

template <class TQ>
unsigned qset<TQ>::QDim(const QVec &qq, int rank) const {

   unsigned i,d0, l=0, n=1;
   if (qq.isEmpty()) wblog(FL,"ERR got empty QVec");

   for (i=0; l<this->len && i<this->len; i++) {
      if (rank<=0) { if (i>=qq.len) wblog(FL,
         "ERR QVec index out of bounds (%d/%d; %d)",i,qq.len,rank);
         d0=qq[i].qlen();
      }
      else d0=qq[i%qq.len].qlen();

      n*=qq[i].QDim(this->data+l);
      l+=d0;
   }

   if (i>qq.len) {
      if (rank<=0 || i%qq.len) 
      wblog(FL,"ERR incompatible QVec (%d/%d; %d)",i,qq.len,rank);
   }

   if (l!=this->len || i>l || i%qq.len) wblog(FL, 
   "ERR qset/QType mismatch (%d/%d, %d/%d)",i,qq.len,l,this->len);

   return n;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init(const char *F, int L,
   const CRef<TQ> &A_, const ctrIdx &ica,
   const CRef<TQ> &B_, const ctrIdx &icb,
   wbperm *Pcgd
){

   if (!A_.cgb || !B_.cgb) {
      if (!A_.isAbelian() || !B_.isAbelian()) wblog(F_L,
         "ERR %s() invalid CRef data\n%s\n%s",
         FCT,STR(A_),STR(B_));
      return init();
   }

   const CDATA_TQ &A=(*A_.cgb), &B=(*B_.cgb);

   unsigned na,nb, i,i_,j, l=0, nq=A.t.qlen(),
      ra=A.rank(F_L), rb=B.rank(F_L), 
      nab=ra+rb;
   const wbperm &pa=A_.cgp, &pb=B_.cgp;
   char ma[nab], *mb=ma+ra; memset(ma,0,nab*sizeof(char));

   char sa=((A_.conj!=0) ^ (ica.conj!=0) ? -1 : +1),
        sb=((B_.conj!=0) ^ (icb.conj!=0) ? -1 : +1), xd=-sa*sb;

   if (A.t!=B.t || !nq) wblog(F_L,
      "ERR %s() got symmetry mismatch (%s, %s)",FCT,QSTR(A),QSTR(B));

   if (pa.len && !pa.isValidPerm(ra)) wblog(F_L,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pa.len,ra);
   if (pb.len && !pb.isValidPerm(rb)) wblog(F_L,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,pb.len,rb);

   if (ica.len!=icb.len) wblog(F_L,"ERR %s() "
      "invalid set of contraction indices (%d/%d)",FCT,ica.len,icb.len);
   if (A.qdir.len!=ra) wblog(F_L,
      "ERR %s() invalid qdir set (A: %d/%d)",FCT,A.qdir.len,ra);
   if (B.qdir.len!=rb) wblog(F_L,
      "ERR %s() invalid qdir set (B: %d/%d)",FCT,B.qdir.len,rb);

   for (i=0; i<ica.len; ++i) { j=pa.el1(ica[i]);
      if (j>=ra) wblog(F_L, 
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,ica.len);
      if ((++ma[j])>1) wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,ica.len
      );
   }
   for (i=0; i<icb.len; ++i) { j=pb.el1(icb[i]);
      if (j>=rb) wblog(F_L, 
         "ERR %s() index out of bounds (%d/%d)",FCT,j+1,icb.len);
      if ((++mb[j])>1) wblog(F_L,
         "ERR %s() index not unique (%d/%d)",FCT,j+1,icb.len
      );
      if (A.qdir[pa.el1(ica[i])] != xd*B.qdir[j]) wblog(F_L,
         "ERR invalid CRef contraction [conj=-(%d,%d)(%d,%d)=%d]\n"
         "having %s @ [%s](%s) <> %s @ [%s](%s)\n"
         "hint: only pairs of {in/out} indices accepted",
         A_.conj, ica.conj, B_.conj, icb.conj, xd,
         A.qdir.toStr().data, STR(pa), STR(ica),
         B.qdir.toStr().data, STR(pb), STR(icb)
      );
   }

   na=ra-ica.len;
   nb=rb-icb.len; nab=na+nb;

   t=A.t; qdir.init(nab); qs.init(nq*nab);
   if (Pcgd) Pcgd->init(nab);

   if (!nab) return *this;
   if (nab>127) wblog(F_L,
      "ERR %s() char index out of bounds (%d)",FCT,nab);

   const char *a=A.qdir.data, *b=B.qdir.data;
   const TQ *qa=A.qs.data, *qb=B.qs.data;
   char *c=qdir.data; TQ *qc=qs.data;

   if (Pcgd) {
      for (i=0; i<ra; ++i) if (!ma[i]) { ma[i]=-char(l); ++l; }
      for (i=0; i<rb; ++i) if (!mb[i]) { mb[i]=-char(l); ++l; }

      if (l!=nab) wblog(F_L,
         "ERR %s() %d != %d + %d",FCT,l,na,nb);
      l=0;

      if (na) { for (i=0; i<ra; ++i) {
         j=pa.el1(i); if (ma[j]<=0) { Pcgd->data[l++]=-ma[j]; }
      }}
      if (nb) { for (i=0; i<rb; ++i) {
         j=pb.el1(i); if (mb[j]<=0) { Pcgd->data[l++]=-mb[j]; }
      }}
      l=0;
   }

   if (nq==1) {
      if (ra) { for (i_=0; i_<ra; ++i_) { i=pa.el1(i_); if (ma[i]<=0) {
         c[l]=sa*a[i]; qc[l]=qa[i]; ++l;
      }}}
      if (rb) { for (i_=0; i_<rb; ++i_) { i=pb.el1(i_); if (mb[i]<=0) {
         c[l]=sb*b[i]; qc[l]=qb[i]; ++l;
      }}}
   }
   else {
      if (ra) { for (i_=0; i_<ra; ++i_) { i=pa.el1(i_); if (ma[i]<=0) {
         c[l]=sa*a[i];
         for (j=0; j<nq; ++j) { qc[j]=qa[i*nq+j]; }
         qc+=nq; ++l;
      }}}

      if (rb) { for (i_=0; i_<rb; ++i_) { i=pb.el1(i_); if (mb[i]<=0) {
         c[l]=sb*b[i];
         for (j=0; j<nq; ++j) { qc[j]=qb[i*nq+j]; }
         qc+=nq; ++l;
      }}}
   }

   if (l!=nab) wblog(FL,"ERR %s() %d/%d",FCT,l,nab); 

   return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init(const QType &t_, const TQ* q, const char *d) {

   unsigned i=0, r=0; 

   if (!d) wblog(FL,"ERR %s() got null qdir string",FCT);
   for (; d[i]; ++i) { if (d[i]!='+' && d[i]!='-') ++r; }

   if (i<2 || r) wblog(FL,"ERR %s() invalid qdir string '%s'",FCT,d);
   r=i; 

   qs.init(r,q);
   qdir.init2val(r,+1); for (i=0; i<r; ++i) { if (d[i]=='-') qdir[i]=-1; }

   t=t_; return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType &t_, const qset<gTQ>&qs_,
   unsigned iout, 
   char isref
){
   if (isref)
        { qs.init2ref(qs_); }
   else { qs=qs_; }

   unsigned n=t_.qlen(), r=qs.len/n;
   if (qs.len%n) wblog(FL,
      "ERR %s() invalid qset (%d @ %d)",FCT,qs.len,n);

   qdir.init_iout(FL, iout? r : 0, iout);

   t=t_; return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType t_, const TQ* qs_, unsigned r, unsigned N,
   unsigned iout 
){
   unsigned n=t_.qlen();
   t=t_; qs.init(r*n);

   if (qs_ && qs.len) Wb::cpyStride(qs.data,qs_,n,NULL,r,-1,N);

   qdir.init_iout(FL, iout? r : 0, iout);

   return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init(
   const QType t_, const TQ* qs_, unsigned r, unsigned N,
   const iTags &itags
){
   unsigned n=t_.qlen();

   t=t_; qs.init(r*n);
   if (qs_ && qs.len) Wb::cpyStride(qs.data,qs_,n,NULL,r,-1,N);

   if (itags.len!=r) {
      if (!itags.len) wblog(FL,
         "ERR %s() got empty itags (%d/%d)",FCT,itags.len,r);
      else wblog(FL,"ERR %s() invalid itags "
         "[%d/%d: %s]",FCT,itags.len,r,STR(itags)
      );
   }
   qdir.init(itags);

   return *this;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::init_str(const char *F, int L, const char *s_) {

   unsigned i=0, nd=0, nc=0, r=0, rs, nw=0, n, e=0;
   char *s, bflag=0, cflag=0;
   wbstring S(s_); 

   s=S.data; n=(S.len? S.len-1 : 0);

   while (isspace(s[i])) { ++s; --n; }

   for (i=0; i<n; ++i) { if (isspace(s[i])) { s[i]=0; break; }}
   if (i>=n || i<2 || i>6) wblog(F_L,
      "ERR %s() invalid QSet specification (QType)\ns='%s'",FCT,s);

   t.init(s); ++i; s+=i; n-=i;

   for (i=0; isspace(s[i]); ++i);

   if (index("({[",s[i])) {
      bflag=s[i]; for (++i; isspace(s[i]); ++i);
   }; if (i) { s+=i; n-=i; }

   for (i=0; i<n; ++i) {
      if (isspace(s[i])) { ++nw; }
      else {
         if (isalnum(s[i])) { ++nd;     } else
         if (s[i]==',') { ++r;          } else
         if (s[i]==';') { ++r; ++cflag; } else
         if (s[i]=='*') { ++nc;         } else
         if (bflag) {
            if ((s[i]==')' && bflag=='(') ||
                (s[i]=='}' && bflag=='{') ||
                (s[i]==']' && bflag=='[')
            ) { bflag=0; continue; } else e=1;
         }
         else e=2;
         if (e) wblog(F_L,"ERR %s() "
            "invalid qset string\n%s (%d/%d; e=%d)",FCT,s,i+1,n,e
         );
      }
   }

   if ((nc && cflag) || cflag>1) wblog(FL, 
      "ERR %s() invalid QSet '%s' (%d,%d)",FCT,s_?s_:"",nc,cflag);
   if (!nd) wblog(FL,"ERR %s() got invalid QSet (nd=0)",FCT);

   rs=t.qrank(); ++r; cflag=0;

   if (nw==0 && nd%r==0 && r*rs==nd) { unsigned l=0;
      qdir.init2val(r,+1); qs.init(nd);
      for (r=i=0; i<n; ++i) {
         if (!CG::cstr2qset(s+i,qs.data+l,1)) { ++l; } else
         if (s[i]==',') { ++r; } else
         if (s[i]==';') { cflag=(++r); } else
         if (s[i]=='*') { qdir[r]=-1; }
      }
   }
   else { unsigned l=0;
      qdir.init2val(r,+1); qs.init(r*rs);

      for (r=i=0; i<n; ++i) {
         if (isalpha(s[i])) { e=1; break; }
         if (isdigit(s[i])) { if (l>=qs.len) { e=2; break; }
            if (sscanf(s+i,"%d",qs.data+(l++))!=1) { e=3; break; }
            for (++i; isdigit(s[i]); ++i);
         }
         if (s[i]==',') { ++r; } else
         if (s[i]==';') { cflag=(++r); } else
         if (s[i]=='*') { if (r<qdir.len) qdir[r]=-1; else { e=4; break; }}
      }

      if (e) wblog(F_L, 
         "ERR %s() invalid QSet (%d,%d,e=%d; %d = %d*%d)\n%s",
         FCT,nc,nw,e, nd,qdir.len,rs, s_);
   }

   if (cflag) { for (i=cflag; i<qdir.len; ++i) { qdir[i]=-qdir[i]; }}

   return *this;
};

template <class TQ>
bool QSet<TQ>::isScalar() const {

   unsigned i=0, r=rank(FL), r2=r/2, m=t.qlen(), n2=r2*m;

   if (r<2) wblog(FL,"ERR %s() got r=%d QSet",FCT,r);
   if (qdir.len!=r || qs.len!=r*m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,STR_(this));

   for (; i<r2; ++i) { if (qdir[i]<=0 || qdir[i+r2]>=0) return 0; }

   for (i=0; i<n2; ++i) { if (qs[i]!=qs[i+n2]) return 0; }

   if (t.isNonAbelian() && r%2) {
   for (i=2*n2; i<qs.len; ++i) { if (qs[i]!=0) return 0; }}

   return 1;
};

template <class TQ>
bool QSet<TQ>::is1JSymbol() const {

   unsigned r=rank(FL), m=t.qlen();

   if (r<2) wblog(FL,"ERR %s() got r=%d QSet",FCT,r);
   if (qdir.len!=r || qs.len!=r*m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,STR_(this));
   if (r!=2) return 0;

   for (unsigned i=0; i<r; ++i) { if (qdir[i]<=0) return 0; };

   return t.isDual(qs.data,qs.data+m);
};

template <class TQ>
bool QSet<TQ>::gotRCData() const {

   wbstring fs;
   if (gStore.get_file_name(FL,fs,*this,"cgd")>0 &&
      Wb::getFileSize(fs.data)>0) return 0;
   return -1;
};

template <class TQ>
int QSet<TQ>::isZero() const {

   const TQ *q=qs.data;
   unsigned i=0, i1=0, i2=1, id=0, l,m, r=rank(FL);

   if (!r) { return 0; } 
   m=qs.len/r;

   if (r==1) {
      return (qs.norm2() ? 1 : 0);
   }
   else if (r==2) { l=qs.len-1;
      if (qdir=="+-") { id=1; } else
      if (qdir=="++" || qdir=="--") id=2; else
      wblog(FL,"ERR %s() invalid (unsorted) qdir='%s'",FCT,STR(qdir));
   }
   else if (r==3) {
      for (i=0; i<r; ++i, q+=m) {
         for (l=0; l<m; ++l) { if (q[l]) break; }
         if (l==m) break; 
      }; q=qs.data;

      if (i<r) { QDir qd(2);
         if (i==0) { i1=1; i2=2; } else
         if (i==1) {       i2=2; } 

         qd[0]=qdir[i1]; qd[1]=qdir[i2];
         if (qd=="+-") { id=1; } else
         if (qd=="++" || qdir=="--") id=2; else
         wblog(FL,"ERR %s() invalid (unsorted) qdir='%s'",FCT,STR(qdir));
      }
   }

   if (id) { i1*=m; i2*=m;
      if (id==1) {
         for (i=0; i<m; ++i) {
         if (q[i+i1]!=q[i+i2]) { return (10*r+1); }}
      }
      else {
         if (!t.isDual(q+i1,q+i2)) { return (10*r+2); }
      }
      return 0; 
   }

   if (r==1) { 
      if (qs.norm2()==0) return 0;
      wblog(FL,"ERR %s() got rank-%d QSet %s",FCT,r,STR_(this));
   }

   if (r!=3) { return gotRCData() ? 0 : -1; }

   if (qdir.len!=3 || qdir[0]<=0 || qdir[1]<=0) wblog(FL,
      "ERR %s() got invalid (unsorted) QSet %s",FCT,STR_(this));
   m=qs.len/r;

   int n3; qset<TQ> J1,J2,J3,J12;

   J1.init2ref(m,qs.data );  J3 .init2ref(m,qs.data+2*m);
   J2.init2ref(m,qs.data+m); J12.init2ref(2*m, qs.data );

   n3=gCS.find_map3(t,J12); 
   if (n3<=0) {
      n3=gStore.load_Std3(0,0,t,J1,J2,0); 
   }

   if (n3<=0) {
      return gotRCData() ? 0 : -1;
   }

   char rflag=0;    
   if (qdir[2]>0) { 
      if (t.isSUN()) rflag=1; else
      if (!t.isSelfDual()) wblog(FL,"ERR %s() "
         "got non-self-dual %s to checkout mp3 data",FCT,STR(t)
      );
   }
   if (gCS.valid_mp3_data(t,J1,J2,J3,rflag)) { return 0; }

   return (10*r+3);
};

template <class TQ> 
QSet<TQ>& QSet<TQ>::Sort(wbperm *cgp, char *conj, char iflag) {

   if (isEmpty() || isSorted()) {
      if (cgp) { cgp->init(); }; if (conj) { (*conj)=0; }
      return *this;
   }

   unsigned i=0, l=0, m=t.qlen(), r=(m ? qs.len/m : 0);
   wbMatrix<TQ> dQ(r,1+m); 
   char cflag=0, *d=qdir.data;
   TQ *dq=dQ.data;
   wbperm P;

   if (!m || r*m!=qs.len || qdir.len!=r || r>127) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,STR_(this));

   if (r<2) {
      if (r==1) {
         if (d[0]<=0) { cflag=1; d[0]=(-d[0]); } 
         if (conj) { (*conj)=cflag; }
         if (cgp ) { cgp->init();   };
         return *this;
      }
      wblog(FL,"ERR %s() got rank-%d QSet",FCT,r);
   }

   for (; i<r; ++i, dq+=dQ.dim2) {
      if (d[i]>0) ++l; else
      if (!d[i]) wblog(FL,"ERR %s() invalid qdir\n%s",FCT,STR_(this));
      dq[0]=(d[i]>0 ? 0 : 1); 
      MEM_CPY<TQ>(dq+1, m, qs.data+i*m);
   }

   if (l<(r+1)/2) {
      cflag=1; dq=dQ.data;
      for (i=0; i<r; ++i, dq+=dQ.dim2) { dq[0]=(!dq[0]); }
   }
   dQ.SortRecs(P);

   if (2*l==r) { int q=0;
      dq=dQ.data+1;
      for (i=0; i<l; ++i, dq+=dQ.dim2) {
         if ((q=Wb::cmpRange(dq,dq+l*dQ.dim2,m))) {
            if (q>0) { P.Rotate(l); cflag=(2+!cflag); }
            else { i=r; }; break;
         }
      }
   }

   if (cflag%2) {
      if (conj) { (*conj)=1; }
      for (i=0; i<r; ++i) d[i]=-d[i]; 
   }
   else {
      if (conj) { (*conj)=0; }
   }

   qdir.Permute(P);

   if (P.isIdentityPerm()) {
      if (cgp) { cgp->init(); }
      return *this;
   }
   else if (cgp) {
      if (iflag) {
         P.invert(*cgp);
      }
      else P.save2(*cgp);
   }

   dq=dQ.data+1; if (cflag<2) l=0;
   for (i=0; i<r; ++i) {
      MEM_CPY<TQ>(qs.data+i*m, m, dq+((i+l)%r)*dQ.dim2);
   }

   return *this;
};

template <class TQ>
bool QSet<TQ>::isSorted() const {

   if (isEmpty() || !qdir.isSorted()) { return 0; }
   if (!qdir.len) {
      if (qs.len) wblog(FL,
         "ERR %s() got invalid empty QSet\n%s",FCT,STR_(this));
      return 1; 
   }

   unsigned i=0, l=0, m=t.qlen(), r=(m ? qs.len/m : 0);
   if (!m || qs.len%m) wblog(FL,
      "ERR %s() got invalid QSet %s",FCT,STR_(this));

   if (r<2) { 
      if (r==1) { 
         return (qdir[0]>0 ? 1 : 0);
      }
      wblog(FL,"ERR %s() got rank-%d QSet %s",FCT,r,STR_(this));
   }

   for (; i<r; ++i) { if (qdir[i]<=0) break; }
   if (i>1) {
      wbMatrix<TQ> X(i,m,qs.data,'r'); 
      if (!X.isSorted(+1)) return 0;
   }

   for (l=i; i<r; ++i) {
      if (!qdir[i]) wblog(FL,
         "ERR %s() got invalid qdir=%s",FCT,STR(qdir));
      if (qdir[i]>0) return 0;
   }

   if (l+1<r) {
      wbMatrix<TQ> X(r-l,m,qs.data+l*m,'r'); 
      if (!X.isSorted(+1)) return 0;
   }

   if (l<(r+1)/2) return 0;
   if (2*l==r) { m*=l;
      if (Wb::cmpRange(qs.data, qs.data+m,m)>0) return 0;
   }

   return 1;
};

template <class TQ>
int QSet<TQ>::checkQ_abelian(const char *F, int L) const {

   unsigned i=0; TQ q=0;

   if (qdir.len!=qs.len) { if (F) wblog(FL,
      "ERR %s() qdir inconsistency (%s)",FCT,STR_(this));
      return -3;
   }
   if (!qdir.len) { return 0; } 

   for (; i<qdir.len; ++i) {
      if (!qdir[i]) wblog(FL,"ERR %s() got qdir=%s",FCT,STR(qdir));
   }

   switch (t.type) {
      case QTYPE_U1:
      case QTYPE_ZN:

         for (i=0; i<qdir.len; ++i) {
            if (qdir[i]>0) { q+=qs[i]; } else { q-=qs[i]; }};
         if (q && t.type==QTYPE_ZN) { q=int(q)%t.sub; }
         break;

      case QTYPE_P :
         q=qs.prod()-1; break;

      default:
         if (F) wblog(FL,"ERR %s() got %s [%s]",FCT,STR(t),STR(qs));
         return -1;
   }

   if (q) {
      if (F) wblog(F,L,
         "ERR %s() got broken %s symmetry [%s]",FCT,STR(t),STR(qs));
      return 2;
   }
   return 0;
};

template <class TQ>
int QSet<TQ>::checkQ_SU2(const char *F, int L) const {

   if (t.type!=QTYPE_SUN || t.sub!=1) { if (F) wblog(F,L,
      "ERR %s() got invalid type %s",FCT,STR(t));
      return -1;
   }
   if (qdir.len<2) { if (F) wblog(FL,
      "WRN %s() got rank-%d QSet (%s)",FCT,qdir.len,STR_(this));
      return -2;
   }
   if (qdir.len!=qs.len) { if (F) wblog(FL,
      "ERR %s() qdir inconsistency (%s)",FCT,STR_(this));
      return -3;
   }

   unsigned i=0, nin=0;

   for (; i<qdir.len; ++i) {
      if (qdir[i]>0) { ++nin; } else
      if (!qdir[i]) wblog(FL,
     "ERR %s() got undetermined qdir=%s",FCT,STR(qdir));
   }

   if (!nin || nin==qdir.len) {
      return 2;
   }

   if (qdir.len==2) {
      if (qs[0]!=qs[1]) { if (F) wblog(FL,
         "ERR %s() got scalar CGC data %s",FCT,STR_(this));
         return 3;
      }
   }
   else if (qdir.len==3) {
      unsigned i1,i2,i3=0;
      if (nin==1)
         for (i=0; i<qdir.len; ++i) { if (qdir[i]>0) { i3=i; break; }}
      else
         for (i=0; i<qdir.len; ++i) { if (qdir[i]<0) { i3=i; break; }}

      if (i3==0) { i1=1; i2=2; } else
      if (i3==1) { i1=0; i2=2; }
      else       { i1=0; i2=1; }

      if (qs[i1]<qs[i2]) { i=i1; i1=i2; i2=i; } 

      if (qs[i3] > qs[i1]+qs[i2] || qs[i3] < qs[i1]-qs[i2]) {
         if (F) wblog(F,L,"ERR %s() "
            "invalid SU2 addition rule (%g + %g = %g!)\n%s",FCT,
            double(qs[i1]),double(qs[i2]),double(qs[i3]),STR_(this));
         return 4;
      }
   }

   return 0;
};

template <class TQ>
QSet<TQ>& QSet<TQ>::permute(
   QSet<TQ> &B, const wbperm &P, char iflag) const {

   if (P.len) {
      unsigned r=rank(FL);
      if (!P.isValidPerm(FL,r)) wblog(FL,
         "ERR invalid permute on CData (%d/%d)",P.len,r);

      B.t=t;
      qdir.permute(B.qdir,P,iflag);
      qs.blockPermute(P,B.qs,iflag);
   }
   else { B=*this; } 

   return B;
};

template <class TQ>
wbstring QSet<TQ>::QStrS(const wbperm *cgp, char s_) const {

   wbstring s; 

   if (!qs.len && !qdir.len) {
      if (!cgp || !cgp->len) { s=""; return s; } else
      wblog(FL,"ERR %s() got cgp=[%s] on empty QSet",FCT,STR_(cgp));
   }

   unsigned k=0, m=t.qlen();

   QSet<TQ> Q;
   if (cgp)
        { permute(Q,*cgp); }
   else { Q.init2ref(*this); }

   k=Q.qdir.isSorted();
   if (!k) wblog(FL,"WRN %s() got unsorted QSet",FCT);

   if (CG::gotQAlpha(qs.data,qs.len)) {
      unsigned i=0, l=0, r=qdir.len; TQ *q=Q.qs.data;
      if (r*m!=qs.len || !r) wblog(FL,
         "ERR %s() got invalid QSet (%d*%d=%d)",FCT,r,m,qdir.len);
      s.init(qs.len+qdir.len); k-=2;

      for (i=0; i<r; ++i, q+=m) {
         if (CG::qset2cstr(q,s.data+l,m)) wblog(FL, 
            "ERR %s() failed to obtain compact qset string",FCT);
         l+=m; s[l++] = (i!=k ? ',' : s_);
      }
      s[l-1]=0; 
   }
   else {
      s=Q.qs.toStrf("", m<2 ? "":" ", t.qlen(),";");
      if (k) {
         unsigned i=0, r=0, n=s.len; char *c=s.data;
         for (--k; i<n; ++i) if (c[i]==';') {
            if ((++r)!=k) { c[i]=','; } else { c[i]=s_; }
         }
      }
   }

   return s;
};

template <class TQ>
wbstring QSet<TQ>::QStr() const {
   wbstring sout; 

   unsigned i=0, j=0, l=0, n, m=t.qlen();
   const TQ *q=qs.data; char *s;

   if (qs.len!=m*qdir.len) {
      if (qdir.len) wblog(FL,"ERR %s() "
         "invalid QSet (%d*%d=%d)",FCT,m,qdir.len,qs.len);
      if (!m || (!t.isUnknown() && !t.isAbelian())) wblog(FL,
         "ERR %s() QSet inconsistency (%s @ qlen=%d!=%d*%d)",
         FCT,STR(t),qs.len,m,qdir.len
      );
      sout=qs.toStrf(""," ",m,";");
      return sout;
   }
   if (!qdir.len) { sout.init(1); return sout; }

   if (CG::gotQAlpha(qs.data,qs.len)) {  
      n=qs.len+2*qdir.len; 
      sout.init(n+1); s=sout.data;
      for (; i<qdir.len; ++i, q+=m) {
         if (CG::qset2cstr(q,s+l,m)) wblog(FL, 
            "ERR %s() failed to obtain compact qset string",FCT);
         l+=m; if (qdir[i]<0) { s[l++]='*'; }
         s[l++]=',';
      }
      s[l-1]=0; 
   }
   else {
      wbstring fmt; fmt.init2Fmt(qs[0]);
      n=5*qs.len+qdir.len; sout.init(n+1); s=sout.data;

      for (; i<qdir.len && l<n; ++i, q+=m) { if (i) { s[l++]=','; }
         for (j=0; j<m && l<n; ++j) {
            if (j) { s[l++]=' '; }
            l+=snprintf(s+l,n-l,fmt.data,q[j]);
         }
         if (qdir[i]<0 && l<n) { s[l++]='*'; }
      }

      if (l>n) { s[n]=0; wblog(FL,
         "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,s); }
      sout.len=l+1;
   }

   return sout;
};

template <class TQ>
wbstring QSet<TQ>::toStr(const char *istr) const {
   wbstring sout; 

   unsigned m=t.qlen(), n=16+4*qs.len;

   if (istr && !istr[0]) istr=0;

   if (qs.len!=m*qdir.len) {
      if (qdir.len) wblog(FL,"ERR %s() severe QSet inconsistency "
         "(len=%d!=%d*%d)",FCT,qs.len,m,qdir.len);
      if (!m || (!t.isUnknown() && !t.isAbelian())) wblog(FL,
         "ERR %s() QSet inconsistency (%s @ qlen=%d!=%d*%d)",
         FCT,STR(t),qs.len,m,qdir.len
      );
   }
   sout.init(n);

   char *s=sout.data;
   unsigned l=snprintf(s,n,"%5s (",STR(t));

   if (l<n) {
      l+=snprintf(s+l,n-l,"%s",QStr().data);
   }
   if (l<n) {
      if (istr)
           { l+=snprintf(s+l,n-l,") %s",istr); }
      else { s[l]=')'; s[++l]=0; }
   }

   if (l>=n) { s[n]=0; wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,s); }
   sout.len=l+1;

   return sout;
};

template <class TQ>
wbstring QSet<TQ>::toTag() const { wbstring s(64); 

   unsigned l=0, n=s.len-1;

   char sep[2]=",";
   if (t.qlen()<=1 || qs.wbvector<TQ>::allIn(0,9)) sep[0]=0;

   if (!t.validType() || !qdir.len || !qs.len) wblog(FL,
      "ERR %s() got empty CData\n%s",FCT,STR_(this));

   l=snprintf(s.data,n,"%s[%s: %s]",
     t.toStr('t').data, STR(qdir),
     qs.wbvector<TQ>::toStrf("",sep,t.qlen(),",").data
   );

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};

template <class TD>
template <class TQ>
cdata<TD>& cdata<TD>::init(const CRef<TQ> &R, unsigned k, char full) {

   if (!R.cgb) wblog(FL,"ERR %s() got NULL cref",FCT);
   if (!R.cgw) wblog(FL,"ERR %s() got empty cgw",FCT);

   if (R.isRefInit()) { 
      R.LoadRef(FL);
   }

   unsigned r=R.cgb->rank(FL), m=R.cgb->gotOM(FL);
   unsigned m1,m2; m1=R.wdim12(FL,m2);

   if (m<m1 || m1<m2) wblog(FL,
      "ERR %s() inconsistent OM m=%d / %d / %d",FCT,m2,m1,m);
   if (k>=m) wblog(FL,"ERR %s() index out of bounds (k=%d/%d)",FCT,k,m);

   if (m<=1) { init(R.cgb->cgd) *= R.wget0();
      if (m) {
         if (r && this->SIZE[r]==1) {
            this->SkipTrailingSingletons(FL,r-1);
            if (this->SIZE.len!=r-1) { wblog(FL,
               "ERR %s() failed to remove OM singleton dimension",FCT);
            }
         }
         else { wblog(FL,"WRN %s() got "
            "invalid OM=%d CData (%s /%d)",FCT,m, SSTR(R.cgb->cgd), r);
         }
      }
   }
   else {
      wbvector<double> wk;
      wk.initp(m, R.cgw.col(k), m1, 1); 
      R.cgb->cgd.cgsparray::contract(FL,r,wk, (cgsparray&)(*this));
   }

   if (full) {
      if (R.cgp.len) { Permute(R.cgp); }
      if (R.conj) { Conj(); }
   }
   return *this;
};

template <class TD>
cdata<TD>& cdata<TD>::initIdentity(
   const QType &t, SPIDX_T d, TD dval
){
   if (d!=1 || t.isNonAbelian())
        wbsparray<TD>::initIdentity(d,dval);
   else wbsparray<TD>::initScalar(dval);

   return *this;
};

template <class TD>
cdata<TD>& cdata<TD>::initIdentity(const wbvector<SPIDX_T> &S) {

   if (S.len==1 && S[0]==1) {
      wbsparray<TD>::initScalar(1.); }
   else if (S.len==2 && S[0]==S[1]) {
      wbsparray<TD>::initIdentity(S[0]);
   }
   else wblog(FL,
      "ERR wbsparse::%s() non-square (%s)",FCT,SSTR(S));

   return *this;
};

template <class TD>
unsigned cdata<TD>::getOM(
   const char *F, int L, const char *istr) const {

   double x2=double(this->norm2()), e=std::fabs(x2-std::round(x2));
   unsigned m=(x2+0.5), r=this->SIZE.len;

   if (e>1E-12) wblog(F_L,"ERR %s() "
      "got non-integer |cgc|^2 (%s%g @ %.3g)",FCT,istr?istr:"",x2,e);
   if (!r) {
      if (this->isDiag()) { 
         if (m!=1) wblog(FL, "ERR %s() got |cgc|^2=%g\n"
            "for diagonal CData (D.len=%d)",FCT,x2,this->D.len); }
      else if (m) wblog(F_L, 
         "ERR %s() got om=%g for empty CData",FCT,x2);
      return 0;
   }

   if (m<=1) return 0;

   const SPIDX_T *sd=this->SIZE.data;

   for (unsigned s=1, i=r-1; i<r; --i) { 
      if ((s*=sd[i])==m) {
         if (r-i>1) wblog(FL, 
            "WRN %s() got %d OM dimensions",FCT,r-i);
         return (r-i);
      }
      if (s>m) wblog(F_L,"ERR %s() OM cgn2=%g "
         "incompatbile with SIZE=[%s]",FCT,x2,this->sizeStr().data
      );
   }

   return 0;
};

template <class TD>
TD cdata<TD>::NormSignC( 
   const char *F, int L,
   unsigned m, 
   TD eps1, TD eps2
){
   TD nrm; 

   if (this->D.len<=1) { if (!this->D.len) return 0; else {
      TD a=ABS(this->D[0]); 
      if (m) wblog(FL,"ERR %s() got m=%d for scalars",FCT,m);
      if (a<eps1) {
         if (a>eps2 && F) wblog(F_L, 
            "WRN %s() CGC noise %.3g (1: %g, %g)",
            FCT, double(this->D[0]), double(eps1), double(eps2));
         if (double(a)<CG_SKIP_EPS2) { this->D[0]=0; return 0; }
      }
      nrm=this->D[0]; this->D[0]=1;
      return nrm;
   }}

   this->Compress(F_L,CG_SKIP_EPS2); 

   SPIDX_T i, n=this->D.len;
   unsigned r=this->rank(), l=this->SIZE.len;
   TD *d=this->D.data; 

   if (m>2 || m>r || r<2 || (r!=l && l && r!=2)) wblog(FL,
      "ERR got OM=%d (r=%d/%d)",m,r,l);

   nrm=this->D.norm2(); 

   if (m) {
      nrm/=Wb::maxRange(this->SIZE.data+(r-m),m);
   }

   nrm=SQRT(nrm);

   if (nrm<eps1) {
      if (nrm>eps2) wblog(FL,
         "WRN %s() CGC noise %.3g (2: %g; %g)",
         FCT, double(nrm), double(eps1), double(eps2));
      return 0;
   }

   for (i=0; i<n; ++i) if (ABS(d[i])>=eps1) {
      if (d[i]<0) { nrm=-nrm; }; break;
   }
   if (i==n) wblog(FL,
      "ERR %s() failed to determine sign (%g)",FCT,double(eps1));

   (this->D)/=nrm;

   CG::FixRational(FL, this->D.data, this->D.len, 4);
   CG::FixRational(FL, &nrm, 1, 4);

   return nrm;
};

template <class TD>
cdata<TD>& cdata<TD>::permute(
   cdata<TD> &B, const wbperm &P0, char iflag
) const {

   if (isScalar() || !P0.len) { B.init(*this); } else
   if (!this->SIZE.len && this->D.len>1) {
      if (!this->isDiag() || P0.len!=2) wblog(FL,"ERR %s() "
         "invalid diagonal (%s, D=%d)",FCT,STR(P0),this->D.len);
      B.init(*this);
   }
   else {
      const unsigned l=this->SIZE.len, r=P0.len;
      if ((l<=2 && l!=r) || (l>2 && (l<r || l>r+1))) wblog(FL,
         "ERR %s() invalid permutation %s [%s @ %.3g]",
         FCT,STR(P0),SSTR_(this),double(this->norm()));

      wbperm P(P0,this->SIZE.len);
      wbsparray<TD>::permute(P,B,iflag);
   }

   return B;
};

template <class TD> 
TD cdata<TD>::contract(const char *F, int L,
   const ctrIdx &ica, const cdata<TD> &B,
   const ctrIdx &icb, cdata<TD> &Cin,
   cdata<TD> *Cx, wbperm *P, char normalize
 ) const {

   TD nrm=1; 
   cdata<TD> C;

   if (this->D.len<=1 && B.D.len<=1) { 
      if (this->contract_scalar(F_L,ica,B,icb,C)) {
         if (Cx) wblog(FL,"ERR %s() got Cx for scalar C",FCT);
         if (!Cin.isEmpty()) wblog(FL, 
            "ERR %s() got non-empty Cin (%s)",FCT,this->sizeStr().data);
         if (normalize) nrm=C.NormSignC(F,L); 
         C.save2(Cin); return nrm;
      }
      else wblog(FL,"WRN %s() " 
         "got non-scalar contraction for D.len=%d/%d",
         FCT,this->D.len, B.D.len
      );
   }

   unsigned ma=this->getOM(FL,"A: "), mb=B.getOM(FL,"B: ");

   this->cgsparray::contract(F,L,ica, B, icb, C, P ? *P : wbperm());

   if (normalize) {
      nrm=C.NormSignC(F,L,ma+mb);
   }

   if (Cx) { Cx->init(); } 

   if (Cin.isEmpty()) { C.save2(Cin); }
   else {
      if (C.SIZE!=Cin.SIZE) {
         if (!C  .SIZE.len) C  .diag2reg(); else
         if (!Cin.SIZE.len) Cin.diag2reg();
         if (C.SIZE!=Cin.SIZE) {
            MXPut(FL,"ix").add(C,"c").add(Cin,"c0").add(nrm,"nrm");
            wblog(FL,"ERR size inconsistency in C.cgs[] (%s; %s)",
            SSTR(C), SSTR(Cin));
         }
      }

      if (fabs(double(nrm))<1E-14) return nrm;

      double d=MAX(SPIDX_T(1),C.numel()), e,e2=std::sqrt(C.normDiff2(Cin)/d);
      if (e2>CG_EPS1) {
         if ((e=double(Cin.norm2()))<CG_EPS2) C.save2(Cin);
         else if (Cx) C.save2(*Cx);
         else {
            MXPut(FL,"i_").add(C,"c").add(Cin,"c0").add(nrm,"nrm");
            wblog(FL,"ERR CGC difference: e=%g\n%s (%s) * %s (%s) => %s",
            e2, this->sizeStr().data, STR(ica+1),
            B.sizeStr().data, STR(icb+1), C.sizeStr().data);
         }
      }
      else if (e2>CG_EPS2) { wblog(FL,"WRN CGC difference: e=%g",e2); }
   }

   return nrm;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::init(const CRef<TQ> &R, unsigned l) {

   cgd.init(R,l); 

   cstat.t=CGD_FROM_CGR;
   if (R.cgp.len)
        R.cgb->QSet<TQ>::permute((QSet<TQ>&)(*this),R.cgp);
   else (*this)=(QSet<TQ>&)(*R.cgb); 

   if (R.conj) this->QSet<TQ>::Conj();

   return *this;
};

template <class TQ, class TD>
int CData<TQ,TD>::RefInit(
   const char *F, int L, const CData<TQ,TD> &B,
   char lflag, 
   char bare
 ){
   int err=0; 

   if (cstat.ID && !cstat.sameID(B.cstat)) {
      gStore.rclog(this->t, PF_L,1,
        "WRN CData::%s() got ID mismatch (#%05x <> #%05x)"
        "%N   %-32s %s%N   %-32s %s%N", FCT, cstat.ID, B.cstat.ID,
         STR_(this),STR2(cstat,'V'), STR(B),STR2(B.cstat,'V'));
      if (!lflag) wblog(FL,"ERR see warning / error above");

      cstat.init(); err=1; 
   }
   cstat=B.cstat;

   cstat.t=( B.cstat!=CGD_EXPLICIT ? CGD_REF_INIT : CGD_EXPLICIT );

   (*this)=(QSet<TQ>&)B; 

   if (bare || !isSymmmetric()) {
      wbvector<SPIDX_T> S; B.getSize(S,bare);
      cgd.wbsparray<RTD>::init(S,0); 
   }
   else { 
      wbvector<SPIDX_T> S;
      wbvector<RTD> cgt; B.getSize(S); B.trace(FL,cgt);
      RefInit_auxtr(F_L,S,cgt);
   }

   return err;
};

template <class TQ, class TD>
int CData<TQ,TD>::Reduce2Ref(const char *F, int L,
   char force
 ){

   if (cstat.t==CGD_REF_INIT) { return 0; }
   if (cstat.t==CGD_EXPLICIT) { return 0; } 

#ifdef QS_USING_OMP
   CG::Guard qLK, qLX;
   if (force) {
      qLK.acquire(FL,*this,"buf");    
      qLX.acquire(FL,*this,NULL,"?");
   }
   else { int i;
      i=qLK.acquire(FL,*this,"buf","1"); if (i<=0) return -1; 
      i=qLX.acquire(FL,*this,NULL,"?1"); if (i<=0) return -1;
   }
#endif

   if (cstat.t==CGD_REF_INIT) { return 0; } 

   if (gotuser_BUF()) { double x2=double(cgd.norm2());
      if (x2<1E-12) {
         wblog(F_L,"ERR %s() got |BUF->CData|^2=%g for %s\n%s",
         FCT,x2,STR_(this));
      }
   }

   setuser_BUF_ref(2); 

#if defined(DBG_QSX_BUF) && (DBG_QSX_BUF & 1)
   BFF.blogf(FL,"--> %-10s %3ld %p  %s%s %p", FCT, gCS.BUF.size(), this,
   STR_(this), cstat.u2Str(" ").data, this->qs.data);
#endif

   wbvector<SPIDX_T> S; getSize(S);
   if (CG_VERBOSE>6) wblog(F_L," *  %s() %s",FCT,STR_(this));

   if (!isSymmmetric()) {
      cgd.wbsparray<RTD>::init(S,0); 
      cstat.t=CGD_REF_INIT;
   }
   else {
      wbvector<RTD> cgt; trace(FL,cgt);
      RefInit_auxtr(FL,S,cgt);
   }

   setuser_BUF_ref(0); 
   return 1;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::RefInit_auxtr(
   const char *F, int L,
   const wbvector<SPIDX_T> &S, const wbvector<RTD> &cgt
){
   unsigned r=this->qdir.len, m=0;

   if (S.len<r || S.len>r+1) wblog(FL,
      "ERR %s() unexpected size [%s] (r=%d)",FCT,STR(S),r);

   if (!cgt.len) {
      if (!S.len) { cgd.init(); return *this; } else
      if (isSymmmetric()) wblog(FL,
         "ERR %s() got empty cgt (S=[%s])",FCT,STR(S));
   }
   else {
      m=(S.len>r ? S[r] : 1);
      if (cgt.len!=m) wblog(FL,"ERR %s() "
         "got OM mismatch (cgt.len=%d/%d)",FCT,cgt.len,m);
      if (m>99) wblog(F_L,"WRN %s() got unexpected OM=%d",FCT,m);
   }

   cgd.wbsparray<RTD>::init(S,m);
   if (m) { cgd.D=cgt;
      if (m>1) {
         unsigned i=1; SPIDX_T *idx=cgd.IDX.data+2*S.len-1;
         for (; i<m; ++i, idx+=S.len) { idx[0]=i; }
      }
   }

   cstat.t=CGD_REF_INIT;

   return *this;
};

template <class TQ, class TD>
int CData<TQ,TD>::LoadRef(const char *F, int L, char force) {

   if (cstat!=CGD_REF_INIT && !force) {
      return 0;
   }

   CDATA_TQ *C= &(gCS.getBUF(F_L, (QSet<TQ>&)*this, LB_UPD));

   if (C!=this) wblog(FL,"ERR %s() got different CData space"
      "\n > %lx\n > %lx",FCT,this,C);
   if (C->isEmpty() || C->isRefInit()) wblog(FL, 
      "ERR %s() got non-existing CData\n  %s\n> %s",
      FCT,STR_(this), C->toStr().data
   );

   return 1;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::initX3(const char *F, int L,
   const QSet<TQ> &Q0, wbvector<char> &cflags, wbperm &P,
   unsigned loadRC,   
   unsigned char flag 
){
   unsigned m;
   bool store=(flag & CX3_SAVE ? 1 : 0);

   bool force=(flag & CX3_FORCE ? 1 : 0);

   CDATA_TQ &C0=gCS.getBUF(
      force ? (F ? F : __FILE__) : 0,
      force ? (L ? L : __LINE__) : 0, Q0, loadRC);

   if (C0.isEmpty()) {
      if (force && F) wblog(F,L,"ERR %s[%s] got empty C0",FCT,STR(Q0));
      this->init(C0); return *this;
   }
   else if (C0.isRefInit()) { 
      C0.LoadRef(FL);
   }

   C0.permute(*this,P); 
   if (cflags.len) {
      if (cflags.len!=Q0.qdir.len) wblog(FL,
         "ERR %s() size mismatch (%d/%d)",FCT,cflags.len,Q0.qdir.len);
      for (unsigned i=0; i<cflags.len; ++i) {
         if (cflags[P.el1(i)]) this->Conj(i);  
      }
   }

   m=cgd.SIZE.len - Q0.qdir.len;
   if (!cgd.SIZE.len) wblog(FL,
      "ERR %s() got empty cgd [%s]",FCT,SSTR(cgd));
   if (m>1) wblog(FL,
      "ERR %s() got om=%d [%s]\n%s",FCT,m,SSTR(cgd),SSTR(Q0));

   { double x=double(cgd.NormSignC(F,L,m));  
     if (fabs(fabs(x)-1)>1E-12) wblog(FL,
        "ERR %s() norm changed by factor %.5g/%.5g\n  %s "
        "p=(%s) c=[%s] l=%s,%s\n  %s\n> %s",FCT, x, double(C0.cgd.norm()),
        STR(Q0), STR(P), cflags.toStr(-1,"").data, BITS(loadRC),
        cSTR(flag), STR(C0), STR_(this)
     );
   }

   cstat.init(); 
   cstat.init(C0.cstat==CGD_STD3 ? CGD_STD3_X : CGD_GEN_X);

   gStore.rclog(this->t, PF_L, CG_VERBOSE>6,"C3X initialize "
     "rank-3 tensor via 1j+perm trafo\n--> %s P=(%s), c=[%s]",
     STR_(this),STR(P),cflags.toStr(-1,"").data); 

   if (store)
   gStore.save_CData(FL,*this);  

   return *this;
};

template <class TQ, class TD>
unsigned CData<TQ,TD>::rank(const char *F, int L) const {

   unsigned r=0; 
   if (isEmpty()) return r; 

   unsigned i1=-1, i2=-1; 
   unsigned m=this->t.qlen();
   unsigned l=this->qs.len, r_=cgd.SIZE.len, r0=this->qdir.len;

#ifdef QS_USING_OMP
   if (gotuser_BUF_ref() || !r_) {
      CG::Guard qLK(FL,*this,"buf","?");  i1=qLK.used(); 
      CG::Guard qLX(FL,*this, NULL,"?9"); i2=qLX.used();

      l=this->qs.len; r_=cgd.SIZE.len; r0=this->qdir.len; 
   }
#endif

   if (!l && !r_ && !r0) {
      if (gotuser_BUF()) wblog(F_L,
         "WRN %s() got empty CData (r=%d) %s",FCT,r,STR_(this));
      return r;
   }

   if (!m || l%m) wblog(F_L,
      "ERR %s() invalid CData (%d @ %d)",FCT,l,m);
   if ((r=l/m)!=r0) wblog(F_L,
      "ERR %s() inconsistent qdir (%d/%d) %s",FCT,r,r0,STR_(this));

   if (r_==r+1) {
      if (r<3 || !this->t.permitsOM(r)) {
         if (cgd.SIZE[r]!=1) { MXPut(FL,"y").add(cgd,"c");
            wblog(F_L,"ERR %s() got OM for rank-%d CGS\n%s: %s",
            FCT, r, STR_(this), cgd.sizeStr().data);
         }
      }
   }
   else if (r_!=r) {
      if (r==2 && cgd.isDiag()) return r; 
      if (!isAbelian()) {
         MXPut(FL,"Irk") 
           .add(*this,"A").add(cgd,"cgd").add(r_,"r_").add(r,"r")
           .add(i1,"i1").add(i2,"i2").add(cgd.isDiag(),"isd");
         wblog(F_L,"ERR CData %s() rank mismatch (%s r=%d/%d/%d)",
         FCT, qStr().data, r_, cgd.SIZE.len, r);
      }
   }
   return r;
};

template <class TQ, class TD>
template <class T>
wbvector<T>& CData<TQ,TD>::getSize(
   wbvector<T> &S, char bare, const wbperm *cgp) const {

   unsigned i=0, l=(cgp ? cgp->len : 0);

#ifdef QS_USING_OMP
   unsigned l0=cgd.SIZE.len;

   CG::Guard qLK;
   if (gotuser_BUF()) {
      int q=qLK.acquire(FL,*this,"buf","?");
      if (q<=0) qLK.acquire(FL,*this,NULL,"?8");
   }
#endif

   S.init(cgd.SIZE.len);

   if (!S.len) {
      cgd.wbsparray<TD>::getSize(S);
      if (S.len==2 && this->qdir.len==2) {
         if (l && l!=2) wblog(FL,
            "ERR %s() invalid cgp=[%s] (2)",FCT,cgp->toStr().data);
         return S;
      }
      if (!S.len && isScalar()) {
         if (l && l!=this->qdir.len) wblog(FL,"ERR %s() "
            "invalid cgp=[%s] (%d)",FCT,cgp->toStr().data,this->qdir.len);
         return S.init2val(this->qdir.len,1);
      }
      if (S.len || l) wblog(FL,
         "ERR %s() got %s",FCT,STR_(this));
      return S;
   }

   const SPIDX_T *sz=cgd.SIZE.data; 

   if (l) {
      if (l+1<S.len || l>S.len) wblog(FL,
         "ERR %s() cgp/SIZE mismatch (%d/%d)",FCT,l,S.len);
      if (!cgp->isIdentityPerm()) {
         for (; i<l; ++i) { S.data[i]=sz[cgp->data[i]]; }
      }
   }

   for (; i<S.len; ++i) { S.data[i]=sz[i]; }

   if (bare) { i=0;
      if (S.len==this->qdir.len+1) {
         if (S.len>3) S.len=this->qdir.len; else i=1;
      }
      else if (S.len!=this->qdir.len) i=2;

      if (i) wblog(FL,
         "ERR %s() got invalid OM setting (%d/%d)\n%s",
         FCT,S.len,this->qdir.len, STR_(this)
      );
   }

#ifdef QS_USING_OMP 
   if (sz!=cgd.SIZE.data) {
      wbstring s1=CG::Guard::Status(*this,"buf");
      wbstring s2=CG::Guard::Status(*this);
      wblog(FL,"ERR %s() cgd.SIZE changed location\n"
      "%p -> %p (len=%d -> %d)%N%N%s%N%s%N",
      FCT,sz,cgd.SIZE.data,l0,cgd.SIZE.len,s1.data,s2.data);
   }
#endif

   return S;
};

template <class TQ, class TD>
SPIDX_T CData<TQ,TD>::dim() const {

   if (this->qdir.len!=2) wblog(FL,
      "ERR %s() got %s (r=%d)",FCT,STR_(this),this->qdir.len);

   if (cgd.SIZE.len) {
      const SPIDX_T *s=cgd.SIZE.data;
      if (cgd.SIZE.len!=2 || s[0]!=s[1] ||
         (this->t.isAbelian() && s[0]!=1)) wblog(FL,
         "ERR %s() unexpected scalar %s (%s)",FCT,STR_(this),SSTR(cgd));
      return s[0];
   }

   if (this->t.isAbelian()) {
      if (cgd.D.len) {
         if (cgd.D.len>1) wblog(FL,
            "ERR %s() got %s (D.len=%d)",FCT,STR_(this),cgd.D.len);
         else wblog(FL,"WRN %s() "
            "got %s (D[0]=%.4g)",FCT,STR_(this),double(cgd.D[0])
         );
      }
      return 1;
   }
   else {
      if (!cgd.D.len) wblog(FL,
         "ERR %s() got %s (empty D)",FCT,STR_(this));
      return cgd.D.len;
   }
};

template <class TQ, class TD>
int CData<TQ,TD>::checkConsistency(const char *F, int L) const {

   if (isEmpty()) return 1;

   unsigned r,d, i=0, l=0, n=this->t.qlen(); int e=0;

   if (!n || !this->qs.len || this->qs.len%n) { e=11; }
   else { r=this->qs.len/n;
      if (!cgd.SIZE.len) { if (!r) return 1;
         if (r==2 && cgd.isDiag()) {
            for (; i<r; ++i, l+=n) {
               d=gRS.qdim(FL,this->t,this->qs.data+l);
               if (d!=cgd.D.len) { e=31+i; break; }
            }
         }
         else return 30;
      }
      else if (cgd.SIZE.len<r || cgd.SIZE.len>r+1) { e=12; } else
      for (; i<r; ++i, l+=n) {
         d=gRS.qdim(FL,this->t,this->qs.data+l);
         if (d!=cgd.SIZE[i]) { e=21+i; break; }
      }
   }
   if (!e && r && cgd.SIZE.len>r) {
      if (int(cgd.SIZE[i])>MAX_OM && CG_VERBOSE) { char s[128];
         snprintf(s,128,"%s => max(OM)=%ld",STR_(this),cgd.SIZE[i]);
         if (cgd.SIZE[i]<256)
              { gStore.rclog(this->t, PF_L,CG_VERBOSE>5 && F,"MOM %s",s); }
         else { gStore.rclog(this->t, PF_L,CG_VERBOSE   && F,"WRN %s",s); }

         while (MAX_OM<int(cgd.SIZE[i]) && MAX_OM<9999) { MAX_OM+=16; }
      }
      if (cgd.SIZE[i]>=9999) e=9999;
   }

   if (e && F) {
      MXPut(FL).add(*this,"cg"); wblog(F_L,
     "ERR %s() invalid CData (e=%d)\n%s",FCT,e,toStr('v').data);
   }

#ifndef WB_SKIP_ASSERT
   if (this->t.isAbelian()) { 
      if (!cgd.isScalar()) { e=999;
          if (F) wblog(FL,"ERR %s() invalid %s CData (%s)",
          FCT, qStr().data, cgd.sizeStr().data);
      }
   }
#endif

   if (!e) {
      if (this->t.isAbelian()) return this->checkQ_abelian(F,L);
      if (this->t.isSU2()    ) return this->checkQ_SU2(F,L);
   }

   return e;
};

template <class TQ, class TD>
int CData<TQ,TD>::checkNormSign(
   const char *F, int L, char xflag) const {

   if (isRefInit() || isAbelian()) { return 0; }

   unsigned r=rank(), l=cgd.SIZE.len, m=getOM(F_L);

   if (m>1 && !this->t.permitsOM(r)) {
       if (F) wblog(F_L,
          "ERR %s() got m=%d for rank-%d/%d CGC in %s",
          FCT,m,r,l,STR(this->t));
      return 1;
   }
   if ((!cgd.isDiag() && (r+1<l || r>l)) || !r) {
      if (F) wblog(F_L, 
         "ERR %s() got m=%d for rank-%d/%d CGC",FCT,m,r,l);
      return 2;
   }
   if (!m || !cgd.D.len) {
      if (F) wblog(FL,
         "ERR %s() got empty CData\n%s",FCT,STR_(this));
      return 3;
   }
   if (r==2 && m>1) {
      if (F) wblog(FL,
         "ERR %s() got invalid rank-2 CGC (m=%d)",FCT,m);
      return 4;
   }
   if (CG::signFirstVal(F_L,cgd.D.data, cgd.D.len)<0) {
      if (F) wblog(F_L,"ERR %s() "
         "got unconventional sign of CData\n%s",FCT,STR_(this));
      return 5;
   }

   if (xflag>3 || xflag<0) { 
      if (xflag=='X') { xflag=2; } else {
         if (xflag!='x' && xflag!='!') wblog(FL,
            "WRN %s() unexpected xflag=%s -> 1",FCT,cSTR(xflag));
         xflag=1;
      }
   }

   if (m<2 || !xflag) {
      RTD cn2=cgd.D.norm2(); if (m>1) { cn2/=m; }
      double e=std::fabs(double(cn2-1));

      if (e>1E-12) {
         if (F) wblog(F_L,"ERR %s() "
            "got |cgc|^2 = %g (r=%d; OM=%d)",FCT,double(cn2),r,m);
         return 11;
      }
   }
   else if (xflag<2) {
      wbvector<RTD> cv2; cgd.norm2vec(r,cv2);
      RTD cn2=cv2.avg(); double e=std::fabs(double(cn2-1));

      if (cv2.len!=m) {
         if (F) wblog(FL,
            "ERR %s() got OM mismatch %d/%d",cv2.len,m);
         return 21;
      }
      if (!cv2.allEqual2(RTD(1),CG_SKIP_EPS2)) {
         if (!cv2.allEqual2(cn2,CG_SKIP_EPS2)) {
            if (F) wblog(FL,"ERR %s() got "
               "non-const normalization (%g)",FCT,double(cn2));
            return 22;
         }
         else {
            if (F) wblog(FL,"ERR %s() got "
               "non-normalized CData (%g)",FCT,double(cn2));
            return 23;
         }
      }
      if (e>1E-12) {
         if (F) wblog(FL,"ERR %s() got "
            "unexpected normalization %g (@ %g)",FCT,double(cn2),e);
         return 24;
      }
   }
   else {
      cdata__ E;
      ctrIdx ia; ia.Index(r); 
      cgd.cgsparray::contract(F_L,ia,cgd,ia, E);
      if (!E.isIdentityMatrix()) { if (F) {
         MXPut(F_L,"Ix").add(cgd,"cgd").add(E,"E");
         wblog(F_L,"ERR %s() CRef got non-orthonormal CData",FCT); }
         return 99;
      }
   }

   return 0;
};

template <class TQ, class TD>
int CData<TQ,TD>::checkQ(
   const char *F, int L, const QSet<TQ> &Q, const char *istr
 ) const {

   if ((*this)!=Q) {
      if (F || istr) wblog(F_L,
         "ERR %s() CRef inconsistency %s %s\n(%s <> %s)",
         FCT, istr? "in":"", istr? istr : "",
         STR_(this), STR(Q)
      );
      return 1;
   }

   return 0;
};

template <class TQ, class TD>
bool CData<TQ,TD>::checkAdditivityZ(const char *F, int L,
   const wbMatrix<double> &Z1, 
   const wbMatrix<double> &Z2, const wbMatrix<double> &Z,
   TD eps
 ) const {

   unsigned nz=Z1.dim2; double e=0;
   wbvector<double> z1,z2,z3,zz; 
   const SPIDX_T *I=cgd.IDX.data;

   checkConsistency(FL);
   if (cgd.numel()<=1) {
      if (!Z1.dim1 && !Z2.dim1 && !Z.dim1) return 1;
      if (!Z1.dim1 || !Z2.dim1 || !Z.dim1) wblog(FL,"ERR %s()",FCT);
   }

   if (cgd.SIZE.len<3 || cgd.SIZE.len>4) wblog(FL,
      "ERR CData::%s() inconsistency (cgd: %s)",
       FCT,cgd.sizeStr().data);

   for (SPIDX_T n=cgd.IDX.dim1, m=cgd.IDX.dim2, i=0; i<n; ++i, I+=m) {
       if (fabs(cgd.D[i])<eps) continue;

       z1.init2ref(nz,Z1.ref(I[0])); zz =z1;
       z2.init2ref(nz,Z2.ref(I[1])); zz+=z2;
       z3.init2ref(nz,Z .ref(I[2])); zz-=z3;

       if ((e=zz.norm2())>CG_EPS2) {
          MXPut(FL,"q").add(*this,"A").add(wbvector<SPIDX_T>(m,I)+1,"I")
           .add(z1,"z1").add(z2,"z2").add(z3,"z3").add(zz,"zz");
          wblog(F_L,
            "ERR z-labels not additive (%d @ %.3g => %d,%d,%d @ %.3g)",
             i+1, double(cgd.D[i]), I[0]+1, I[1]+1, I[2]+1, e
          );
       }
   }
   return 1;
};

template <class TQ, class TD>
bool CData<TQ,TD>::sameType(
   const CData<TQ,TD> &S, const char *F, int L
 ) const {

   if (this->t!=S.t) {
      if (F) wblog(F_L,"ERR %s() got different type (%s; %s)",
          FCT, qStr().data, S.qStr().data);
      return 0;
   }

   if (this->qs.len!=S.qs.len) {
      if (F) wblog(F,L,"ERR %s() "
         "got different length in qs (%d/%d)",FCT,this->qs.len,S.qs.len);
      return 0;
   }

   if (cgd.SIZE!=S.cgd.SIZE) {
      if (F) wblog(F_L,"ERR %s() got different size in cgd\n"
         "(%s; %s)",FCT, cgd.sizeStr().data, S.cgd.sizeStr().data);
      return 0;
   }

   return 1;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::initOM( const char *F, int L) {

   unsigned r=this->rank(F_L); 
   wbvector<SPIDX_T> &S=this->cgd.SIZE;

   this->cgd.wbsparray<TD>::checkSize(F_LF); 
   if (!r) wblog(FL,
      "ERR %s() got empty CGDdata (%s)",FCT,sizeStr().data);

   if (S.len==r) { S.Append(1); } else
   if (S.len!=r+1) wblog(FL,
      "ERR %s() invalid size (len=%d/%d)",FCT,S.len,r);

   return *this;
};

template <class TQ, class TD>
int CData<TQ,TD>::completeOM_DegQ(
   const char *F, int L, const wbperm &p,
   unsigned level 
) {

   unsigned dm=0, m=1, r=this->qdir.len; double e;

#ifdef QS_USING_OMP
   CG::Guard qLK, qLX;
   if (gotuser_BUF()) {
      qLK.acquire(FL,*this,"buf");
      qLX.acquire(FL,*this,NULL,"?");
   }
#endif

   if (p.len!=r || p.isIdentityPerm()) wblog(FL,"ERR %s() %s %s p=%s",
      FCT,STR_(this), p.len!=r? "invalid":"trivial", STR(p));
   if (!cgd) { return dm; }

   cdata<TD> x3;
   ctrIdx ica(p,'*'); ctrIdx icb; icb.Index(r);

   if (isRefInit()) { LoadRef(FL); }
   cgd.wbsparray<TD>::contract(FL,ica,cgd,icb,x3);

   if (level>8) {
      MXPut(FL,"Ix").add(x3,"x3").add(*this,"A").add(level,"level");
      wblog(FL,"ERR %s() reached recursive level l=%d",FCT,level);
   }

   if (x3.SIZE) { m=x3.SIZE[0];
      if (!m || x3.SIZE.len!=2)
      wblog(FL,"ERR %s() x3.SIZE = %s",FCT,SSTR(x3));
   }

   if (m<=1) {
      TD x2=x3.norm2(); e=ABS(double(x2-m));
      if (e<CG_SKIP_EPS1) {
         return dm;
      }

      cdata<TD> a(cgd); a.Permute(p);
      if ((e=double(x2))>CG_SKIP_EPS1) {
         if (x3.D.len!=1) wblog(FL,"ERR %s() D.len=%d",FCT,x3.D.len);

         a.Plus(FL,cgd,-x3.D[0]); 

         if (ABS(e=double(a.NormSignC(FL)))<1E-8) wblog(FL,
           "ERR %s() got small OM weight (cfac=%g)",FCT,e);
      }

      this->AddMultiplicity(FL,a); 
      this->cstat.t=CGD_FROM_DEC; ++dm;
   }
   else { 
      cdata<TD> a,b,E; wbvector<TD> x; TD n2;
      char update=1;

      if (r<3) wblog(FL,"ERR %s() got r=%d",FCT,r);

      for (unsigned k=0; k<m; ++k) {
         if (update) { if (k) {
            cgd.wbsparray<TD>::contract(FL,ica,cgd,icb,x3); }
            x3.wbsparray<TD>::contract(FL,"2",x3,"2*",E);
            if (E.isIdentityMatrix(CG_SKIP_EPS1)) {
               break; 
            }
         }
         n2=E.getDiag(FL,k);

         if ((e=ABS(double(n2-1)))<CG_SKIP_EPS1) { update=0; continue; }
         if (e<1E-8) wblog(FL,"ERR %s() e=%g",FCT,e);

         getMultiplicity(k,a).Permute(p);
         x3.getRow(k,x);

         cgd.wbsparray<TD>::contract(FL,r,x, b); 
         a-=b;

         if (ABS(e=double(a.NormSignC(FL)))<1E-8) {
            MXPut(FL,"Idb").add(x3,"x3").add(k,"k").add(x,"x")
            .add(E,"E").add(*this,"C").add(n2,"n2").add(k,"k").add(e,"e");
            wblog(FL,"ERR %s() got small OM weight (cfac=%g)",FCT,e);
         }
         update=1; ++dm;

         this->AddMultiplicity(FL,a); 
         this->cstat.t=CGD_FROM_DEC;
      }
   }

   if (dm) {
      dm+=completeOM_DegQ(F,L,p,level+1); 
   }

   if (!level) {
      if (dm && gotuser_BUF()) { gStore.save_CData(FL,*this); }
   }

   return dm;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::Project(const char *F, int L,
   cdata<TD> c,     
   wbvector<TD> &w  
) {
   TD x; double e=0;
   char bflag=gotuser_BUF();

   if (this->t.isAbelian()) wblog(FL,
      "WRN %s() got Abelian CData %s",FCT,STR_(this));

   if (!bflag) wblog(FL,"WRN %s() got non-BUF object",FCT);

   if (cgd.isEmpty()) {

      x=c.NormSignC(F,L); 

      if ((e=fabs(double(x)))<1E-3) {
         MXPut(FL,"a").add(c,"c").add(double(x),"x");
         wblog(FL,"WRN %s() got small new coefficient (%.3g)",FCT,e);
      }
      if (w.len) wblog(FL,
         "ERR %s() got w.len=%d for empty CData",FCT,w.len);

      w.init(1,&x);
      c.save2(cgd); 

      if (bflag) {
         gStore.rclog(this->t, PFL,
            (CG_VERBOSE>8 || (CG_VERBOSE>6 && rank(FL)<4)) && F,
            "[+] CBUF[%d] new %s",gCS.BUF.size(),STR_(this));
         gStore.save_CData(FL,*this); 
      }
      return *this;
   }

   unsigned r=rank(), m=getOM();
   int largeD=(c.D.len>(1<<26)); 

   if (r<=2) { w.init(1);
      if (!r || (r==1 && this->qs.norm2())) wblog(FL, 
         "ERR %s() got CData %s",FCT,STR_(this));
      if (m!=1 || c.sameUptoFac(cgd, w.data)!=0) wblog(FL,
         "ERR %s() invalid scalar CData\n%s",FCT,STR_(this));
      return *this;
   }

   if (r<2) wblog(FL,"ERR %s() invalid CData %s",FCT,STR_(this));
   if (c.SIZE.len!=r) wblog(FL,
      "ERR %s() got OM in input CData (%d/%d)",FCT,c.SIZE.len,r);

   checkNormSign(F_L,'x'); 

   if (m>1) { cdata__ X;
      ctrIdx ia, i1(1), im(1);

      ia.Index(r);  
      im[0]=r;      

      if (cgd.SIZE.len!=r+1) wblog(FL,
         "ERR %s() invalid CData (%d/%d)",FCT,cgd.SIZE.len,r);

      if (largeD) { Wb::MemStat(FL); } 

      for (unsigned it=0; it<2; ++it) {

         c.cgsparray::contract(F_L,ia,cgd,ia, X);
         if (!X.isVector()) wblog(FL,
            "ERR %s() got rank-%d object (%s)",FCT,X.sizeStr().data);

         cgd.cgsparray::contract(F_L,im,X,i1,c,wbperm(),-1,1);

         if (it==0) { w.init(X);
            e=double(c.norm()); if (e<CG_SKIP_EPS2) { break; }
         }
      }
   }
   else {
      for (unsigned it=0; it<2; ++it) {
         x=c.cgsparray::dotProd(F_L,cgd);

         c.Plus(FL,cgd,-x);
         if (it==0) { w.init(1,&x);
            if ((e=double(c.norm()))<CG_SKIP_EPS2) { break; }
         }
      }
   }

   if (e<1E-8) {
      if (e>CG_SKIP_EPS2) { char estr[64];
         snprintf(estr,64,"CGC ortho @ %.3g / %.3g",e,CG_SKIP_EPS2);
         if (e>1E-16) wblog(FL,"ERR %s() %s",FCT,estr);
         else wblog(FL,"WRN %s() %s\nw=[%s]",FCT,estr,STR(w));
      }
      return *this;
   }

   x=c.NormSignC(F,L); 

   if ((e=fabs(double(x)))<1E-3) wblog(FL,
      "WRN %s() got small new coefficient (%.3g)",FCT,e);
   w.Append(x);

   if (largeD) { Wb::MemStat(FL); } 

   AddMultiplicity(FL,c); 
   cstat.t=CGD_FROM_DEC;

   checkConsistency(FL);

   if (r==2) { 
      if (!cgd.isProptoId()) wblog(FL,"ERR %s() "
         "got non-scalar rank-%d CGC (%s)",FCT,r,SSTR_(this));
      else wblog(FL,"ok. %s() got scalar rank-%d CGC\n"
         "%s @ %s",FCT, r, SSTR_(this), STR_(this)
      );
   }

   if (bflag) {
      gStore.rclog(this->t, PFL, CG_VERBOSE>6,
       "[u] CBUF[%d] project #%05x %s",gCS.BUF.size(),cstat.ID,STR_(this));
      gStore.save_CData(FL,*this); 
   }

   return *this;
};

template <class TQ, class TD>
double CData<TQ,TD>::normDiff(
  const char *F, int L, const CData<TQ,TD> &S, double eps
) const {

   double e, emax=0; sameType(S,F_L);

   e=this->qs.normDiff(S.qs);
      if (e>eps) return e;
      if (e>emax) emax=e;
   e=cgd.normDiff2(S.cgd); e=std::sqrt(std::fabs(e));

   if (e>eps) return e;
   if (e>emax) emax=e;

   return emax;
};

template <class TQ, class TD>
double CData<TQ,TD>::norm2(unsigned k) const { 

   if (int(k)<0) { return double(this->cgd.D.norm2()); }

   unsigned m=getOM(FL);

   if (k>=m) wblog(FL,"ERR index out of bounds (%d/%d)",k,m);
   if (cgd.D.len!=cgd.IDX.dim1) wblog(FL, 
      "ERR %s() size mismatch (%d/%d)",FCT,cgd.IDX.dim1,cgd.D.len);

   if (m<2) { return double(cgd.D.norm2()); }
   else {
      RTD x2=0; const size_t *i4=cgd.IDX.data+3;
      for (SPIDX_T i=0; i<cgd.D.len; ++i, i4+=cgd.IDX.dim2) {
          if ((*i4)==k) x2+=(CONJ(cgd.D[i])*cgd.D[i]);
      }
      return double(x2);
   }
};

template <class TQ, class TD>
int CData<TQ,TD>::cmpOM(
   const CData<TQ,TD> &B) const { 

   unsigned r=this->qdir.len, mA=0, mB=0;

   if (!r || r!=B.qdir.len) { return -10; }
   if (r<=2) {
      if (  cgd.SIZE.len &&   cgd.SIZE.len!=r) { return -11; }
      if (B.cgd.SIZE.len && B.cgd.SIZE.len!=r) { return -12; }
   }
   else {
      if (  cgd.SIZE.len<r ||   cgd.SIZE.len>r+1) { return -13; }
      if (B.cgd.SIZE.len<r || B.cgd.SIZE.len>r+1) { return -14; }

      if (  cgd.SIZE.len>r) mA=  cgd.SIZE[r]; 
      if (B.cgd.SIZE.len>r) mB=B.cgd.SIZE[r];
   }

   return (mA==mB ? 0 : (mA<mB ? -1 : +1));
};

template <class TQ, class TD>
template <class DB>
bool CData<TQ,TD>::sameSizeR(
   const CData<TQ,DB> &B, unsigned *r_, const char *F, int L
 ) const {

   unsigned i=0, ra=this->rank(F_L), rb=B.rank(F_L);
   const SPIDX_T *sa=this->cgd.SIZE.data, *sb=B.cgd.SIZE.data;

   if (ra!=rb || ra<2 || this->t!=B.t) {
      if (F) wblog(F_L,
         "ERR %s() invalid input CG sets (%s <> %s; %s <> %s; %d/%d)",
          FCT, this->qStr().data, B.qStr().data,
          this->cgd.sizeStr().data, B.cgd.sizeStr().data, ra, rb);
      return 0;
   }

   for (i=0; i<ra; ++i) { if (sa[i]!=sb[i]) return 0; }

   if (r_) { (*r_)=ra; }

   return 1;
};

template <class TQ, class TD>
template <class T2>
int CData<TQ,TD>::sameSizeR(
   const char *F, int L, const wbvector<T2> &S
 ) const {

   if (cgd.isEmpty()) {
      if (isAbelian() || !S.len) { return 1; }
      return 0;
   }

   unsigned i=0, ra=this->rank(FL), m=gotOM(FL), m2=1;
   const SPIDX_T *sa=this->cgd.SIZE.data;

   if (S.len==ra+1) { m2=S[ra]; } else
   if (!ra || S.len<ra || S.len>ra+1) { 
      if (F) wblog(F,L,
         "ERR %s() CData rank mismatch (%s: %s <> %s; %d/%d)",
          FCT, this->qStr().data,
          this->cgd.sizeStr().data, STR(S), ra, S.len);
      return -1;
   }

   if (sa) {
      for (; i<ra; ++i) { if (S.data[i]!=sa[i]) { return 0; }}
   }
   else { 
      SPIDX_T s=cgd.dim();
      if (ra!=2) wblog(FL,
         "ERR %s() invalid CData\n%s",FCT,STR_(this));
      for (; i<ra; ++i) { if (S.data[i]!=s) { return 0; }}
   }

   if (!m) { m=1; }

   if (m2>m) { return 2; } else
   if (m2<m) { return 3; } else { return 1; }
};

template <class TQ, class TD>
template <class DB>
bool CData<TQ,TD>::sameSizeR(
   const cdata<DB> &b, unsigned *r_, const char *F, int L
 ) const {

   unsigned i=0, ra=this->rank(F_L), rb=b.rank(F_L);
   const SPIDX_T *sa=this->cgd.SIZE.data, *sb=b.SIZE.data;

   if (ra!=rb || ra<2) {
      if (F) wblog(F_L,
         "ERR %s() invalid input CG sets (%s: %s <> %s; %d/%d)",
          FCT, this->qStr().data,
          this->cgd.sizeStr().data, b.sizeStr().data, ra, rb);
      return 0;
   }

   for (i=0; i<ra; ++i) { if (sa[i]!=sb[i]) return 0; }

   if (r_) { (*r_)=ra; }

   return 1;
};

template <class TQ, class TD>
bool CData<TQ,TD>::isSymmmetric() const {

   unsigned r=this->qdir.len;

   if (r%2) { return 0; }
   if (!r) {
      if (cgd.SIZE.len || cgd.IDX.dim1 || cgd.D.len>1)
         wblog(FL,"ERR %s() got invalid empty CGC",FCT);
      return 1;
   }

   unsigned i=0, r2=r/2;
   const char *qd=this->qdir.data;

   for (; i<r2; ++i) { if (qd[i]<=0 || qd[i+r2]>=0) return 0; }

   r2=this->qs.len/2;
   if (memcmp(this->qs.data, this->qs.data+r2, r2*sizeof(TQ))) return 0;

   return 1;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::AddMultiplicity(
   const char *F, int L, const cdata<TD> &c
){
   SPIDX_T i,l,m; unsigned n, r=0;
   wbsparray<TD> X;
   const cdata<TD> &a=this->cgd;

   a.checkSize(F_LF,"A");
   c.checkSize(F_LF,"C");

   if (!sameSizeR(c,&r)) wblog(F_L,
      "ERR %s() got incompatible CData\n%s: %s <> %s",FCT,
      STR_(this), a.sizeStr().data, c.sizeStr().data);
   if (r<3) wblog(FL,"ERR %s() unexpected rank-%d CData",FCT,r);

   wbvector<SPIDX_T> s(r+1);

   for (i=0; i<a.SIZE.len; ++i) s[i]=a.SIZE[i];
   for (; i<s.len; ++i) s[i]=1;
   m=(s[r]++); 

   X.init(s, a.IDX.dim1 + c.IDX.dim1);

   MEM_CPY<TD>(X.D.data,         a.D.len, a.D.data);
   MEM_CPY<TD>(X.D.data+a.D.len, c.D.len, c.D.data);

   for (n=a.IDX.dim2, i=0; i<a.IDX.dim1; ++i) {
      X.IDX.recSetP(i,a.IDX.rec(i),n);
   }
   for (n=c.IDX.dim2, l=i, i=0; i<c.IDX.dim1; ++i, ++l) {
      X.IDX.recSetP(l,c.IDX.rec(i),n);
      X.IDX(l,r)=m;
   }

#ifndef WB_SKIP_ASSERT
   { ctrIdx icm; icm.Index(r);
     wbsparray<TD> E; X.contract(FL,icm,X,icm,E);
     TD nrm=TD(1); 
     if (!E.isProptoId(nrm, CG_EPS1)) {
        MXPut(FL,"a").add(*this,"A")
        .add(c,"c").add(X,"X").add(E,"E").add(icm,"ic");
        wblog(FL,"ERR %s() got non-orthonormal OM space\n%s\n%s",
        FCT,STR_(this),STR2(cstat,'V'));
     }
   }
#endif

   X.save2(this->cgd);
   cstat.update_m(); 

   gStore.rclog(this->t, PFL, CG_VERBOSE>6, 
   " +  %s() #%05x %s",FCT,cstat.ID,STR_(this));

   return *this;
};

template <class TQ, class TD> 
cdata<TD>& CData<TQ,TD>::getMultiplicity(unsigned im, cdata<TD> &c) const {

   unsigned m=gotOM(FL), r=rank(FL);

   if (im && im>=m) wblog(FL,
      "ERR %s() OM index out of bounds (%d/%d)",FCT,im,m);
   if (r<3) wblog(FL,"ERR %s() got r=%d",FCT,r);
   cgd.checkSize(FLF,"cgd");

   if (m>1) { cgd.splitLast(im,c); }
   else {
      c.init(cgd); if (m) {
         if (c.SIZE[r]!=1) wblog(FL,"ERR %s() %s /%d",FCT,SSTR(c),r);
         c.SIZE.len=r;
      }
   }

   if (!c.D) wblog(FL,
      "ERR %s() got empty OM component l=%d/%d\n%s -> S=%s",
      FCT,im,m,STR_(this),SSTR(cgd));
   return c;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::save2(CData<TQ,TD> &B) {

   if (&B==this) { return B; }
   int cflag=0;

   if (!B.gotuser_BUF() || !B.qdir.len || !B.qs.len) {
      B.t=this->t;
      this->qs.swap(B.qs);
      this->qdir.swap(B.qdir);

      cgd.swap(B.cgd); 
      cflag=-1; 
   }
   else {
      if ((QSet<TQ>&)B != *this) wblog(FL,
         "ERR %s() got change in BUF data\n%s -> %s",
         FCT,STR(B),STR_(this)
      );

      if (this->t.isAbelian()) { 
         if (cgd!=B.cgd) wblog(FL,
            "ERR %s() BUF change in abelian cdata %s\n%s -> %s",
            FCT, STR_(this), SSTR(B.cgd), SSTR(cgd));
         cflag=-1; 
      }
   }

   if (cflag>=0) { 
      unsigned r=this->qdir.len;
      int q=cmpOM(B);

      if (!r || q<-1) { wblog(FL,
         "ERR %s() BUF unexpected CData (r=%d, q=%d)\n%s / %s",
         FCT,r,q, STR_(this), SSTR(cgd));
      }
      if (q<0) { wblog(FL,
         "ERR %s() BUF reduces rank/OM (q=%d, r=%d)\n%s -> %s",
         FCT,q,r, SSTR(B.cgd), SSTR(cgd));
      }

      if (cgd.SIZE.len==B.cgd.SIZE.len) { if (q) { 
         SWAP(cgd.SIZE[r],B.cgd.SIZE[r]); }
         cgd.IDX.swap(B.cgd.IDX);
         cgd.D  .swap(B.cgd.D);

         cflag=(cstat!=CGD_REF_INIT && B.cstat.t!=CGD_REF_INIT);
      }
      else if (cgd.isDiag()) { cgd.diag2reg(); 
         if (cgd.SIZE!=B.cgd.SIZE) { wblog(FL,
            "ERR %s() got SIZE mismatch %s / %s (r=%d)",
             FCT, SSTR(B.cgd), SSTR(cgd), r);
         }
         cgd.IDX.swap(B.cgd.IDX);
         cgd.D  .swap(B.cgd.D); cflag=2;
      }
      else { cgd.swap(B.cgd); }
   }

   B.cstat=cstat; 

   if (cflag>0 && cgd!=B.cgd) {
      unsigned r=this->qdir.len; double e=0;
      ctrIdx ic; ic.Index(r);
      cdata<TD> X,E; TD x;

      B.cgd.cgsparray::contract(FL,ic,cgd,ic,X);

      if (X.SIZE.len==2) { unsigned m=X.SIZE[1];
         if (!m || m>X.SIZE[0]) wblog(FL,"ERR %s() "
            "got reduced OM in BUF\n%s -> %s",FCT,SSTR(B.cgd),SSTR(cgd));
         X.getBlock(FL,0,m-1,0,m-1,E);

         if (!E.isProptoId(x,CG_EPS1)) {
            MXPut(FL,"Idb").add(cgd,"A").add(B.cgd,"B").add(X,"X").add(E,"E");
            wblog(FL,"ERR %s() got altered existing OM space",FCT);
         }
         e=ABS(double(x-1));
      }
      else {
         if (X.SIZE.len>2 || !X.D) wblog(FL,
            "ERR %s() X=%s (n=%d)",FCT,SSTR(X),X.D.len);
         e=ABS(double(X.D[0]-1)); 
      }
      if (e>1E-20) wblog(FL,
         "ERR %s() got CG difference @ %.3g",FCT,sqrt(e));
   }

   this->init();  
   return B;
};

template <class TQ, class TD>
bool CData<TQ,TD>::isBasicCG(const char *F, int L) const {

   if ((this->rank(F_L))!=3 || this->qdir!="++-") {
      if (F) wblog(F_L,
         "ERR invalid basic CG set (%s; %s)",
         this->sizeStr().data, STR(this->qdir));
      return 0;
   }
   return 1;
};

template <class TQ, class TD>
bool CData<TQ,TD>::isScalar(char dflag) const {

   if (!cgd.SIZE.len) {
      if (cgd.D.len==1 || (!cgd.D.len && CGD_ABELIAN)) return 1;
   }
   else if (cgd.isScalar()) { return 1; }

   if (cgd.D.len<1) {
      if (!isRefInit()) wblog(FL,
         "WRN %s() got empty CRef data\n%s",FCT,STR_(this));
      return 0;  
   }

   if (this->qdir.len!=2 || this->qdir[0]==this->qdir[1]) return 0;

   if (dflag) { double x;
      if (isRefInit()) { 
         if (cgd.D.len!=1) wblog(FL,
        "ERR %s() invalid CRef data\n%s",FCT,STR_(this));
      }
      else { RTD q; if (!cgd.isProptoId(q,CG_EPS1))
         wblog(FL,"ERR %s() got invalid scalar CData\n%s\ncdata: %s [%d]",
         FCT, STR_(this), cgd.info2Str().data, cgd.isDiag(FL));
      }
      if (fabs(x=double(cgd.D[0]))<CG_EPS1) wblog(FL,
         "WRN %s() got small scalar CData (%.4g)",FCT,x);
      return 1;
   }
   return 0;
};

template <class TQ, class TD>
SPIDX_T CData<TQ,TD>::getBasicCGSize(unsigned *m) const {

   unsigned r=rank(FL); SPIDX_T n=this->cgd.SIZE.len;
   const SPIDX_T *s=this->cgd.SIZE.data;

   if (n==r) { if (m) (*m)=1; }
   else {
      if (n!=r+1) wblog(FL,"ERR %s() got %d/%d",FCT,n,r+1);
      if (m) { (*m)=s[r]; }
   }

   if (r) { n=s[0]; for (unsigned i=1; i<r; ++i) n*=s[i]; }
   else { n=0; }

   return n;
};

template <class TQ, class TD>
wbsparray<TD>& CData<TQ,TD>::getBasicCG(
   unsigned k, wbsparray<TD> &a
 ) const {

   unsigned om=this->getOM(FL), r=this->rank(FL);

   if (k>=om) wblog(FL,"ERR OM out of bounds (k=%d/%d)",k+1,om);
   if (this->cgd.SIZE.len==r) {
      if (om!=1) wblog(FL,"ERR cdata inconsistency (%d)",om);
      a.init(this->cgd);
   }
   else {
      const cdata<TD> &c = this->cgd;
      wbvector<size_t> S(r,c.SIZE.data);

      const size_t *im=c.IDX.data+r;
      SPIDX_T i=0, l=0, N=c.D.len; unsigned m=c.IDX.dim2;

      for (; i<N; ++i, im+=m) { if ((*im)==k) { ++l; }}
      a.init(S,l);

      im=c.IDX.data+r;
      for (l=i=0; i<N; ++i, im+=m) { if ((*im)==k) {
         a.setRecP(l++, im-r, c.D[i]);
      }}
   }
   return a;
};

template <class TQ, class TD>
int CData<TQ,TD>::getCG_set(
   const char *F, int L, wbIndex &Idx, wbsparray<TD> &a
 ) const {

   const wbvector<SPIDX_T> &S2=this->cgd.SIZE;
   unsigned nq=this->t.qlen(), r=this->qdir.len, r2=S2.len, rm=r2-r;

   if (!r2) { r2=this->cgd.rank(); rm=r2-r; } 

   if (!r || !r2 || r2<r || !nq || this->qs.len/nq!=r || this->qs.len%nq)
      wblog(F_L,"ERR %s() got invalid/inconsistent QSet\n{%d,%d,%d}: %s",
      FCT, r,r2,nq, ((QSet<TQ>*)this)->toStr().data
   );

   if (Idx.isEmpty()) {
      if (r2>r) { Idx.init(rm,S2.data+r); }
      else { a.init(cgd);
         if (cgd.D.len) {
            widx_t s=1; Idx.init(1,&s); ++Idx;
            return 1;
         }
         else return 0;
      }
   }
   else {
      if ((r2!=r && r2!=r+Idx.len) || (r2==r && Idx.len!=1)) wblog(F_L,
         "ERR %s() got rank mismatch in OM (%d+%d <> %d)",
         FCT,r,Idx.len,r2
      );
   }

   if (!(++Idx)) { a.init(); return 0; }
   Idx.checkValid(FL); 

   const wbMatrix<SPIDX_T> &IDX=this->cgd.IDX;
   const wbvector<TD> &D=this->cgd.D;
   wbvector<SPIDX_T> S(r,S2.data);

   if (IDX.dim1!=D.len || IDX.dim2!=S2.len) wblog(FL,
      "ERR %s() cdata mismatch (%dx%d <> %dx%d)",
      FCT,IDX.dim1,IDX.dim2,D.len,S2.len);
   if (!IDX.dim1 || !IDX.dim2) wblog(FL,
      "ERR %s() got empty cdata (%dx%d <> %dx%d)",
      FCT,IDX.dim1,IDX.dim2,D.len,S2.len);

   widx_t
      i1=Wb::findfirst_sorted(NULL,0,
         Idx.data, IDX.data+r, rm, IDX.dim1,IDX.dim2,-1),
      i2=Wb::findlast_sorted(NULL,0,
         Idx.data, IDX.data+r, rm, IDX.dim1,IDX.dim2,-1),
      i=0, l=i2-i1+1;

   if (i1>i2) wblog(FL,
      "ERR %s() got unsorted data (%d,%d)",FCT,i1,i2);

   if (i2>=IDX.dim1) {
      if (i1<IDX.dim1) wblog(FL,"ERR %s() i1=%d, i2=%d",FCT,i1,i2);
      a.init(S,0); return 1;
   }

   a.init(S,l);
   for (; i<l; ++i, ++i1) {
      a.setRecP(i, IDX.rec(i1), D[i1]);
   }

   return 1;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::reduceto1JSymbol(
   const char *F, int L, CData<TQ,TD> &C) const {

   unsigned i=0, m=this->qs.len/3, l=2*m;
   const cdata__ &c=this->cgd;

   if (!this->isStd3()) wblog(F_L,"ERR %s() "
      "got invalid CData %s",FCT, ((QSet<TQ>*)this)->toStr().data);
   if (c.SIZE.len!=3 || c.SIZE[2]!=1) wblog(FL,"ERR %s() "
      "incompatible cgdata (%s)",FCT, c.sizeStr().data);
   if (c.IDX.dim2!=3) wblog(FL,"ERR %s() invalid cgdata (%s; %s)",
      FCT, c.IDX.sizeStr().data, c.sizeStr().data);
   for (; i<m; ++i) { if (this->qs.data[l+i]) {
      wblog(FL,"ERR %s() incompatible CData\n%s",
      FCT, ((QSet<TQ>*)this)->toStr().data); }
   }

   C=(const QSet<TQ>&)(*this); {
      C.cstat.init(CGD_1JSY_ST3);
      C.qs.len=l;
      C.qdir.len=2;
      cgd.skipTrailingSingletons(FL,C.cgd,2);
   }
   return C;
};

template <class TQ, class TD>
CData<TQ,TD>& CData<TQ,TD>::Conj(unsigned k) {

   unsigned m=this->t.qlen(); ctrIdx ia,ib;
   cdata<TD> c;

   if (this->qdir.len*m!=this->qs.len) wblog(FL,"ERR %s "
      "having %d*%d =? %d",FCT,STR_(this),this->qdir.len,m,this->qs.len);
   if (k>=this->qdir.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k,this->qdir.len);

   wbperm P;
   CRef<TQ> R;
   gCS.getIdentity1J(FL,R,this->t,this->qs.data+k*m,-1,LB_GEN); 

   ia.init1(k);
   ib.init1(R.cgp.isEmpty() ? 0 : 1, this->qdir[k]>0 ? '*' : 0);

   cgd.contract(FL,ia,R.cgb->cgd,ib,c,NULL,NULL,0);  
   { double e=R.cgw[0]; unsigned d=R.cgb->cgd.dim(); e=fabs(e*e-d);
     if (e>1E-12) wblog(FL,
        "ERR %s() unexpected cgw^2 = %g / %d @ %.3g",FCT,R.cgw[0],d,e);
     c*=SQRT(RTD(d));
   }

   P.initLastTo(k, cgd.SIZE.len ? cgd.SIZE.len : this->qdir.len); 
   c.permute(this->cgd,P); 

   QSet<TQ>::Conj(k);
   return *this;
};

template <class TQ, class TD>
int CData<TQ,TD>::gotDiff(const CData &B, char lflag) const {

   if (lflag<=0) { return ((*this)==B ? 0 : 1); }
   if (this->t!=B.t || this->qs!=B.qs || this->qdir!=B.qdir) { return 2; }

   if (lflag>2) {
      if (lflag=='l') lflag=1; else 
      if (lflag=='L') lflag=2; else
      wblog(FL,"WRN %s() unexpected lflag=%d",FCT,lflag);
   }

   if (lflag<2) { int got=0;
      if (!((cstat==CGD_REF_INIT) ^ (B.cstat.t==CGD_REF_INIT))) {
         if (!checkOM(FL)) { if (cgd!=B.cgd) { got|=4; }}
         if (cstat.t!=B.cstat.t) { got|=8; }
      }
      if (got) return got;
   }

   return (cstat.sameAs(B.cstat,2) ? 0 : 16); 
};

template <class TQ, class TD>
void CData<TQ,TD>::info(const char *istr, const char *F, int L) const {

   unsigned e=0, r=rank(F_L); 

   checkConsistency(F_L); 

   if (F) wblog(F,L,
      "CG Space %s :: %s",istr && istr[0] ? istr : "",qStr().data);
   else printf(
      "\nCG Space %s :: %s\n",istr && istr[0] ? istr : "",qStr().data);
   printf("\n");

   printf("  Q=[%s] (rank=%d)\n",this->QStr().data,r);

   if (cgd.SIZE.len)  cgd.info();

   if (e) wblog(FL,"ERR (%d)",e);
   printf("\n");
};

template <class TQ, class TD>
wbstring CData<TQ,TD>::toStr(char vflag) const { 

   wbstring s_; 
   if (isEmpty()) { return (s_="(empty)"); }

   s_.init(128);
   unsigned n=s_.len-1; char *s=s_.data;

   unsigned l=0, m=1; 
   m=cgd.SIZE.prod(this->qdir.len, -1, 1); 

   if (cstat.t>=CGD_NUM_TYPES) wblog(FL,"ERR %s() "
      "ctype out of range (%d/%d)",FCT,cstat.t,CGD_NUM_TYPES
   );

   l=snprintf(s,n,"%s",STR_((const QSet<TQ>*)this));
   if (!vflag) { if (m>1 && l<n) l+=snprintf(s+l,n-l," @ %d",m); }
   else {
      if (cstat.t!=CGD_UNKNOWN && l<n) l+=snprintf(s+l,n-l,
         " %s",CGD_TYPE_STR[cstat.t]);
      if (l<n) 
         l+=snprintf(s+l,n-l,"%s",cstat.u2Str(" ").data);
   }

   if (vflag&4 && l<n) {
      for (; l<32; ++l) { s[l]=' '; } 
      l+=snprintf(s+l,n-l," %s",SSTR(cgd));
   }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s_;
};

template <class TQ, class TD>
wbstring CData<TQ,TD>::sizeStr(char xflag) const {

   unsigned r=rank(FL);

   if (!r) {
      if (cgd.SIZE.len) wblog(FL,
         "ERR %s() got r=%d/%d",FCT,r,cgd.SIZE.len);
      return "(empty)";
   }

   unsigned l=0, n=128; char s[n];

   if (xflag && r>2) { 
      if (cgd.SIZE.len==r+1) {
         l+=snprintf(s+l,n-l,"%s", cgd.sizeStr().data);
      }
      else if (cgd.SIZE.len==r) {
         wbvector<SPIDX_T> s1(1); s1[0]=1; 
         wbvector<SPIDX_T> S(cgd.SIZE,s1);
         l+=snprintf(s+l,n-l,"%s", SSTR(S));
      }
      else wblog(FL,"ERR %s() invalid rank %d/%d",FCT,cgd.SIZE.len,r);
   }
   else {
      if (cgd.SIZE.len!=r+1) {
         l+=snprintf(s+l,n-l,"%s", cgd.sizeStr().data);
      }
      else {
         wbvector<SPIDX_T> S(cgd.SIZE); --S.len;
         l+=snprintf(s+l,n-l,"%s",SSTR(S)); if (l<n) {
         l+=snprintf(s+l,n-l," @%ld", cgd.SIZE[r]); }
      }
   }
   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};

template <class TQ> inline
CRef<TQ>& CRef<TQ>::initBase(
   const CDATA_TQ* cgr_, unsigned m,
   const CDATA_TQ* x 
) {

   if (x) {
      if (!cgr_) wblog(FL,"ERR %s() got NULL cgr for %s",FCT,SSTR_(x));
      if (!cgr_->isEmpty() && (*cgr_) != (const QSet<TQ>)(*x))
      wblog(FL,"ERR %s() QSet mismatch\n%s <> %s",FCT,STR_(cgr_),STR_(x));
   }
   else { x=cgr_; }

   if (x && x->t.isNonAbelian()) {
      cgb=(CDATA_TQ*)cgr_; 

      if (int(m)<0) { m=cgb->getOM(FL); } else {
      if (!cgb->isEmpty() && m>(cgb->getOM(FL))) wblog(FL,"ERR %s() "
         "OM out of bounds (%d/%d) %s",FCT,m,cgb->getOM(FL),STR_(cgb));
      }
      cgw.initIdentity(m);
   }
   else {
      cgb=NULL; if (x) x->checkQ_abelian(FL);
      if (m!=1) wblog(FL,"ERR %s() invalid abelian cref (%d)",FCT,m);
      cgw.init();
   }

   cgp.init(); conj=0; rtype=0;

   return *this;
};

template <class TQ> inline
CRef<TQ>& CRef<TQ>::initBase(const CDATA_TQ &C, unsigned m) {
   QSet<TQ> Q(C);
   const CDATA_TQ &Cb=gCS.getBUF(FL,Q,LB_REF); 
   return initBase(&Cb,m);
};

template <class TQ> inline 
unsigned CRef<TQ>::qdim(unsigned k,
   const QType *qt 
 ) const {

   if (qt) {
      if (isAbelian()) { if (!qt->isAbelian()) wblog(FL,
      "ERR %s() abelian symmetry mismatch (%s)",FCT,STR_(qt)); }
      else if ((*qt)!=cgb->t) wblog(FL,
      "ERR %s() symmetry mismatch (%s/%s)",FCT,STR(cgb->t),STR_(qt));
   }

   return Size(k); 
};

template <class TQ>
CRef<TQ>& CRef<TQ>::initIdentityR(
   const char *F, int L, const QType &t, const TQ *qs,
   unsigned dim, 
   char xflag    
){
   char isa=t.isAbelian(); cgp.init(); conj=0;

   if (isa && !xflag) {
      cgw.init(); cgb=NULL; rtype=CGR_ABELIAN;
   }
   else {
      cgb=&gCS.getIdentityC(F_L,t,qs,dim);
      cgw.init(1,1);

      if (!cgb->cgd.D.len) {
         if (!isa) wblog(FL,
            "ERR %s() got invalid CData\n%s",FCT,cgb->toStr().data);
         cgw[0]=1;
      }
      else {
         cgw[0]=double(1/cgb->cgd.D.data[0]); 
      }
   }

   return *this;
};

template <class TQ>
CRef<TQ>& CRef<TQ>::initIdentity1J(
   const char *F, int L, const QType &t, const TQ *qs, unsigned dim,
   char xflag
){
   cgp.init(); conj=0;

   if (t.isAbelian() && !xflag)
        { cgw.init(); cgb=NULL; rtype=CGR_ABELIAN; }
   else { gCS.getIdentity1J(F_L,*this,t,qs,dim); }

   return *this;
};

template <class TQ>
CRef<TQ>& CRef<TQ>::Reduce2Identity(char xflag) {

   if (!cgb) {
      if (cgp || cgw) { 
         wblog(FL,"ERR %s() invalid CRef\n%s",FCT,STR_(this)); }
      return *this;
   }

   unsigned r=rank(FL);

   if (r<=2) { 
      if (r<2)
           { wblog(FL,"ERR %s() got rank-%d CRef",FCT,r); }
      else { wblog(FL,"WRN %s() got rank-%d CRef",FCT,r); }
      return *this;
   }
   if (isScalar()) {
      if (!cgb->cgd.SIZE.len && cgb->cstat!=CGD_ABELIAN) wblog(FL,
         "ERR %s() invalid CData\n%s",FCT,STR_(cgb));
   }

   double cfac=1.;
   unsigned d=1, n=cgb->t.qlen(), k=(cgp.len ? cgp[0] : 0);
   QSet<TQ> Q(*this);

   if (!Q.isScalar() || !wscalar()) {
      MXPut(FL,"a").add(*this,"R").add(*cgb,"C").add(cgb->cgd,"c");
      wblog(FL,"ERR %s() unexpected scalar CRef\n%s",FCT,STR_(this));
   }

   if (!cgb->cgd.SIZE.len) {
      cfac=double(cgb->getScalar());
   }
   else { RTD x=1;
      d=cgb->cgd.SIZE.el(k);
      if (cgb->cgd.isIdentity(cgp,&x)!=0) {
         MXPut(FL,"a").add(*cgb,"cgr").add(cgb->cgd,"cgd");
         wblog(FL,"ERR %s() failed to reduce to Id (%g)",
         FCT,cgb->cgd.D.len? double(cgb->cgd.D[0]) : -1.);
      }
      cfac=double(x);
   }

   cfac*=cgw[0]; 

   QType t=cgb->t;
   qset<TQ> qs(n, cgb->qs.data+k*n);

   initIdentityR(FL,t,qs.data, d, xflag); 
   cgw*=cfac;

   return *this;
};

template <class TQ, class TD>
wbvector<RTD>& CData<TQ,TD>::trace(
   const char *F, int L, wbvector<RTD> &cgt) const {

   if (cgd.isEmpty() && cstat==CGD_ABELIAN) {
      TD x=1; return cgt.init(1,&x);
   }

   unsigned r=rank();

   if (cgd.SIZE.len==r || (!cgd.SIZE.len && r==2)) { 
      TD x=cgd.trace(); 
      return cgt.init(1,&x);
   }
   else {
      if (cgd.SIZE.len!=r+1) wblog(F_L,"ERR %s() "
         "invalid OM data (%d/%d)",FCT,cgd.SIZE.len, r);
      return cgd.trace(r,cgt);
   }
};

template <class TQ>
bool CRef<TQ>::isSymmmetric(const char *F, int L) const {

   if (!cgb) wblog(F_L,"ERR %s() got %s",FCT,STR_(this));

   if (!cgb->isSymmmetric() || cgp.len%2) { return 0; }

   if (cgp.len) {
      unsigned i=0, r2=cgp.len/2;
      if (cgp[0]<cgp[r2])
           { for (; i<r2; ++i) { if (cgp[i]+r2!=cgp[r2+i]) return 0; }}
      else { for (; i<r2; ++i) { if (cgp[i]!=cgp[r2+i]+r2) return 0; }}
   }

   return 1;
};

template <class TQ>
bool CRef<TQ>::isAbelian(const char *F, int L) const {

   if (cgb) 
        { if (!cgb->isAbelian() ) { return 0; }}
   else { if (rtype!=CGR_ABELIAN) { return 0; }}

   if (cgp.len) {
      if (cgp.len!=cgb->qdir.len) wblog(F_L,
         "ERR %s() invalid ablian CRef (cgp=[%s])",FCT,STR(cgp));
      if (!cgb || F) wblog(F_L,
         "WRN %s() CRef with p=%s",FCT,STR(cgp));
   }
   if (cgw) {
      if (!wscalar1()) wblog(F_L,
         "ERR %s() invalid ablian CRef (cgw=[%s])",FCT,STR(cgw));
      if (!cgb || F) wblog(F_L, 
         "WRN %s() CRef with cgw=%s",FCT,STR(cgw));
   }
   return 1;
};

template <class TQ>
bool CRef<TQ>::isScalar(char dflag) const {

   if (!cgb) { unsigned n=wnumel();
      if (rtype<CGR_ABELIAN ||  
         (rtype!=CGR_CTR_SCALAR && rtype!=CGR_CTR_ZERO)) { wblog(FL,
         "ERR %s() got unexpected rtype=%d (cgw %s)",FCT,rtype.t,SSTR(cgw));
      }
      if (cgp || conj || n>(rtype==CGR_CTR_SCALAR ? 1:0)) {
         char s[32]; if (n==1)
              { snprintf(s,32,"cgw=%g, %s",cgw[0],SSTR(cgw)); }
         else { snprintf(s,32,"|cgw|=%g, %s",cgw.norm(),SSTR(cgw)); }
         wblog(FL,"ERR %s() got invalid scalar (%s)",FCT,s);
      }
      if (n==1 && cgw[0]!=1) wblog(FL, 
         "ERR %s() unexpected cgw=%g for cgb=0",FCT,cgw[0]);
      return 1;
   }
   else return cgb->isScalar(dflag); 
};

template <class TQ>
int CRef<TQ>::checkQ( const char *F, int L,
   const QSet<TQ> &Q 
 ) const {

   if (!cgb) {
      if (!Q.t.isAbelian()) {
         if (F) wblog(FL,"ERR %s() missing cgb data (%s)",FCT,STR(Q.t));
         return 1;
      }
   }
   else if (!sameQSet(Q)) {
      if (F) wblog(FL,
         "ERR %s() QSet mismatch\n   %s\n<> %s",FCT,STR_(this),STR(Q));
      return 2;
   }
   return 0;
};

template <class TQ>
unsigned CRef<TQ>::numel() const {

   if (!cgb) {
      if (rtype==CGR_ABELIAN || rtype==CGR_CTR_SCALAR) return 1;
      else {
         wblog(FL,"WRN %s() got rtype=%d",FCT,rtype.t);
         return 0;
      }
   }
   return cgb->cgd.numel();
};

template <class TQ>
unsigned CRef<TQ>::Size(unsigned k, unsigned r) const {

   if (!cgb || cgb->cgd.isScalar()) return 1;

   const cdata__ &c=cgb->cgd;
   const wbvector<unsigned> &S=c.SIZE;

   if (!S.len && c.isDiag()) {
      if (k>2 || (int(r)>=0 && r!=2)) wblog(FL,
         "ERR %s() invalid rank-%d object (k=%d)",FCT,r,k);
      if (cgp.len && cgp.len!=2) wblog(FL,
         "ERR %s() invalid cgp.len=%d/2",FCT,cgp.len);
      return c.D.len; 
   }

   if (int(r)<0) r=cgb->qdir.len;
   if (k>r) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,r);

   if (r<2 || (S.len!=r && S.len!=r+1)) wblog(FL,
      "ERR %s() rank mismatch (%s; %d/%d/%d; d=%d)",
       FCT,c.sizeStr().data, k+1,r,S.len, c.D.len
   );

   if (k<r && cgp.len) {
      if (cgp.len!=r) wblog(FL,"ERR %s() "
         "invalid cgp=(%s/%d)",FCT,STR(cgp),r);
      k=cgp.data[k];
   }

   return (k<S.len ? S.data[k] : 1);
};

template <class TQ>
unsigned CRef<TQ>::getOM(const char *F, int L, char wflag) const {
   unsigned m=0; 

   if (cgb) {
      if (wflag==1) { m=wdim1(F_L); } else
      if (wflag==2) { m=wdim2(F_L); }
      else {
         m=cgb->getOM(); 
         if (!cgw || m<wdim1()) wblog(F_L,     
            "ERR %s() invalid cgw (%d / %s)",FCT,m,SSTR(cgw));
         if (wflag) wblog(FL,"WRN %s() got wflag=%s",FCT,cSTR(wflag));
      }
   }
   else if (isAbelian()) { m=1; }
   else { wblog(F_L,"ERR %s() unknown OM [%s]",FCT,STR_(this)); }

   return m;
};

template <class TQ>
unsigned CRef<TQ>::rankS(char lflag) const { 
   if (cgb) {
      unsigned r=cgb->rankS();
      if (!lflag) { unsigned m=cgb->getOM();
         if (!cgw || wdim1()>m) wblog(FL,
            "ERR %s() invalid cgw (%s /%d)",FCT,SSTR(cgw),m);
         if (cgp.len && cgp.len!=r) wblog(FL,
            "ERR %s() invalid cgp (len=%d/%d)",FCT,cgp.len,r
         );
      }
      return r;
   }
   else {
      if (!lflag) {
         if (cgw || cgp) wblog(FL,"ERR %s() "
         "cgw or cgp out of bounds (%s, %d/0)",FCT,SSTR(cgw),cgp.len);
      }
      return 0;
   }
};

template <class TQ>
unsigned CRef<TQ>::rank(const char *F, int L, char lflag) const {
   if (cgb) {
      unsigned r=cgb->rank();
      if (!lflag) { unsigned m=cgb->getOM();
         if (!cgw || wdim1()>m) wblog(F_L,
            "ERR %s() invalid cgw (%s /%d)",FCT,SSTR(cgw),m);
         if (cgp.len && cgp.len!=r && cgp.len!=(r-1)) wblog(F_L,
            "ERR %s() invalid cgp (len=%d/%d)",FCT,cgp.len,r
         );
      }
      return r;
   }
   else {
      wblog(F_L,"ERR %s() unknown rank (since cref=NULL)",FCT);
      return 0;
   }
};

template <class TQ>
template <class T>
wbvector<T>& CRef<TQ>::getSize(wbvector<T> &S, char bare) const {

   if (cgb) {
      cgb->getSize(S, bare, cgp.len ? &cgp : NULL);
   }
   else {
      if (cgp || cgw) wblog(FL,"ERR %s() "
         "invalid scalar/empty CRef (%s; %d)",FCT,SSTR(cgw),cgp.len);
      S.init();
   }

   return S;
};

template <class TQ>
int CRef<TQ>::HConjOpScalar() {

   if (cgp.len) { wblog(FL,
      "WRN %s() expecting scalar operator\ngot %s with cgp=[%s]",
         STR_(this),FCT,STR(cgp));
      cgp.init();
   }

   if (!cgb || cgb->isEmpty()) {
      if (conj) {
         wblog(FL,"WRN %s() got conj=%d for scalar",FCT,conj);
         conj=0;
      }
      return 0; 
   }

   conj=0; 

   unsigned i=0, r=cgb->qdir.len, n=(r ? cgb->qs.len/r : 0);

   if (r<2 || cgb->qs.len%r || !n) wblog(FL,
      "ERR %s() got rank-%d QSpace (qlen=%d)",FCT,r,n);

   if (r!=3) {
      if (r!=2) {
         wblog(FL,"WRN %s() got rank-%d CGR data",FCT,r);
         return 1;
      }
      return 0; 
   }

   const char *d=cgb->qdir.data;
   const TQ *q=cgb->qs.data;

   if (d[0]!=-d[1]) { return 2; }

   if (memcmp(q,q+n,n)) { return 12; } 

   for (q+=2*n; i<n; ++i) { if (q[i]) { return (30+i); }}

   return 0;
};

template <class TQ>
int CRef<TQ>::checkAbelian(const char *F, int L) const {

   if (cgb) {
      if (!cgb->isAbelian()) {
         if (F) wblog(F,L,"ERR %s() "
            "got cref with Abelian symmetry (%s)",FCT,cgb->toStr().data);
         return 1;
      }
      if (cgp.len && cgp.len!=cgb->qdir.len) wblog(FL,
         "ERR %s() got invalid abelian CRef\n%s",FCT,STR_(this));
   }

   int e=0;

   if (cgp.len && cgp.len!=cgb->qdir.len) { e=2; } else
   if (!cgw) { if (rtype!=CGR_ABELIAN) e=3; } else
   if (!wscalar1()) { e=4; }

   if (e && F) wblog(F,L,
      "ERR %s() invalid scalar CRef (e=%d)\n%s",FCT,e,STR_(this));
   return e;
};

template <class TQ>
int CRef<TQ>::check(const char *F, int L) const {

   int e=0; 

   if (!cgb) {
      if (rtype==CGR_ABELIAN) {
         if (cgp || cgw) { e=1; }
      }
      else if (rtype==CGR_CTR_SCALAR) {
         if (cgp || !wscalar1()) { e=2; }
      }
      else if (cgp || cgw) { e=3; }

      if (e && F) { wblog(F_L,
         "ERR %s() invalid scalar CRef (e=%d)\n%s",FCT,e,STR_(this));
      }
   }
   else {
      if (!cgw) e=11; else
      if (wdim1()>cgb->getOM()) e=12;

      if (e && F) { wblog(F_L,
         "ERR %s() invalid CRef data (e=%d, m=%d)\n%s",
         FCT,e,cgb->getOM(), STR_(this));
      }
   }

   return e;
};

template <class TQ>
int CRef<TQ>::check(const char *F, int L,
   const QType &q, const QDir &qdir
 ) const {

   int e=check(F,L); if (e) { return e; }

   if (q.isAbelian()) {
      if (cgb)  {
         if (sameQDir(FL,qdir)) { return (e=1); }
      }
      if (cgp) { if (F) wblog(FL,
         "ERR %s() got cgp=[%s] for %s",FCT,STR(cgp), STR(q));
         return (e=2);
      }
      if (cgw) {
         if (!wscalar1()) { if (F) wblog(FL,
            "ERR %s() got cgw=[%s] for %s",FCT,STR(cgw), STR(q));
            return (e=3);
         }
      }
   }
   else {
      if (!cgb) {
         if (rtype==CGR_CTR_SCALAR && wscalar1()) {
            return e;
         }
         if (qdir.len) { 
            if (F) wblog(F,L,"ERR %s() got NULL cref for %s: %s (%s)",
               FCT,STR(q),STR(qdir),STR(cgw));
            return (e=12);
         }
         return (e=11);
      }
      if (sameQDir(FL,qdir)) { return (e=12); }
   }

   return e;
};

template <class TQ>
QSet<TQ>& CRef<TQ>::adapt(QSet<TQ> &Q, char iflag) const {
   if (cgp.len) Q.Permute(cgp,iflag ? 0 : 'i'); 
   if (conj) Q.qdir.Conj();
   return Q;
};

template <class TQ>
ctrIdx& CRef<TQ>::adapt(ctrIdx &I, char iflag) const {

   if (cgp.len && !cgp.isIdentityPerm()) {
      unsigned i=0; wbvector<unsigned> x(I.len,I.data);
      if (I.len>cgp.len) wblog(FL,"ERR %s() "
         "unexpected index length (%d/%d)",FCT,I.len,cgp.len);
      if (iflag==0) {
         for (; i<I.len; ++i) { I[i]= cgp.el(x[i]); }}
      else { wbperm cpi(cgp,'i');
         for (; i<I.len; ++i) { I[i]= cpi.el(x[i]); }
      }
   }
   if (conj) I.Conj();
   return I;
};

template <class TQ>
bool CRef<TQ>::wscalar(const char *F, int L) const {
   if (!cgw.data) wblog(F_L,"ERR %s() got empty cgw",FCT);
   if (F || L)
        { return (isw2(F_L)==1); }
   else { return cgw.numel()==1; } 
};

template <class TQ>
double CRef<TQ>::wel(unsigned k) const { 
   double rval=1; if (cgw || k) {
      if (cgw.data && k<cgw.numel()) { rval=cgw.data[k]; }
      else wblog(FL,
      "ERR %s() index out of bounds (%d / %s)",FCT,k,SSTR(cgw));
   }
   return rval;
};

template <class TQ>
unsigned CRef<TQ>::wdim(const char *F, int L) const {
   if (!isw2() || !cgw.isSMatrix()) {
      if (F || L) wblog(F_L,
         "ERR %s() invalid cgw (s=%s)",FCT,SSTR(cgw));
      return 0;
   }
   return cgw.SIZE[0];
};

template <class TQ>
unsigned CRef<TQ>::wdim12_(const char *F, int L, unsigned &d2) const {
   unsigned d1=0; 
   if (cgw) { d1=isw2(F_L); d2=cgw.SIZE[1]; }
   else if (!F)
        { d2=0; }
   else { wblog(F,L,"ERR %s() got empty wdim: %s",FCT,STR_(this)); }
   return d1;
};

template <class TQ>
unsigned CRef<TQ>::wdim12(const char *F, int L, unsigned &d2,
   unsigned d0 
 ) const {

   unsigned d1=d0; 
   if (cgw) { d1=isw2(F_L); d2=cgw.SIZE[1]; }
   else {
      if (!isAbelian(F,L)) wblog(FL,
         "ERR %s() non-abelian CRef got empty cgw\n%s",FCT,STR_(this));
      d2=d0;
   }
   return d1;
};

template <class TQ>
unsigned CRef<TQ>::wdim3(const char *F, int L) const {
   if (!cgw || cgw.SIZE.len!=3) { if (F || L) { wblog(F_L,
      "ERR %s() invalid cgw size %s",FCT,SSTR(cgw)); }
      return -1;
   }
   return cgw.SIZE[2];
};

template <class TQ>
double CRef<TQ>::wget0() const {
   double x=1; 
   if (cgw) { unsigned l=cgw.numel();
      if (l==1)
           { x=cgw.data[0]; }
      else { wblog(FL,"ERR %s() got cgw of size %s",FCT,SSTR(cgw)); }
   }
   return x;
};

template <class TQ>
double CRef<TQ>::wget1() const {

   unsigned r=cgw.isScalar();
   if (r<2 || r>3) {
      if (!cgw)
           wblog(FL,"ERR %s() got empty cgw",FCT);
      else wblog(FL,"ERR %s() got cgw of size %s",FCT,SSTR(cgw));
   }
   return cgw.data[0];
};

template <class TQ>
size_t CRef<TQ>::isw2(const char *F, int L) const {
   const size_t *s=cgw.SIZE.data;

   if (cgw.SIZE.len==2 && s[0]>=s[1] && s[1]) {
      return s[0]; 
   }
   if (F || L) wblog(F_L,
      "ERR %s() invalid size(cgw) = %s",FCT,SSTR(cgw));
   return 0;
};

template <class TQ>
int CRef<TQ>::cgw_check(const char *F, int L) const {
   int q=0; 
   if (cgw) {
      if (cgw.SIZE.len!=2) { q|=1; if (F) wblog(F,L,
         "ERR %s() invalid cgw (S=%s)",FCT,SSTR(cgw)); }
      if (!cgw.isFinite()) { q|=2; if (F) wblog(F,L,
         "ERR %s() invalid cgw (inf/nan values)",FCT); }
   }; return q;
};

template <class TQ>
unsigned CRef<TQ>::Reduce_w3Id(unsigned d_, char lflag) {

   unsigned d=1;
   if (cgw.SIZE.len==3) { d=cgw.SIZE[2];
      if (int(d_)>=0) {
         if (d<d_) { d=d_; } else 
         if (d>d_) wblog(FL,"ERR %s() OM inconsistency (m=%d/%d)",FCT,d,d_);
      }
      cgw.initIdentity(d);
   }
   else if (cgw.SIZE.len==2) {
      if (!wscalar()) wblog(FL,
         "ERR %s() invalid x3 %s for %s",FCT,SSTR(cgw),STR_(this));
      cgw.data[0]=1; 
   }
   else if (cgb || cgw) wblog(FL,
      "ERR %s() invalid x3 CRef %s: %s (%p)",FCT,STR_(this),SSTR(cgw),cgb);
   else if (!isAbelian()) {
      if (lflag || rtype==CGR_CTR_ZERO) { d=0; }
      else wblog(FL,"ERR %s() invalid x3 CRef %s",FCT,STR_(this));
   }

   return d;
};

template <class TQ>
int CRef<TQ>::wisId() const {
   int rval=0;
   if (cgw.SIZE.len==2) { double x=CG_SKIP_DEPS2;
      if (cgw.isIdentityMatrix(&x)) {
         rval=wdim(0,0); 
         if (!rval) { rval=-int(cgw.SIZE[0]); }
      }
   }
   return rval;
};

template <class TQ>
bool CRef<TQ>::cgw_exists(unsigned i, unsigned j) {
   bool q=0; 
   if (cgw) {
      if (cgw.SIZE.len!=2) wblog(FL,
         "ERR %s() invalid cgw size %s",FCT,SSTR(cgw));
      if (i<cgw.SIZE[0] && j<cgw.SIZE[1]) { q=1; }
   }; return q;
};

template <class TQ>
int CRef<TQ>::cgw_check_std3(
   const char *F, int L, unsigned d3, unsigned im) {

   int e=0;

   int fflag=(L || F); if (L<0) { L=0; }

   if (!cgw.isSMatrix()) { e|=1; if (fflag) wblog(F_L,
      "ERR %s() expecting square cgw (%s)",SSTR(cgw)); }
   if (!cgw.isDiagMatrix()) { e|=2; if (fflag) wblog(F_L,
      "ERR %s() expecting diagonal cgw data\n%s",FCT,STR(cgw)); }
   if (im>=cgw.SIZE[0]) { e|=4; if (fflag) wblog(F_L,
      "ERR %s() cgw index out of bounds (%d/%d)",im,cgw.SIZE[0]); }

   double x=cgw(im,im); x*=x;
   if (x<1 || fabs(x-round(x))>1E-12) { e|=8; if (fflag)  wblog(F_L,
      "ERR %s() got cgw(%d,%d)=%.3g",FCT,im,im,cgw(im,im)); }

   if (int(d3)<=0) {
      if (cgb && cgb->cgd.SIZE.len>=3) { d3=cgb->cgd.SIZE[2]; }
      else d3=0;
   }

   if (d3) {
      if (fabs(x-d3)>1E-12 && (x!=1 || im)) { e|=16; if (fflag)
         wblog(F_L,"ERR %s() got cgw[%d]=%.3g [expecting ±sqrt(%d)]",
         FCT,im,cgw(im,im),d3);
      }
   }

   return e;
};

template <class TQ>
bool CRef<TQ>::wSame(
   const wbarray<double> &cwB, char lenient, double eps) const {

   unsigned na=cgw.numel();
   unsigned nb=cwB.numel();

   if (na<=1 && nb<=1) {
      if (na && nb)
           { return fabs(cgw[0]-cwB[0])<=eps; }
      else { return (na==nb); }
   }

   if (cgw.sameAs(cwB,eps)) { return 1; }
   if (!lenient || cgw.rank()!=2 || cwB.rank()!=2) { return 0; }
   else { 
      unsigned i,j=0,
      d1a=cgw.SIZE[0], d1b=cwB.SIZE[0], d1=MIN(d1a,d1b),
      d2a=cgw.SIZE[1], d2b=cwB.SIZE[1], d2=MIN(d2a,d2b);
      const double *a=cgw.data, *b=cwB.data;

      for (   ; j<d2; ++j) {
      for (i=0; i<d1; ++i) { 
         if (fabs(a[i+j*d1a] - b[i+j*d1b])>eps) { return 0; }
      }}
   }
   return 1;
};

template <class TQ>
bool CRef<TQ>::isDiagCSC(RTD eps) const {

   if (isScalar()) return 1; 

   if (!isRefInit(FL)) {
      const cdata__ &c=cgb->cgd;
      if (c.numel()>1 && !c.isDiagMatrix(eps)) return 0;
   }
   else {
      unsigned m=getOM(FL);
      if (cgb->cgd.D.len!=m) wblog(FL,"ERR %s() "
         "got invalid CGD_REF data (%d/%d)",FCT,cgb->cgd.D.len,m);
      return (cgb->cgd.isSMatrix());
   }
   return 1;
};

template <class TQ>
bool CRef<TQ>::isIdentityCG(
   wbvector<double> *nrm, 
   double eps
 ) const {

   bool rval=0;

   if (isScalar()) { 
      if (nrm && !nrm->data) { nrm->init2val(1,1.); }
      return (rval=1);
   }

   if (!cgb) { wblog(FL,"ERR %s() got null cgb", FCT); }
   if (!cgw) { wblog(FL,"ERR %s() got empty cgw",FCT); }
   if (!cgb->isSymmmetric()) { return rval; }

   unsigned r=rank(FL);

   if (!r || r!=cgb->qdir.len) wblog(FL,
      "ERR %s() got r=%d/%d",FCT,r,cgb->qdir.len);
   if (r%2 || (r<3 && wnumel()!=1)) { return rval; }

   wbvector<RTD> cgt; 
   wbindex ip;
   cgb->qdir.findGT(0,ip);

   if (!ip.len || 2*ip.len!=cgb->qdir.len) { return rval; }

   if (cgb->cstat!=CGD_REF_INIT) {

      trace(FL,cgt);

      if (cgb->cstat==CGD_REF_INIT) wblog(FL,
         "ERR %s() got mixed ref_init\n%s",FCT,STR_(this));
      if (cgt.len && cgt.len!=wdim()) wblog(FL,
         "ERR %s() got OM mismatch (%d; %s)",FCT,cgt.len,SSTR(cgw));

      if (cgt.len==1) { RTD x=RTD(1)/cgt[0];
         rval=cgb->cgd.isProptoId(x,eps); 
      }
      else {
         rval=1; 
      }
   }
   else {
      if (!isRefInit() || !cgb->cgd.SIZE.len) wblog(FL,
         "ERR %s() got invalid CGD_REF data\n%s",FCT,STR_(this));
      if (nrm) {
         cgt.init2ref(cgb->cgd.D); 
         if (!cgt) { nrm=NULL; }
      }
      rval=1;
   }

   if (!nrm) { return rval; }

   unsigned i=0, j, m=wdim1(), n=wdim2(); 
   wbvector<double> x_(n); double *x=x_.data, ti;

   if (cgt.len!=m || !m) wblog(FL,
      "ERR %s() got OM mismatch (%d / %s)",FCT,cgt.len,SSTR(cgw));

   for (; i<m; ++i) { ti=double(cgt[i]);
      if (fabs(ti)<eps) wblog(FL,
          "ERR %s() got small Id fac (%d/%d: %g)",FCT,i+1,m,ti);
      for (j=0; j<n; ++j) { x[j]+=cgw(i,j)/ti; }
   }

   if (!nrm->data || !n) { x_.save2(*nrm); }
   else if (n==1) { (*nrm)*=x[0]; }
   else { nrm->tensorProd(x_,*nrm); } 

   return rval;
};

template <class TQ>
double CRef<TQ>::NormStd(const char *F, int L, unsigned r)  {

   double nrm=1; 
   if (isScalar()) { return nrm; } 
   if (!cgb) wblog(F_L,"ERR %s() got null cgb",FCT);

   wbindex ip,in;
   unsigned i;

   if (int(r)>=0 && r!=cgb->qdir.len) wblog(F_L,
      "ERR %s() got unexpected CRef of rank %d/%d",FCT,r,cgb->qdir.len);

   cgb->qdir.findGT(0,ip,in);

   if (in.len==1) { i=in[0]; } else
   if (ip.len==1) { i=ip[0]; }
   else {
      if (F) { wblog(F_L,
         "ERR %s() got unexpected CRef %s",FCT,STR_(this)); }
      return nrm;
   }

   i=cgb->cgd.dim0(i);
   if (i!=1) {
      if (int(i)<=0) wblog(F_L,"ERR %s() invalid dim=%d",FCT,i);
      nrm=sqrt(double(i)); cgw*=nrm; 
      nrm=1./nrm;
   }

   return nrm;
};

template <class TQ>
double CRef<TQ>::normExt(const char *F, int L) const {

   double rval=1; 
   unsigned r,d;

   if (isScalar()) { return rval; } 
   if ((r=rank(F_L))>3) { return rval; } 

   if (!cgb) wblog(FL,"ERR %s() got NULL cgb",FCT);
   if (r!=cgb->qdir.len) wblog(FL,
      "ERR %s() rank mismatch (%d/%d)",FCT,r,cgb->qdir.len);

   if (r==3) {
      int i=cgb->qdir.find1();
          if (i<-99) { return rval; } 
      d=cgb->cgd.dim0(abs(i));
   }
   else { d=cgb->cgd.maxDim(r); } 

   if (int(d)<=0) wblog(FL,"ERR %s() invalid qdim=%d",FCT,d);

   return (rval=sqrt(double(d)));
};

template <class TQ>
double CRef<TQ>::norm2(char checks) const { 

   double w2=0; 

   if (!cgb || isAbelian()) {
      unsigned n=cgw.numel();
      if (n>1 || (n==1 && fabs(cgw[0]-1)>CG_SKIP_DEPS1)) wblog(FL,
         "ERR %s() invalid CRef %s (n=%d)",FCT,STR_(this),n);
      return (w2=1);
   }

   if (!cgw) wblog(FL,"ERR %s() got empty cgw (%s)",FCT,SSTR(cgw));
   if (!cgw.isOrthoCols(&w2,'t',CG_SKIP_DEPS1)) wblog(FL,
      "ERR %s() got non-orthogonal cgw",FCT);

   if (!isfinite(w2)) { wblog(FL,"ERR %s() w2=%g",FCT,w2); }

   if ((checks&1) && fabs(w2-1)>CG_SKIP_DEPS1) {
      double x=normExt(FL); 
      if (w2<1 || fabs(w2-x*x)>CG_SKIP_DEPS1) { wblog(FL,
         "ERR %s() got CRef with |cgw|^2=%g (expected %g)",FCT,w2,x*x);
      }
   }

   if (cgb) { 
      if (checks>1 && !isRefInit()) {
         cgb->checkNormSign(FL, checks<<1);
      }
      else {
         unsigned r=cgb->rank(FL); check(FL);
         if (r<2) {
            if (r!=1 || cgb->qs.norm2()) 
            wblog(FL,"ERR %s() got r=%d",FCT,r);
         }
      }
   }

   return w2;
};

template <class TQ> inline
bool CRef<TQ>::isSortedDegQ(QSet<TQ> *QS, wbperm &pxt, char &cxt) const {

   if (!cgb) {
      if (pxt.len) { pxt.init(); }; cxt=0;
      return 1;
   }

   QSet<TQ> Q(*cgb);
      adapt(Q,'i'); if (QS) QS->init(Q); 
      Q.Sort(&pxt, &cxt,'i');

   if (cgp.sameAs(pxt) && conj==cxt) {
      return 1;
   }
   return 0;
};

template <class TQ>
int CRef<TQ>::SortDegQ(const char *F, int L, QSet<TQ> *QS) {

   wbperm pxt; char cxt=0;
   if (isSortedDegQ(QS,pxt,cxt)) { return 0; }

   if (cgb) {
      if (cgb->t.isAbelian()) {
         if (cgb->cgd.D.len && (cgb->cgd.D.len>1 || cgb->cgd.D[0]!=1))
            wblog(F_L,"ERR %s() invalid abelian CRef data\n%s",FCT,
            STR_(this)
         );
      }
      else if (cgb->qdir.len==2 && cgb->qdir.prod()<0) {
         if (isRefInit()) {
            const wbvector<unsigned> &S=cgb->cgd.SIZE;
            if (S.len!=2 || S[0]!=S[1] || cgb->cgd.D.len>1) wblog(F_L,
               "ERR %s() invalid ref CData (%d,%d)\n%s",
               FCT, cgb->cgd.SIZE.len, cgb->cgd.D.len, STR_(this)
            );
         }
         else {
            if ((cgb->cgd.SIZE.len && cgb->cgd.SIZE.len!=cgb->qdir.len) ||
               !cgb->cgd.D.len || !wscalar()) { wblog(F_L,
               "ERR %s() invalid scalar CData (%d,%d)\n%s",
            FCT, cgb->cgd.SIZE.len, cgb->cgd.D.len, STR_(this)); }
         }
      }
      else {
         CRef<TQ> R=*this;
         pxt.save2(R.cgp); R.conj=cxt;

         if (cgb->qdir.len<2) wblog(FL,"ERR %s()",FCT,STR_(this));

         gXS.contractDegQ(FL,*this,R); 

         R.save2(*this);
         return 2;
      }
   }

   pxt.save2(cgp); conj=cxt;
   return 1;
};

template <class TQ>
double CRef<TQ>::trace(const char *F, int L,
   const ctrIdx &i1, const ctrIdx &i2, CRef<TQ> *Rt) const {

   double q=0; wblog(FL,"ERR %s() to be cont'd",FCT);
   return q;
};

template <class TQ>
wbvector<RTD>& CRef<TQ>::trace(
   const char *F, int L, wbvector<RTD> &cgt) const {

   if (!cgb) { return cgt.init(); }
   if (!cgb->isSymmmetric()) {
      if (F) wblog(F_L,"ERR %s() got %s",FCT,STR_(this));
      return cgt.init();
   }
   else {
      if (isRefInit())
           { return cgt.init(cgb->cgd.D); }
      else { return cgb->trace(F_L, cgt); }
   }
};

template <class TQ>
double CRef<TQ>::trace() const {

   if (!cgb) {
      if (cgw) wblog(FL,
         "ERR %s() invalid cgw %s (cgb=0)",FCT,SSTR(cgw));
      return 1;
   }
   if (!isSymmmetric()) wblog(FL,
      "ERR %s() got %s\ncgp=%s",FCT,STR_(this),STR(cgp));

   unsigned i,j,m, n=cgb->gotOM(FL);
   wbvector<RTD> cgt;

   if (!n) n=1; else 
   if (n<2) wblog(FL,"ERR %s() got OM=%d",FCT,n);

   if (isRefInit())
        cgt.init2ref(cgb->cgd.D); 
   else cgb->trace(FL,cgt);

   if ((m=wdim1())>n || cgt.len!=n) wblog(FL,"ERR %s() "
      "got OM size mismatch (%s %d/%d)",FCT,SSTR(cgw),cgt.len,n);
   n=wdim2();
   wbvector<double> x1(n); double tj, *x=x1.data;

   for (j=0; j<m; ++j) { tj=double(cgt[j]);
   for (i=0; i<n; ++i) { x[i]+=(cgw(i,j)*tj); }}

   for (j=1; j<m; ++j) { if (fabs(x[j]-x[0])>1E-12) wblog(FL,
      "ERR %s() got varying trace %g @ %.3g",FCT,x[j],fabs(x[j]-x[0]));
   }

   return x[0];
};

template <class TQ> 
CRef<TQ>& CRef<TQ>::Permute(
   const wbperm &P, char iflag,
   char isnew 
) {

   if (!P.len) { return *this; }
   if (!cgb) { 
      if (cgw || cgp) wblog(FL,   
         "ERR %s() invalid CRef\n%s",FCT,STR_(this));
      return *this;
   }

   unsigned r=cgb->rank(FL);

   if (isnew && cgp.len) wblog(FL,
      "ERR %s() got existing permutation cgp=[%s]",FCT,STR(cgp));

   if (cgb->cgd.isScalar()) { if (!wscalar()) wblog(FL,
      "ERR %s() got cgw %s having %s",FCT,SSTR(cgw),cgb->QStr().data);
   }

   cgp.Permute(P,iflag,r);

   if (cgp.isIdentityPerm()) {
      cgp.init();
   }
   else {
      SortDegQ(FL);
   }

   return *this;
};

template <class TQ>
double CRef<TQ>::NormSignW(
   const char *F, int L,
   char useExt, 
   double eps, double eps2
){
   double wnrm=1; 
   if (!cgb || cgb->isAbelian()) { checkAbelian(F_L); return wnrm; }

   unsigned n=wdim1(), m=cgb->getOM(F_L);

   if (n==1 && m==1) { 
      if (fabs(cgw[0])<1-eps) wblog(F_L,"ERR %s() "
         "got cgw=%.3g without OM\nhaving %s",FCT,cgw[0],STR_(this));
      if (useExt) 
           { wnrm=cgw[0]/normExt(F_L); cgw[0]/=wnrm; }
      else { wnrm=cgw[0]; cgw[0]=1; }

      return wnrm;
   }
   if (n>m) wblog(F_L,
      "ERR %s() cgw out of bounds (%d/%d)",FCT,SSTR(cgw),m);

   unsigned i;
   double a, w2=0;

   if (!cgw.isOrthoCols(&w2,'t',CG_SKIP_DEPS1) || w2<1E-8) wblog(F_L,
      "ERR %s() got non-orthogonal cgw (w2=%.3g)",FCT,w2);

   wnrm=sqrt(w2);
   if (useExt) { wnrm/=normExt(F_L); } 

   if (wnrm<=eps) wblog(F_L,
      "ERR %s() got small wnrm %.3g (%s)",FCT,wnrm,STR_(this));

   for (n=wnumel(), i=0; i<n; ++i) { if (fabs(cgw[i])>eps) {
      if (cgw[i]<0) { wnrm=-wnrm; }
      break;
   }}

   for (i=0; i<n; ++i) { if ((a=fabs(cgw[i]))<eps) {
      if (a>eps2) wblog(F_L,
         "WRN %s() cgw noise %.3g (%g, %g; i=%d %s)",
         FCT,cgw[i],eps,eps2,i,SSTR(cgw));
      cgw[i]=0;
   }}

   if (wnrm!=1) { cgw /= wnrm; }

   return wnrm;
};

template <class TQ>
char CRef<TQ>::sameUptoFac(const CRef &B,
   double *fac_, 
   double eps) const {

   if ((cgb==NULL) ^ (B.cgb==NULL)) wblog(FL,
      "ERR CRef::%s() inconsistency in cgb",FCT);
   if (cgb!=B.cgb) { return 1; }

   if (!cgb) {
      if (cgw || B.cgw) wblog(FL,"ERR %s() "
         "invalid CRef data (%s / %s; cgb=0)",FCT,SSTR(cgw),SSTR(B.cgw));
      if (fac_) { (*fac_)=1; }
      return 0; 
   }

   if (!cgw && !B.cgw) {
      if (fac_) (*fac_)=1;
      return 0;
   }

   if (!cgb->cstat.sameAs(B.cgb->cstat,'L')) { return 2; }

   if (anyTrafo() || B.anyTrafo()) {
      if (!sameQSet(B)) return 3;
   }

   unsigned i=0, j=0, n1=wdim1(), n2=B.wdim1(), N1,N2;
   const size_t *sa=cgw.SIZE.data, *sb=B.cgw.SIZE.data;
   unsigned m=cgb->getOM();
   double a=cgw.aMax(i,j, &B.cgw);

   if (!n1 || !n2 || n1>m || n2>m) wblog(FL,
      "ERR %s() cgw out of bounds (%s/%s/%d)",FCT,SSTR(cgw),SSTR(B.cgw),m);
   if (a<=eps) wblog(FL,
      "ERR %s() got small cgw (%.3g)",FCT,cgw.norm());
   if (ABS(B.cgw(i,j))<=eps) { return 4; } 

   double fac=B.cgw(i,j)/cgw(i,j); if (fac_) { (*fac_)=fac; }

   n1=MIN(sa[0],sb[0]); N1=MAX(sa[0],sb[0]);
   n2=MIN(sa[1],sb[1]); N2=MAX(sa[1],sb[1]);

   for (i=j=0; j<N2; ++j, i=0) { 
      if (j<n2) {
         for (; i<n1; ++i) {
            a=ABS( fac*cgw(i,j) - B.cgw(i,j) );
            if (a>eps) { return 10; }
         }
      }
      if ((j<n2 && n1<sa[0]) || (j>=n2 && n2<sa[1]))
           { for (; i<N1; ++i) { if (ABS(  cgw(i,j)) > eps) return 11; }}
      else { for (; i<N1; ++i) { if (ABS(B.cgw(i,j)) > eps) return 12; }}
   }

   return 0; 
};

template <class TQ>
char CRef<TQ>::sameAs(const CRef &B, double eps) const {

   char q=0; double fac=0.;
   q=sameUptoFac(B,&fac,eps); 
   if (!q && fabs(fac-1)>eps) { q=6; }

   return q;
};

template <class TQ>
char CRef<TQ>::sameAs_fix(
   const char *F, int L, CRef<TQ> &B, double eps) {

   char q=0; double fac=0.;

   if (cgb!=B.cgb) {
	  if ((cgb && B.cgb) || !isAbelian() || !B.isAbelian()) {
         if (F || L) wblog(F_L,
            "ERR CGR data mismatch (%d,%d)\n%s <> %s",STR_(this), STR(B));
         q=7;
      }
      else {
         if (!cgb)
              { cgb=B.cgb; if (!  cgw) { cgw=B.cgw; }}
         else { B.cgb=cgb; if (!B.cgw) { B.cgw=cgw; }}
      }
   }

   if (!q) { q=sameUptoFac(B,&fac,eps); } 
   if (!q && fabs(fac-1)>eps) { q=6; }
   return q;
};

template <class TQ>
char CRef<TQ>::got3(const char *F, int L,
   const QType &q, const qset<TQ> &J12, const qset<TQ> &J3
 ) const {

   if (!cgb) wblog(F_L,"ERR %s() got empty cgb",FCT);

   const CDATA_TQ &C = *cgb;

   if (q!=C.t) {
      if (F) wblog(FL,"ERR %s() symmetry mismatch (%s <> %s)",
         FCT,STR(q),STR(C.t));
      return 1;
   }
   if (C.qs.len!=J12.len + J3.len || J12.len!=2*J3.len) wblog(FL,
      "ERR %s() qset inconsistency (%d == 2*%d; %d)",FCT,
      J12.len,J3.len,C.qs.len);

   if (cgp.len) {
      if (cgp.len!=C.qdir.len) wblog(FL,"ERR %s() "
         "invalid cgp=[%s] %d",FCT,STR(cgp),C.qdir.len);

      unsigned i=0, l=J12.len;
      qset<TQ> qx; C.qs.blockPermute(cgp,qx);
      const TQ *qs=qx.data;

      for (; i<J12.len; ++i) {
         if (qs[i]!=J12.data[i]) {
            if (F) wblog(FL,
               "ERR %s() qset mismatch: %s @ %s <> (%s | %s); i=%d",
               FCT, STR(qx), STR(cgp),
               J12.toStr().data, STR(J3), i+1);
            return 2;
         }
      }
      for (; i<J3.len; ++i) {
         if (qs[i+l]!=J3.data[i]) {
            if (F) wblog(FL,
               "ERR %s() qset mismatch (%s @ %s <> %s ; %s; %d)",
               FCT, STR(qx), STR(cgp),
               J12.toStr().data, STR(J3), i+l+1);
            return 3;
         }
      }
   }
   else {
      const TQ *qs=C.qs.data;
      if (Wb::cmpRange(qs, J12.data, J12.len)) { return 2; }; qs+=J12.len;
      if (Wb::cmpRange(qs, J3 .data, J3 .len)) { return 3; }
   }

   return 0;
};

template <class TQ>
bool CRef<TQ>::sameQDir(const char *F, int L, const QDir &qd) const {

   if (!cgb) return 1;

   if (cgb->qdir.len!=qd.len) wblog(F_L,"ERR %s() QDir length mismatch "
      "(%s <> %s)",FCT, STR(cgb->qdir), STR(qd));
   if (cgp.len && cgp.len!=qd.len) wblog(F_L,
      "ERR %s() got length mismatch (%d/%d/%d)",
      FCT,cgb->qdir.len, qd.len, cgp.len
   );

   char e=0, cflag=gotConj(), pflag=gotPerm(); 
   if (!cflag && !pflag) return qd==cgb->qdir;

   unsigned i=0;
   const char *d=cgb->qdir.data, *d0=qd.data;

   if (pflag) {
      const wperm_t *p=cgp.data;

      if (cflag) {
         for (; i<qd.len; ++i) { if (d[p[i]]!=-d0[i]) { e=1; break; }}
      }
      else {
         for (; i<qd.len; ++i) { if (d[p[i]]!= d0[i]) { e=2; break; }}
      }
   }
   else { 
      for (; i<qd.len; ++i) { if (d[i]!=-d0[i]) { e=3; break; }}
   }

   if (e) { char s[128]; snprintf(s,128,
      "%s <> %s (e=%d)", STR_(this), STR(qd), e);
      if (F) wblog(F,L,"ERR %s() got\n%s",FCT,s);
      else   wblog(FL, "WRN %s() got\n%s",FCT,s);
   }

   return (e? 0:1);
};

template <class TQ>
bool CRef<TQ>::sameQSet(const CRef &B) const {

   bool ipa=gotPerm(FL), ipb=B.gotPerm(FL);

   if (cgb!=B.cgb) { return 0; } 
   if (!cgb) { return 1; }       

   if (!ipa && !ipb) { return (conj==B.conj); } 

   const QSet<TQ> &Q=(QSet<TQ>&)(*cgb);
   unsigned i=0, r=Q.qdir.len;

   if ((cgp.len && cgp.len!=r) || (B.cgp.len && B.cgp.len!=r)) {
      wblog(FL,"ERR %s() got cgp.len=%d,%d/%d",
      FCT,cgp.len,B.cgp.len,r);
   }

   if (Q.qs.blockCompare(cgp, Q.qs, B.cgp)) { return 0; }

   const wperm_t *pa=cgp.data, *pb=B.cgp.data;
   const char *d=Q.qdir.data;

   if (ipa) {
      if (ipb) {
         if (conj ^ B.conj)
              for (; i<r; ++i) { if (d[pa[i]]!=-d[pb[i]]) return 0; }
         else for (; i<r; ++i) { if (d[pa[i]]!=+d[pb[i]]) return 0; }
      }
      else {
         if (conj ^ B.conj)
              for (; i<r; ++i) { if (d[pa[i]]!=-d[i]) return 0; }
         else for (; i<r; ++i) { if (d[pa[i]]!=+d[i]) return 0; }
      }
   }
   else {
      if (ipb) {
         if (conj ^ B.conj)
              for (; i<r; ++i) { if (d[i]!=-d[pb[i]]) return 0; }
         else for (; i<r; ++i) { if (d[i]!=+d[pb[i]]) return 0; }
      }
   }

   return 1;
};

template <class TQ>
int CRef<TQ>::gotSameCData(
   const char *F, int L, const CRef& B, char xflag) const {

   int rval=0; 

   if (cgb!=B.cgb) {
      if (F) wblog(F_L,"ERR %s() got CGC mismatch\n   %s\n<> %s",FCT,
            cgb ? STR_(  cgb) : "(null)",
          B.cgb ? STR_(B.cgb) : "(null)"
      );
      return (rval=-1);
   }

   if (!cgb) {
      if (cgw || B.cgw || rtype!=CGR_ABELIAN || rtype!=B.rtype)
         wblog(FL,"ERR %s() invalid scalar cref (%s/%s; %d/%d)",
         FCT, SSTR(cgw),SSTR(B.cgw), rtype.t, B.rtype.t);
      return (rval=0);
   }

   unsigned r=cgb->rank(FL), m=cgb->getOM(), n1=wdim1(), n2=B.wdim1();

   if (r<2) { if (F) wblog(FL,
      "ERR %s() got rank-%d CData (%s)",FCT,r,cgb->sizeStr().data);
      return (rval=-2);
   }

   if (n1>m || n2>m) {
      if (F) { wblog(FL,
         "ERR %s() CData out of multiplicity range (%s,%s /%d)",
         FCT, SSTR(cgw), SSTR(B.cgw), m);
      }
      return (rval=-3);
   }

   if (!cgp.sameAs(B.cgp) || conj!=B.conj) {
      if (F) wblog(FL,
         "ERR %s() got cgp+conj mismatch (%s%s <> %s%s)",
         FCT, STR(cgp), conj?"*":"", STR(B.cgp), B.conj?"*":"");
      return (rval=-4);
   }

   if (xflag) {
      char ira=isRefInit(), irb=B.isRefInit();
      if (ira || irb) {
         if (ira ^ irb) wblog(FL,
            "ERR %s() got ref mismatch (%d/%d)",FCT,ira,irb);
      }
      else { cgb->checkNormSign(F_L,xflag); }
   }

   return (rval=MAX(n1,n2));
};

template <class TQ>
bool CRef<TQ>::sameQSet(const QSet<TQ> &Q) const {

   if (!cgb || cgb->isEmpty()) {
      if (cgp.len || conj) wblog(FL,
         "WRN %s() invalid scalar (0x%lX)\ncgp.len=%d, conj=%d : %s",
         FCT, cgb, cgp.len, conj, STR(Q));
      return Q.isEmpty();
   }

   if (anyTrafo()) {
      QSet<TQ> X((QSet<TQ>&)(*cgb));
         if (cgp.len) X.Permute(cgp);
         if (conj) X.Conj();
      return X==Q;
   }

   return ((QSet<TQ>&)(*cgb))==Q;
};

template <class TQ>
bool CRef<TQ>::sameQDir(const iTags& b) const {

   if (!cgb) {
      if (cgp.len || conj) {
         wblog(FL,"ERR %s() %s",FCT,STR_(this)); }
      return 1;
   }
   if (cgb->qdir.len!=b.len) { return 0; }
   else {
      unsigned i=0; char a,c;
      if (cgp.len && cgp.len!=b.len) wblog(FL,"ERR %s() "
         "invalid permutation (len=%d/%d)",FCT,cgp.len,b.len);
      for (; i<b.len; ++i) {
         a=cgb->qdir.el(cgp.el1(i));  
         c=b.data[i].isConj() ^ conj; 
         if ((c && a>=0) || ((!c) && a<=0)) { return 0; }
      }
   }

   return 1;
};

template <class TQ>
bool CRef<TQ>::affectsQDir() const {

   if (!cgb || !cgb->qdir.len || !cgp.len) return 0;
   if (cgp.len && cgp.len!=cgb->qdir.len) wblog(FL,
      "ERR %s() got length mismatch (%d/%d)",FCT,cgp.len,cgb->qdir.len
   );

   char cflag=gotConj(), pflag=gotPerm(); 
   if (!cflag && !pflag) return 0;

   if (pflag) {
      const char *d=cgb->qdir.data;
      const wperm_t *p=cgp.data; unsigned i=0;

      if (cflag) {
         for (; i<cgp.len; ++i) { if (d[i]!=-d[p[i]]) return 1; }
      }
      else {
         for (; i<cgp.len; ++i) { if (d[i]!= d[p[i]]) return 1; }
      }
   }
   else { 
      return 1;
   }

   return 0;
};

template <class TQ>
double CRef<TQ>::normDiff2(const char *F, int L, const CRef &B) const {

   double x2=0; 

   unsigned n=gotSameCData(F_L,B,'!'); 
   if (int(n)<=0) {
      if (!n)
           { return x2; } 
      else { wblog(FL,"ERR %s() got e=%d",FCT,n); }
   }

   if (n==1) { 
      x2 = wget0() - B.wget0();
      x2=NORM2(x2);
   }
   else {
      x2=cgw.normDiff2(B.cgw,'!'); 
   }

   return x2;
};

template <class TQ>
wbarray<double>& CRef<TQ>::wProd(const CRef &B, wbarray<double> &x2) const {

   unsigned n=gotSameCData(FL,B,'!');  
   if (int(n)<=0) {
      if (!n) { double x=1;
         x2.init(1,1,&x); 
      }
      else { wblog(FL,"ERR %s() got e=%d",FCT,n); }
   }
   else if (n==1) {
      double x=wget0()*B.wget0();
      x2.init(1,1,&x); 
   }
   else {
      unsigned ma=wdim1(), mb=B.wdim1();
      if (ma==mb) {
         Wb::MatProd(cgw,B.cgw,x2,'T');
      }
      else {
         wbvector<size_t> S; wbarray<double> wx;
         if (ma>mb) {
            S=cgw.SIZE; S[0]=mb;
            cgw.resize(S,wx); Wb::MatProd(wx,B.cgw,x2,'T');
         }
         else {
            S=B.cgw.SIZE; S[0]=ma;
            B.cgw.resize(S,wx); Wb::MatProd(cgw,wx,x2,'T');
         }
      }
   }
   return x2;
};

template <class TQ>
double CRef<TQ>::safeCpy(
   const char *F, int L, const CRef<TQ> &B, CRef<TQ> &C) const {

   double x=1; unsigned e=0;

   if (this->isEmpty() && B.isEmpty()) {
      if (!C.isEmpty()) wblog(F_L,"ERR CRef inconsistency");
      return x;
   }

   if (!(e=sameUptoFac(B,&x))) { C=(*this); }
   else {
      MXPut(FL).addP(this->toMX(),"a")
      .addP(B.toMX(),"b").add(C,"c").add(x,"x").add(e,"e");
      wblog(F_L,"ERR %s() CRef not the same (e=%d)",FCT,e);
   }

   return x;
};

template <class TQ>
wbstring CRef<TQ>::toStr(char lflag) const {

   wbstring sout; 
   unsigned l=0, n; char *s;

   if (!cgb) { sout.init(32); s=sout.data; n=sout.len;

      if (!cgw) {
         if (rtype==CGR_ABELIAN)
              { l=snprintf(s,n,"(abelian @ [1])"); }
         else { l=snprintf(s,n,"(%s @ [])", rtype.tostr()); }
      }
      else {
         if (wnumel()>1 && rtype!=CGR_CTR_SCALAR) { wblog(FL,"ERR %s() "
            "got cgw %s for %s (cgb=0)",FCT,SSTR(cgw),rtype.tostr());
         }
         l=snprintf(s,n,"%s to %s", 
         rtype==CGR_CTR_SCALAR ? "fully contracted": rtype.tostr(), STR(cgw));
      }

      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   }
   else {
      unsigned m=cgb->getOM();
      int is2=isw2(); 

      if (rtype>=CGR_NUM_TYPES) wblog(FL,"ERR %s() "
         "cgrType out of range (%d/%d)",FCT,rtype.t,CGR_NUM_TYPES);

      if (is2) { sout.init(118+10*wdim2()); } else  
      if (lflag) {sout.init(128); } 
      else { isw2(FL); } 
      s=sout.data; n=sout.len;

      QSet<TQ> Q(*cgb); adapt(Q,'i');
      l=snprintf(s,n,"%s", Q.toStr(statStr().data).data);

      if (l<n) {
         if (wscalar() && m==1) {
            l+=snprintf(s+l,n-l,", w=%g",cgw[0]); if (!is2 && l<n) {
            l+=snprintf(s+l,n-l," %s",SSTR(cgw)); } 
         }
         else if (is2) {
            bool isd=cgw.isDiagMatrix();
            char sm[8]=""; if (wdim1()!=m) snprintf(sm,8,"%d/",m);
            wbvector<double> cc = (isd ? cgw.getDiag() : cgw.normCols());

            l+=snprintf(s+l,n-l,", %s=[%s] %s%s", 
               isd? "w=diag":"|w|",cc.toStrf("%.4g"," ").data,sm,SSTR(cgw)
            );
         }
         else {
            l+=snprintf(s+l,n-l,", |w|=%.4g %s",cgw.norm(),SSTR(cgw));
         }
      }
      if ((cgp.len || conj) && l<n) { wbperm P(cgp,'i');
         l+=snprintf(s+l,n-l,", p=%s%s",
         (P+1).toStrf("%d","").data, conj?"*":"");
      }
      if (l>=n) wblog(FL,"ERR %s() "
      "string out of bounds (%d/%d)%N%N%s%N%N",FCT,l,n,s);
   }
   return sout;
};

template <class TQ>
wbstring CRef<TQ>::sizeStr() const {

   if (cgb) { return cgb->sizeStr(); }

   if (cgw) wblog(FL,
      "ERR %s() inconsistent cgw %s",FCT,SSTR(cgw));
   return (rtype==CGR_ABELIAN ? "(scalar)" : "[]");
};

template <class TQ>
wbstring CRef<TQ>::statStr(char vflag) const {
   wbstring s(32); 
   unsigned l=0, n=s.len;

   if (vflag) {
      l=snprintf(s.data,n,"%s",rtype.tostr());
      if ((cgp.len || conj) && l<n) l+=snprintf(s.data+l,n-l,
         " %s%s",STR(cgp), conj?"*":"");
      if (cgb && l<n) l+=snprintf(s.data+l,n-l," 0x%lX",(long)cgb);
   }

   if (l<n) {
   if (cgb && cgb->cstat.t) {
      if (rtype.t) {
         l+=snprintf(s.data+l,n-l,"%s;%s",
         rtype.tostr(), STR(cgb->cstat));
      } else {
         l+=snprintf(s.data+l,n-l,"%s",STR(cgb->cstat));
      }
   }
   else {
      if (rtype.t) {
         l+=snprintf(s.data+l,n-l,"%s",rtype.tostr());
      }
   }}

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};

template <class TQ, class TD>
double genRG_base<TQ,TD>::normDiff(const char *F, int L,
  const genRG_base<TQ,TD> &B, TD eps
) const {

   double e, emax=0;

   if (Sp.len!=B.Sp.len || Sz.len!=B.Sz.len || !Z.sameSize(B.Z))
       wblog(F_L,"ERR size mismatch (%d/%d; %d/%d; %s <> %s)",
       Sp.len, B.Sp.len, Sz.len, B.Sz.len,
       Z.sizeStr().data, B.Z.sizeStr().data);

   e=Z.normDiff(B.Z); if (e>double(eps)) {
      if (F) {
          MXPut(F_L).add(Z,"Z").add(B.Z,"Z2").add(e,"e");
          wblog(F_L,"ERR %s() %.3g",FCT,e);
      } else wblog(FL,"WRN %s() %.3g",FCT,e);
      return e;
   }
   emax=MAX(emax,e);

   e=J.normDiff(B.J); if (e>double(eps)) {
      if (F) {
         MXPut(F_L).add(Z,"Z").add(B.Z,"Z2")
         .add(J,"J").add(B.J,"J2").add(e,"e");
         wblog(F_L,"ERR %s() J (e=%.3g)",FCT,e);
      } else wblog(FL,"WRN %s() J (e=%.3g)",FCT,e);
      return e;
   }
   emax=MAX(emax,e);

   for (unsigned i=0; i<Sp.len; i++) {
      e=Sp[i].normDiff2(B.Sp[i])/Sp[i].numel();
      e=std::sqrt(ABS(e)); if (e>double(eps)) {
         wblog(F_L,"TST %s() Sp[%d] (e=%.3g)",FCT,i+1,e);
         return e;
      }
      emax=MAX(e,emax);
   }
   if (emax && emax>1E-12) wblog(FL,"TST %s() e=%g",FCT,emax);

   for (unsigned i=0; i<Sz.len; i++) {
      e=Sz[i].normDiff2(B.Sz[i]);
      e=std::sqrt(ABS(e)); if (e>double(eps)) {
         wblog(F_L,"TST %s() Sz[%d] (e=%.3g)",FCT,i+1,e);
         return e;
      }
      emax=MAX(e,emax);
   }
   if (emax && emax>1E-12) wblog(FL,"TST %s() e=%g",FCT,emax);

   return emax;
};

template <class TQ, class TD>
void genRG_base<TQ,TD>::checkOrthoCommRel(
   const char *F, int L, const genRG_base<TQ,TD> &B
 ) const {

   if (&B==this) {
      wblog(FL,"WRN %s() got same genRG_base set",FCT);
      return;
   }

   unsigned i=Sp.len+Sz.len, j=B.Sp.len+B.Sz.len, e=0;
   wbvector< const wbSparrayTD* > S1(i), S2(j);
   wbSparrayTD C; double x; char s[64];

   for (j=i=0; i<  Sp.len; ++i, ++j) S1[j]=&(  Sp[i]);
   for (  i=0; i<  Sz.len; ++i, ++j) S1[j]=&(  Sz[i]);

   for (j=i=0; i<B.Sp.len; ++i, ++j) S2[j]=&(B.Sp[i]);
   for (  i=0; i<B.Sz.len; ++i, ++j) S2[j]=&(B.Sz[i]);

   for (i=0; i<S1.len; ++i) {
   for (j=0; j<S2.len; ++j) {
      S1[i]->comm(FL,*S2[j],C);
      if ((x=double(C.norm()))>1E-10) {
         snprintf(s,64,"[ S(%d), S(%d) ] @ %.3g",i+1,j+1,x);
         e=1; break;
      }
      S1[i]->comm(FL,*S2[j],C,'N','C');
      if ((x=double(C.norm()))>1E-10) {
         snprintf(s,64,"[ S(%d), S(%d)' ] @ %.3g",i+1,j+1,x);
         e=2; break;
      }
   }  if (e) break; }

   if (e) {
      MXPut(FL,"a").addP(s,"info").add(*this,"R1").add(B,"R2")
        .add(*S1[i],"a").add(*S2[j],"b").add(C,"c").add(x,"x");
      wblog(F_L,"ERR %s() got non-commuting symmetries\n%s",FCT,s);
   }
};

template <class TQ, class TD>
genRG_base<TQ,TD>& genRG_base<TQ,TD>::Sort() {

    if (Z.isEmpty()) {
       MXPut(FL,"ans").add(*this,"R");
       wblog(FL,"ERR %s() got empty Z",FCT);
    }

    unsigned i; wbperm P; wbMatrix<double> z2(Z);
    z2.FlipCols(); 
    z2.sortRecs_float(P,-1); 

    Z.Set2Recs(P);

    for (i=0; i<Sp.len; ++i) { Sp[i].MatPermute(P); }
    for (i=0; i<Sz.len; ++i) { Sz[i].MatPermute(P); }

#ifdef CG_CHECK_MW_PERM
    if (P0.isEmpty()) P0=P;
    else P0.Select(P);
#endif

    return *this;
};

template <class TQ, class TD>
genRG_base<TQ,TD>& genRG_base<TQ,TD>::ApplyQFac(
   const wbvector<double> &qfac,
   const wbarray<double> *JM
){
   if (J.len!=Z.dim2) wblog(FL,
      "ERR %s() severe Q-label mismatch (%d/%d)",FCT,J.len,Z.dim2);
   if (qfac.len>1 && qfac.len!=Z.dim2) wblog(FL,
      "ERR %s() qfac: length mismatch (%d/%d)",FCT,qfac.len,Z.dim2);
   if (JM && (!JM->isSMatrix() || JM->SIZE[0]!=Z.dim2)) wblog(FL,
      "ERR %s() J-map size mismatch (%dx%d <> %s)",
       FCT,Z.dim1,Z.dim2, JM->sizeStr().data);

   unsigned i,j;

   if (qfac.len) {
      if (qfac.len==1) { double x=qfac[0]; J*=x; Z*=x; }
      else {
         const double *x=qfac.data;
         for (j=0; j<Z.dim2; ++j) { J[j]*=x[j];
         for (i=0; i<Z.dim1; ++i) { Z(i,j)*=x[j]; }}
      }
   }

   if (JM) {
      wbarray<double> x2, x(1,J.len); x.init_T(J.data);
      Wb::MatProd(x,*JM,x2);

      for (i=0; i<J.len; ++i) {
         J[i]=num2int(FL,x2[i]); 
      }
   }

   SkipTiny(); return *this;
};

template <class TQ, class TD>
void genRG_base<TQ,TD>::compareStdSU2(const char *F, int L) const {

   cgsparray xSp,xSz,xS2,xE;
   wbvector<double> xsz; 
   double e1,e2,e3;

   if (!q.isSU2() || J.len!=1) wblog(F_L,
      "ERR invalid %s qset record [%s]",STR(q), STR(J));
   if (Sp.len!=1 || Sz.len!=1) wblog(F_L,
      "ERR got invalid SU(2) generator set (Sp[%d], Sz[%d])",
      Sp.len,Sz.len);

   CG::get_SU2mat(J[0],xsz,xSp,xSz,xS2,xE);
   xSz*=2; xsz*=2; 

   e1=Sp[0].normDiff(xSp);
   e2=Sz[0].normDiff(xSz);
   e3=Z.normDiff(wbMatrix<double>().init2ref(xsz,'t'));

   if (e1>1E-12 || e2>1E-12 || e3>1E-12) {
      char s[]="got SU(2) inconsistency";
      MXPut(FL,"q").add(Sz[0],"Sz").add(xSz,"Sz_")
         .add(Sp[0],"Sp").add(xSp,"Sp_").add(Z,"Z").add(xsz,"xsz");
      if (e1>1E-12) wblog(F_L,"ERR %s (Sp @ %.3g)",s,e1);
      if (e2>1E-12) wblog(F_L,"ERR %s (Sz @ %.3g)",s,e2);
      if (e3>1E-12) wblog(F_L,"ERR %s ( Z @ %.3g)",s,e3);
   }
};

template <class TQ, class TD>
void genRG_struct<TQ,TD>::initCommRel(
   const char *F, int L,
   const genRG_base<TQ,TD> &R, char vflag
){
   unsigned i,j,k,p,n, np=R.Sp.len, nz=R.Sz.len;
   wbsparray<TD> C; double x;
   TD dz;

   if (!np || !nz) wblog(F_L,"ERR %s() got empty RSet (%d,%d)",FCT,np,nz);

   CR.init(np);

   for (i=0; i<np; ++i) {
      if (R.Sp[i].isEmpty()) wblog(FL,
         "ERR %s() got empty R.Sp[%d/%d]",FCT,i+1,np);

      R.Sp[i].comm(FL,R.Sp[i],C,'N','C');
      if (C.norm2()<(TD)(1E-10)) {
         MXPut(FL,"ans").add(R,"R").add(i+1,"i");
         wblog(F_L,"ERR %s() [Sp,Sp'] has norm 0",FCT);
      }

      CR[i].init(i,i,nz);

      wbvector<unsigned> &kk=CR[i].k;
      wbvector<double> &fac=CR[i].fac;

      for (p=k=0; k<nz; ++k) {
         x=double(C.froNorm2(R.Sz[k])/R.Sz[k].froNorm2(R.Sz[k]));
         if (ABS(x)>1E-10) { kk[p]=k; fac[p++]=x; }
      }
      if (!p) {
         MXPut(FL,"i").add(R,"R").add(CR,"CR").add(C,"C")
          .add(p+1,"p").add(np,"np");
         wblog(F_L,"ERR %s() got no overlap of [Sp,Sp'] with Sz "
         "(%d/%d)",FCT,i+1,np);
      }
      kk.len=p; fac.len=p;
   }

   DZ.init(np,nz);

   for (i=0; i<np; ++i) {
   for (j=0; j<nz; ++j) {
      R.Sz[j].comm(FL,R.Sp[i],C);
      if (C.sameUptoFac(R.Sp[i],&dz)) { DZ(i,j)=double(dz);
         MXPut(FL,"cr").add(R,"R").add(DZ,"DZ").add(i+1,"i").add(j+1,"j");
         wblog(FL,"ERR %s() invalid Sp/Sz operators (%s):\n"
         "[Sz,Sp] ~ Sp not satisfied",FCT,STR(q)); }
      DZ(i,j)=double(dz);
   }}

   checkCommRel(F_L,R); 

   if (vflag) {
      wblog(FL,"<i> summary of %s commutator relations",STR(q));
      for (i=0; i<CR.len; i++) { n=CR[i].k.len;
      for (k=0; k<n; k++) {
         printf("%4d: %2d %2d %2d : %8g\n", i+1,
         CR[i].i+1, CR[i].j+1, CR[i].k[k]+1, CR[i].fac[k]);
      }}
   }

};

template <class TQ, class TD>
double genRG_struct<TQ,TD>::checkCommRel(
  const char *F, int L, const genRG_base<TQ,TD> &R
) const {

   unsigned i,j,p,l,m, n=CR.len, np=R.Sp.len, nz=R.Sz.len;
   wbsparray<TD> C; double e2, r2=0;
   QType qt=(q==QTYPE_UNKNOWN ? R.q : q);

   if (q!=R.q) { char s[64];
      snprintf(s,64,"got QType mismatch: %s <> %s",STR(q),STR(R.q));
      if (q!=QTYPE_UNKNOWN)
           wblog(F_L,"ERR %s() %s",FCT,s);
   }

   if (qt!=QTYPE_UNKNOWN && (np!=nz || nz!=qt.sub)) wblog(F_L, 
      "ERR %s() unexpected Sz[%d], Sp[%d], sub=%d\nhaving symmetry='%s'",
      FCT, nz, np, qt.sub, STR(qt));
   if (!R.Sp.len || !R.Sz.len) wblog(F_L,
      "ERR %s() CR not yet initialized",FCT);

   for (i=0; i<nz; ++i) {
   for (j=i+1; j<nz; ++j) {
      e2=double(R.Sz[i].froNorm2(R.Sz[j]));
      if (e2>1E-10) { if (F) {
          MXPut(FL,"I_").add(R,"R").add(i+1,"i").add(j+1,"j").add(C,"C");
          wblog(F_L,"ERR %s() CR inconsistency:\n"
         "z-ops not mutually orthogonal (@ e=%.3g)",FCT,e2);
      }}
      R.Sz[i].comm(FL,R.Sz[j],C); e2=double(C.norm2()); r2+=e2;
      if (e2>1E-10) { if (F) {
         MXPut(FL,"I_").add(R,"R").add(i+1,"i").add(j+1,"j").add(C,"C");
         wblog(F_L,"ERR %s() CR inconsistency ([Z,Z] @ e=%.3g)",FCT,e2);
      }}
   }}

   for (l=0; l<n; ++l) {
      const wbvector<double> &fac=CR[l].fac;
      const unsigned *k=CR[l].k.data; i=CR[l].i; j=CR[l].j;

      if (i>=np || j>=np || CR[l].k.max()>=nz) wblog(F_L,
         "ERR %s() index out of bounds (%d,%d,%s/%d)",
          FCT, i+1, j+1, STR(CR[l].k+1), n);

      R.Sp[i].comm(FL,R.Sp[j],C,'N','C'); m=fac.len;

      for (p=0; p<m; ++p) { C.Plus(FL,R.Sz[k[p]],-fac[p]); }

      e2=double(C.norm2()); r2+=e2; if (e2>1E-10) { if (F) {
          MXPut(FL,"q").add(R,"R").add(i+1,"i").add(j+1,"j")
             .add(CR[l].k+1,"k").add(fac,"fac").add(C,"C");
          wblog(F_L,"ERR %s() CR inconsistency (%d,%d: e=%.3g)",
          FCT,i+1,j+1,e2);
      }}
   }

   if (DZ.dim1!=np || DZ.dim2!=nz) wblog(FL,"ERR %s() "
      "got empty DZ (%s; %s)",FCT,STR(q),DZ.sizeStr().data);

   for (i=0; i<np; ++i) {
   for (j=0; j<nz; ++j) {
      R.Sz[j].comm(FL,R.Sp[i],C); C.Plus(FL,R.Sp[i],-DZ(i,j));
      e2=double(C.norm2()); r2+=e2;

      if (e2>1E-10) { if (F) {
         MXPut(FL,"q").add(R,"R") 
         .add(R.Sz[j],"Sz").add(R.Sp[i],"Sp").add(DZ(i,j),"sfac")
         .add(*this,"I").add(i+1,"i").add(j+1,"j").add(C,"C");
         wblog(F_L,"ERR %s() CR inconsistency ([Z,Sp]: %.3g)",FCT,e2);
      }}
   }}

   r2/=MAX(widx_t(1), (n+R.Sz.len*(1+R.Sp.len))*R.Z.dim1); 

   e2=sqrt(r2);
   if (e2>cg_eps1) wblog(FL,"ERR %s() CR error = %.3g",FCT,e2); else
   if (e2>cg_eps2) wblog(FL,"WRN %s() CR @ %.3g",FCT,e2);

   return r2;
};

template <class TQ, class TD>
unsigned genRG_struct<TQ,TD>::genTensorProds(
    const qset<TQ> &J, unsigned flag) {

    unsigned i=0, l=0, m=0, n=0;
    wbvector< const qset<TQ>* > qq(RSet.size());

    for (auto I=RSet.begin(); I!=RSet.end(); ++I, ++i) {
       const genRG_base<TQ,TD> &R=I->second;
       if (!R.Sz.isEmpty() && R.Z.dim1 && R.Z.dim1<=32*q.qlen()) {
       qq[m++]=&(I->first); }
    }

    for (i=0; i<m; ++i) {
       l=getTensorProdReps_gen(J,*qq[i], flag | TP3_LDR); 
       if (int(l)>0) n+=l;
    }

    return n;
};

template <class TQ, class TD>
unsigned genRG_struct<TQ,TD>::genTensorProds(
    unsigned dmax, char sdig, char vflag) {

    unsigned i,j, l=0, n3=0; int e=0;
    if (int(dmax)<0) { dmax=MAX(8U,2*q.sub); } 
    unsigned dmax2=dmax*dmax;

    wbvector< const qset<TQ>* > qq(RSet.size());
    wbvector<unsigned> dd(RSet.size());
    Wb::SigHandler SIG(FL); 

    for (auto I=RSet.begin(); I!=RSet.end(); ++I) {
       const genRG_base<TQ,TD> &R=I->second;
       if (!R.Sz.isEmpty() && R.Z.dim1 && R.Z.dim1<=dmax2) {
          if (sdig && I->first.anyGT(9)) continue;
          qq[l]=&(I->first);
          dd[l]=R.Z.dim1; ++l;
       }
    }

    if (!l) wblog(FL,"ERR %s() no RSets found (d<=%d)",FCT,dmax);

    if (vflag) printStatus(FL,"old");

    for (i=0; i<l; ++i)
    for (j=0; j<l; ++j) { if (dd[i]*dd[j]<=dmax2) {
       if ((e=SIG.check911('t'))) {
          wblog(FL,"TST %s(%g,%s,%s) interrupt %d at (%d,%d)/%d",
             FCT,double(dmax),cSTR(sdig),cSTR(vflag),e,i+1,j+1,l);
          if (e>1) { SIG.check911(); }
       }
       n3+=getTensorProdReps_gen(*qq[i],*qq[j],TP3_LDR);
    }}

    if (vflag) printStatus(FL,"new");

    return n3;
};

template <class TQ, class TD>
void genRG_struct<TQ,TD>::printStatus(
   const char *F, int L, const char *istr
){
    unsigned n=0;

    wblog(F_L,"<i> %s status for %s: %d entries",
       istr && istr[0] ? istr:"current", STR(q), RSet.size());
    wbvector<qset<TQ> > qq(RSet.size()); {
       for (auto I=RSet.begin(); I!=RSet.end(); ++I) {
          const genRG_base<TQ,TD> &R=I->second;
          printf("   %4d : [%s] (d=%ld)%s\n",
          ++n, STR(I->first), R.Z.dim1,
          R.Sz.isEmpty() || !R.Z.dim1 ? "  *** EMPTY ***":"");
       }
    }
};

template <class TQ, class TD>
void genRG_struct<TQ,TD>::initSU2(const qset<TQ> &J) {

   if (!q.isSU2() || J.len!=1) wblog(FL,"ERR %s() invalid setting "
      "(%s with J=[%s])",FCT,STR(q),STR(J)
   );

   wbsparray<TD> Sp,Sz,S2,E;
   wbvector<double> sz; 

   genRG_base<TQ,TD> &R=RSet[J];

   CG::get_SU2mat(J[0],sz,Sp,Sz,S2,E);

   if (R.Sp.len || R.Sz.len) {
      if (R.q!=q || R.J!=J) wblog(FL,
         "ERR %s() already got initialized (%s <> %s; [%s] <> [%s])",
         FCT, STR(R.q), STR(q),
         R.J.toStr().data, STR(J));
      if (R.Sp.len!=1 || R.Sz.len!=1) wblog(FL,
         "ERR %s() already got initialized (Sp[%d], Sz[%d])",
         FCT, R.Sp.len, R.Sz.len);
      if (R.Sp[0].SIZE!=Sp.SIZE) wblog(FL,"ERR %s() "
         "already got Sp at %s (vs. %s)",
         FCT,R.Sp[0].sizeStr().data, Sp.sizeStr().data);
      if (R.Sz[0].SIZE!=Sz.SIZE) wblog(FL,"ERR %s() "
         "already got Sz at %s (vs. %s)",
         FCT,R.Sz[0].sizeStr().data, Sz.sizeStr().data
      );

      double e;
      if ((e=R.Sp[0].normDiff(Sp))>1E-12)
         wblog(FL,"ERR got SU(2) inconsistency (Sp @ %.3g)",e);
      if ((e=R.Sz[0].normDiff(Sz))>1E-12)
         wblog(FL,"ERR got SU(2) inconsistency (Sz @ %.3g)",e);
   }
   else {
      R.Sp.init(1); Sp.save2(R.Sp[0]);
      R.Sz.init(1); Sz.save2(R.Sz[0]); R.q=q; R.J=J;

      R.Z.init(sz.len,1);
      for (unsigned i=0; i<sz.len; ++i) { R.Z.data[i]=sz[i]; }
   }

#ifdef CG_CHECK_MW_PERM
   if (R.P0.len) R.P0.init();
#endif
};

template <class TQ, class TD>
int genRG_struct<TQ,TD>::get1J_gen(
   const char *F, int L, const qset<TQ> &J1,
   char flag 
){
   if (J1.len!=q.qlen()) wblog(F_L,"ERR %s() "
      "got invalid qset [%s] having %s",FCT,STR(J1),STR(q));

   qset<TQ> J2(J1.len); {  q.getDual(J1.data,J2.data); }

   int i=(J1<J2), i_=-1;
   const qset<TQ> &J  (i ? J1 : J2);
   const qset<TQ> &J_ (i ? J2 : J1);
   QSet<TQ> Q;

   Q.init1J(q,J.data); 

   CData<TQ,TD> &Z=gCS.getBUF(0,0,Q,
      (flag & TP3_LOAD) ? ((flag & TP3_ITER) | LB_GEN) : 0);

   if (Z.isEmpty() || Z.cstat.t==CGD_REF_INIT) {
      if (flag & TP3_TST) return 0;
   }
   else {
      return 0;
   }

   CG::FileLock flk; 
   if (flk.init(q,J,J_,"c1j")<=0) wblog(FL,
      "ERR %s() failed to setup c1j file lock\n%s: [%s], [%s]",
      FCT, STR(q), STR(J), STR(J_));

   const genRG_base<TQ,TD> &G1=RSet[J], &G2=RSet[J_];

      if (G1.isEmpty()) gStore.load_RSet(F_L,q,J );
      if (G2.isEmpty()) gStore.load_RSet(F_L,q,J_);

   const wbvector< wbsparray<TD> >
         &Sp1=G1.Sp, &Sp2=G2.Sp,
         &Sz1=G1.Sz, &Sz2=G2.Sz;

   if (!Sp1.len || !Sp2.len || !Sz1.len || !Sz2.len) {
      if (flag & TP3_TST) return -1;
      if (q.isSU2()) {
         if (!Sp1.len || !Sz1.len) initSU2(J );
         if (!Sp2.len || !Sz2.len) initSU2(J_);
      }
      else wblog(F_L,
        "ERR %s() Sp/Sz not yet generated for %s\n"
        "[%s]*[%s]: %d/%d; %d/%d",FCT,STR(q),STR(J1),STR(J2),
         Sp1.len, Sp2.len, Sz1.len, Sz2.len
      );
   }

   if (Sp1.len!=Sp2.len || Sz1.len!=Sz2.len) wblog(FL,
      "ERR %s()\ngot Sp/Sz inconsistency (Sp[%d/%d], Sz[%d/%d])",
      FCT, Sp1.len, Sp2.len, Sz1.len, Sz2.len);
   if (Sp1.len!=Sz1.len) wblog(FL,
      "ERR %s(): %s J =[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, STR(q), STR(J1), Sp1.len, Sz1.len);
   if (Sp2.len!=Sz2.len) wblog(FL,
      "ERR %s(): %s J_=[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, STR(q), STR(J2), Sp2.len, Sz2.len
   );

   SPIDX_T d=G1.dim(), rs=G1.Sp.len;
   wbvector< wbSparrayTD > SL(2*rs), SR(2*rs);
   wbSparrayTD HL,HR,X,Xbest,HX;

   if (d!=G2.dim()) wblog(F_L,"ERR %s() invalid dual irrep\n"
      "%s: [%s] <> [%s]",FCT,STR(q),STR(G1.J),STR(G2.J));

   for (i=0; unsigned(i)<rs; ++i) {
       SL[i].init2ref(G1.Sp[i]); G1.Sp[i].htranspose(FL,SL[i+rs]);
       SR[i].init2ref(G2.Sp[i]); G2.Sp[i].htranspose(FL,SR[i+rs]);

       Wb::MatProd(G1.Sp[i],G1.Sp[i], HL,'N','C',TD(0.5),TD(1));
       Wb::MatProd(G1.Sp[i],G1.Sp[i], HL,'C','N',TD(0.5),TD(1));

       Wb::MatProd(G2.Sp[i],G2.Sp[i], HR,'N','C',TD(0.5),TD(1));
       Wb::MatProd(G2.Sp[i],G2.Sp[i], HR,'C','N',TD(0.5),TD(1));
   }

   X.initDiag2(d).setRand_data(1.,'s'); 
   X.Normalize();

   unsigned iter=0, nK=8, nmax=d*d, eflag=0, neps=32;
   wbarray<TD> HK; wbvector<TD> EK;
   double xmin=1E99, x1=1, Ebest=1, e2=0;

   gStore.rclog(q, FL, CG_VERBOSE>6 && F,
      "%N--> %s() J=(%s) %d x %d = %d (K=%d)",FCT,STR2(J1,q), d,d,d*d,nK);

   for (; iter<nmax; ++iter) {
      i=BuildKrylovH(q, HK, nK, X, EK, HL,HR,SL,SR, x1); 

      if (i>1 || i!=i_ || iter%100==0 || (i && iter%10)) {
         gStore.rclog(q, FL, CG_VERBOSE>6 && F,
         " *  %s() %3d/%d E0=%10.4g%s", FCT,iter,d,
         double(EK[0]), i? " ok":""); 
      }; i_=i;

      if (xmin>x1) xmin=x1;
      if (i>1) {
         if (EK.len) Ebest=double(EK[0]);
         break; 
      }
      if (i) {
         if ((++eflag)>neps) { X=Xbest; break; }
      }
      else if (eflag) {
         Xbest=X; 
         if (EK.len) Ebest=double(EK[0]);
      }
   }

   if (!i && d==1) { if (EK.len) Ebest=double(EK[0]); }

   if (!EK.len || ABS(EK[0])>TD(CG_EPS1) || ABS(Ebest)>CG_EPS1) {
      wbarray<double> X_; MXPut(FL,"q1j").add(i,"i").add(d,"d")
        .add(wbarray<double>(HK),"HK")
        .add(wbvector<double>(EK),"EK")
        .add(X.toFull(X_),"X").add(Z,"Z")
        .add(HL,"HL").add(HR,"HR").add(SL,"SL").add(SR,"SR");
      wblog(FL,"ERR %s(%s) at (iter=%d/%d; i=%d)\ngot EK[0]=[%.4g %.4g]",
      FCT,STR2(J1,q),iter,nmax,i, EK.len ? double(EK[0]) : NAN, Ebest);
   }

   Z=Q; Z.cgd.init(X); Z.cstat.init(CGD_1JSY_GEN); {
      e2=Z.SkipTiny(FL);
      Z.cgd.NormSignC(FL);
   }

   gStore.rclog(q, FL, i<2 && CG_VERBOSE>6 && F,
      " *  %s() k=%2d/%d @ E0=%.3g; %.3g, %.3g (i=%d/%d)",
      FCT,iter,d,Ebest,double(xmin),SQRT(e2),i,nK
   );

   i=gStore.save_CData(0,0,Z);

   gStore.rclog(q, PFL, CG_VERBOSE>6 && F,
      "[+]  %c CBUF[%03d] 1J symbol #%05x %s ", i? 'W':'w',
      gCS.BUF.size(), Z.cstat.ID, STR((QSet<TQ>&)Z)
   );

   return 1;
};

template<class TD>
int BuildKrylovH(
   const QType &q,  
   wbarray<TD> &HK, unsigned nk, 
   wbSparrayTD &X,        
   wbvector<TD> &EK,
   const wbSparrayTD &HL,            
   const wbSparrayTD &HR,            
   const wbvector< wbSparrayTD > SL, 
   const wbvector< wbSparrayTD > SR,
   double &xmin_
){
   unsigned rval=0, ik=0, d=HL.dim();
   wbSparrayTD HX, U;
   wbarray<TD> UK;
   wbvector<TD> v;
   TD x=0; double xmin=1E99;

   if (nk<2) wblog(FL,"ERR %s() got invalid nk=%d",FCT,nk);

   if (nk>d) nk=d; 
   HK.init(nk,nk); X.Normalize();

   for (; ik<nk; ++ik) {
      GetHPsi(HX,X,HL,HR,SL,SR);
      HK(ik,ik)=X.dotProd(FL,HX);
      if (ik) {
         TD d=ABS(HK(ik-1,ik-1))+ABS(HK(ik,ik));
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
         U.contract(FL,2,v,X); HX-=X;  
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
   wbSparrayTD &HX, const wbSparrayTD &X,
   const wbSparrayTD &HL, 
   const wbSparrayTD &HR, 
   const wbvector< wbSparrayTD > SL, 
   const wbvector< wbSparrayTD > SR  
){
   wbSparrayTD Q; HX.init();
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
int genRG_struct<TQ,TD>::getTensorProdReps_gen(
   const qset<TQ> &J1_, const qset<TQ> &J2_, unsigned flag
){
   unsigned np=0, n3=0, mp=(J2_<J1_); 
   int n1=0, n2=0;

   const qset<TQ> &J1 (mp ? J2_ : J1_);
   const qset<TQ> &J2 (mp ? J1_ : J2_);
   qset<TQ> Jbar(J1.len);

   qset<TQ> J12(J1,J2);
   qset<TQ> J21(J2,J1);

   unsigned loadR = (flag & TP3_LDR ? MP3_LDR : 0);
   unsigned l=(flag & TP3_ITER);

   if (l<TP3_ITER) { l=((++flag) & TP3_ITER); if (l>2)
      wblog(FL,"WRN %s() got recursive call @ l=%d",FCT,l); }
   else { 
      wblog(FL,"ERR %s() got recursive call @ l=%d",FCT,l+1);
   }

   if (flag & TP3_LDR) flag|=TP3_LDM;

   mp=(J1!=J2 ? 2 : 1); 

   n1=gCS.find_map3(q,J12);
   if (n1<=0 && !(flag & TP3_TST)) {
      gStore.load_Std3(0,0,q,J1,J2, loadR);
      n1=gCS.find_map3(q,J12);
      if (n1>0) { n3+=n1; }
   }

   if (mp>1) {
      n2=gCS.find_map3(q,J21);
      if (n2<=0 && (n1>0 || (flag & TP3_LOAD))) {
         gStore.load_Std3(0,0,q,J2,J1, loadR);
         n2=gCS.find_map3(q,J21);
         if (n2>0) { n3+=n2; }
      }

      if (n1!=n2) { char s[32];
         wblog(FL,"WRN tensorProd() got missing *.mp3 pair (%s,%s)",
            BITS(loadR),BITS(flag));

         if (n1) snprintf(s,32,"%d output multiplet%s",n1,n1!=1?"s":"");
         else strcpy(s,"empty");
         wblog(FL,"WRN (%s,%s) %s",STR(J1),STR(J2),s);

         if (n2) snprintf(s,32,"%d output multiplet%s",n2,n2!=1?"s":"");
         else strcpy(s,"empty");
         wblog(FL,"WRN (%s,%s) %s",STR(J2),STR(J1),s);
      }

      if (n1>0 && n2>0) { return n3; }
   }
   else if (n1>0) { return n3; }

   if ((flag & TP3_TST) && n1<=0 && n2<=0) { return -1; } 

   CG::FileLock flk; 
   if (flk.init(q,J1,J2,"mp3")<=0) { wblog(FL,
      "ERR %s() failed to setup mp3 file lock\n%s: [%s], [%s]",
      FCT,STR(q),STR(J1),STR(J2));
   }

   if (n1<0) { n1=gCS.find_map3(q,J12); }
   if (mp>1) {
      if (n2<0) { n2=gCS.find_map3(q,J21); }
      if (n1!=n2) wblog(FL,"WRN %s() got n=%d/%d multiplets",FCT,n1,n2);
      if (n1>0 && n2>0) { return n3=(n1+n2); }
   }
   else if (n1>0) { return (n3=n1); }

   const genRG_base<TQ,TD> &G1=RSet[J1], &G2=RSet[J2];

   if (G1.isEmpty()) gStore.load_RSet(FL,q,J1);
   if (G2.isEmpty()) gStore.load_RSet(FL,q,J2);

   const wbvector< wbsparray<TD> >
         &Sp1=G1.Sp, &Sp2=G2.Sp,
         &Sz1=G1.Sz, &Sz2=G2.Sz;

   if (!Sp1.len && !(flag & TP3_TST)) { gRS.getR(FL,q,J1.data); }
   if (!Sp2.len && !(flag & TP3_TST)) { gRS.getR(FL,q,J2.data); }

   if (!Sp1.len || !Sp2.len || !Sz1.len || !Sz2.len) {
      if (flag & TP3_TST) return -1;
      if (q.isSU2()) {
         if (!Sp1.len || !Sz1.len) initSU2(J1);
         if (!Sp2.len || !Sz2.len) initSU2(J2);
      }
      else wblog(FL,
        "ERR %s() Sp/Sz not yet generated\n"
        "for %s ([%s]*[%s]: %d/%d; %d/%d)", FCT,
         q.toStr().data, STR(J1), STR(J2),
         Sp1.len, Sp2.len, Sz1.len, Sz2.len
      );
   }

   if (Sp1.len!=Sp2.len || Sz1.len!=Sz2.len) wblog(FL,
      "ERR %s()\ngot Sp/Sz inconsistency (Sp[%d/%d], Sz[%d/%d])",
      FCT, Sp1.len, Sp2.len, Sz1.len, Sz2.len);
   if (Sp1.len!=Sz1.len) wblog(FL,
      "ERR %s(): %s J1=[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, STR(q), STR(J1), Sp1.len, Sz1.len);
   if (Sp2.len!=Sz2.len) wblog(FL,
      "ERR %s(): %s J2=[%s]\ngot Sp/Sz inconsistency (Sp[%d] vs. Sz[%d])",
      FCT, STR(q), STR(J2), Sp2.len, Sz2.len
   );

   np=Sp1.len;

   unsigned i,it,im,m,M,d3, r=q.sub, nz=Sz1.len, nel=0,
         d1=Sp1[0].SIZE[0], d2=Sp2[0].SIZE[0], d12=d1*d2;
   int sr=0, normR3=0, isLargeD=d12>10000; double r2=0,e0=0,e2=0;
   wbvector< wbsparray<TD> > CM, Sp(np), Sz(nz);
   wbvector<genRG_base<TQ,TD> > RR;
   wbvector<double> c2eps;

   wbMatrix<unsigned> iOM;

   wbsparray<TD> X1,E1,E2;
   wbvector< wbsparray<TD> > U;
   wbMatrix<TQ> ssz,ss,sx,sz,Z;
   wbperm P, p213("2,1,3");
   wbindex I;

   const wbMatrix<double> &Z1=G1.Z, &Z2=G2.Z;
   wbvector<unsigned> dd;

   int io3=-1;

   if (G1.q!=G2.q) wblog(FL,"ERR %s() qtype inconsistency "
      "(%s; %s)",FCT,STR(G1.q),STR(G2.q));
   if (q.isNonAbelian() && (np!=r || nz!=r)) wblog(FL,
      "ERR %s()\ninvalid Sp/Sz set for %s (%d/%d)",
      FCT,STR(q),np,nz
   );
   for (i=0; i<np; ++i) {
      if (!Sp1[i].isSMatrix(FL,d1) || !Sp2[i].isSMatrix(FL,d2))
         wblog(FL,"ERR %s()\nmatrix size mismatch (%s <> %s)",
         FCT,Sp1[i].sizeStr().data,Sp2[i].sizeStr().data
      );
   }
   for (i=0; i<nz; ++i) {
      if (!Sz1[i].isSMatrix(FL,d1) || !Sz2[i].isSMatrix(FL,d2))
         wblog(FL,"ERR %s()\nmatrix size mismatch (%s <> %s)",
         FCT,Sz1[i].sizeStr().data,Sz2[i].sizeStr().data
      );
   }
   if (Z1.dim1!=d1 || Z1.dim2!=nz) wblog(FL,
      "ERR %s() invalid Sz: %dx%d <> %dx%d (%s)",
      FCT, Z1.dim1, Z1.dim2, d1,nz, STR(q));
   if (Z2.dim1!=d2 || Z2.dim2!=nz) wblog(FL,
      "ERR %s() invalid Sz: %dx%d <> %dx%d (%s)",
      FCT, Z2.dim1, Z2.dim2, d2,nz, STR(q)
   );

   E1.initIdentity(d1);
   E2.initIdentity(d2);

#ifndef WB_SKIP_ASSERT
   for (i=0; i<np; ++i) {
      if (!Sz1[i].isCompressed()) wblog(FL,"ERR %s() ",FCT); 
      if (!Sz2[i].isCompressed()) wblog(FL,"ERR %s() ",FCT);
   }
   for (i=0; i<nz; ++i) {
      if (!Sz1[i].isCompressed()) wblog(FL,"ERR %s() ",FCT); 
      if (!Sz2[i].isCompressed()) wblog(FL,"ERR %s() ",FCT);
   }
#endif

   for (i=0; i<np; ++i) {
      Sp1[i].kron(FL,E2,Sp[i]); E1.kron(FL,Sp2[i],X1);
      Sp[i].Plus(FL,X1);
   }
   for (i=0; i<nz; ++i) {
      Sz1[i].kron(FL,E2,Sz[i]); E1.kron(FL,Sz2[i],X1);
      Sz[i].Plus(FL,X1);
   }

   gStore.rclog(q, PFL, CG_VERBOSE>4 && isLargeD,
      "CG3 %s() %s (%d*%d)", 
      FCT,STR(QSet<TQ>().init2(q,J1,J2)),d1,d2);
   doflush();

   r2=CG::getSymmetryStates(FL,q,Sp,Sz,U,dd,RR,c2eps,&iOM);

   if (isLargeD) { gStore.rclog(q, PFL, CG_VERBOSE>4,
      "CG3 %s() %s (%dx%d)[%d]", 
      FCT,STR(QSet<TQ>().init2(q,J1,J2)),d1,d2,RR.len);
   }; doflush();

   for (it=0; it<RR.len; ++it) {
      const genRG_base<TQ,TD> &R=RR[it];
      const qset<TQ> J(R.J);

      CData<TQ,TD> X;
      genRG_base<TQ,TD> G, &G0=RSet[J]; d3=dd[it];
      int ioR=-1, ioC=-1, saveR=0, saveC=0, xflag=0;

      X.init3(q,d1,d2,d3); X.cstat.init(CGD_STD3); 
      X.qs.init(J1,J2,J);

      if (!U[it].hasSize(d1*d2,d3)) wblog(FL,
         "ERR %s() size mismatch (%s <> %dx%dx%d)",
         FCT,U[it].sizeStr().data,d1,d2,d3
      );

      U[it].Reshape(d1,d2,d3); 

      X.cgd.wbsparray<TD>::init(FL,d1,d2,d3, U[it]);

      if (CG::signFirstVal(FL,X.cgd.D.data, X.cgd.D.len)<0)
      wblog(FL,"ERR %s() wrong sign!",FCT); 

      X.cgd.NormSignC();

      G.q=G1.q; G.J=J; G.Z=R.Z; G.Sp=R.Sp; G.Sz=R.Sz;
#ifdef CG_CHECK_MW_PERM
      G.P0=R.P0;
#endif

      G.SkipTiny(); 

      try { r2+=(e2=checkCommRel(FL,R)); e2=sqrt(e2); }
      catch (...) {
         MXPut(FL,"I_CR").add(J1,"J1").add(J2,"J2").add(J,"J")
         .add(G1,"G1").add(G2,"G2").add(G,"G");
         wblog(FL,"ERR %s() inconsistency",FCT);
      };

      if (G0.isEmpty()) gStore.load_RSet(0,0,q,J,RCL_V2); 
      if (!G0.Sp.isEmpty() || !G0.Sz.isEmpty()) {
         double e=G.normDiff(0,0,G0,1E-10);
         if (e>1E-10) {
            MXPut(FL).add(J1,"J1").add(J2,"J2").add(J,"J")
              .add(G1,"G1").add(G2,"G2").add(G,"G").add(G0,"G0");
            wblog(FL,"ERR %s()\n"
              "inconsistency with earlier data (e=%.3g)",FCT,e);
         }

         e0=sqrt(checkCommRel(FL,G0));
         if (e2<=0.5*e0) { saveR=2;
            gStore.rclog(q, PFL, CG_VERBOSE>6,
               "NB! replacing %s RStore (%s) [e=%.3g->%.3g]",
               STR(q), R.J.toStr(q).data, e0, e2
            );
         }
      }
      else {
         if (q.isSU2()) G.compareStdSU2(FL);
         saveR=1;
      }

      if (saveR) {
         char sflag=(saveR>1 && G0.istr.data && G0.istr.data[0]);
         unsigned l=0, n=64 + (sflag ? strlen(G0.istr.data) : 0);
         char istr[n];

         if (sflag  ) { l=snprintf(istr,n,"%s; ",G0.istr.data); } else
         if (saveR>1) { l=snprintf(istr,n,"e=%.3g; ", e0); }
         l+=snprintf(istr+l,n-l,"(%s,%s) @ e=%.3g",
             J1.toStr(q).data, J2.toStr(q).data, e2);
         if (l>=n) wblog(FL,"WRN %s() "
            "string out of bounds (%d/%d)\n%s",FCT,l,n,istr);

         G.save2(G0);
         G0.istr=istr; G0.err=e2;

         ioR=gStore.save_RSet(0,0,G0);
      }

      if (G0.q.validType()) {
         gCS.getIdentityC(FL,q,J.data);
      }

      CData<TQ,TD> &Cb=gCS.getBUF(0,0,(QSet<TQ>&)X,LB_UPD);

      if (!Cb.cgd.isEmpty()) {
         if (Cb.isRefInit()) wblog(FL,
            "ERR %s() got REF data\n%s",FCT,STR(Cb));

         if (Cb.t!=X.t || Cb.qs!=X.qs || Cb.qdir!=X.qdir) wblog(FL,
            "ERR %s() got CData inconsistency\n[%s] <> [%s]",
            FCT, Cb.QStr().data, X.QStr().data
         );
      }

      CRef<TQ> &R3=gCS.map3[q][J12][J];
      if (iOM.dim1)
           { im=iOM(it,0); M=iOM(it,1); }
      else { im=0; M=1; }

      if ((normR3=(!R3.cgw))) {
         R3.initBase(&Cb,M,&X);
      }
      else if (!R3.wdim_is(M)) wblog(FL,
        "ERR %s() size mismatch map3->cgw %s /%d",FCT,SSTR(R3.cgw),M);

      if (mp>1) {
         CData<TQ,TD> X2;
         X.permute(X2,p213); 
         sr=CG::signFirstVal(FL,X2.cgd.D.data, X2.cgd.D.len);
      }

      nel+=X.cgd.numel(); 

      m=Cb.getOM(FL); 
      if (im<m) { xflag|=1; }

      if (xflag) {
         if (M>1) { unsigned m2; wbvector<RTD> w;
            Cb.Project(FL,X.cgd,w);  
            R3.cgw.setCol(im,w);     

            m2=Cb.getOM();
            if (m2<m || m2>m+1 || m2>M) { 
               MXPut(FL,"Iq","base")
                 .add(X,"X").add(Cb,"Cb").add(R3,"R3")
                 .add(im,"im").add(m,"m").add(m2,"m2").add(M,"M");
               wblog(FL,"ERR %s() OM out of bounds (%d/%d/%d / %s)",
                  FCT,m,m2,M, SSTR(R3.cgw) );
            }
if (im>1) {
   MXPut(FL,"Iq","base").add(X,"X").add(Cb,"Cb").add(R3,"R3")
     .add(im,"im").add(m,"m").add(m2,"m2").add(M,"M");
   wblog(FL,"ERR %s() check this",FCT);
}
            if (m<m2) { m=m2; n3+=mp; saveC|=1; }
         }
         else {
            double e=X.normDiff(FL,Cb,1E-10);
            if (e<1E-10) { m=1;
               if (CG_VERBOSE>6) wblog(FL, 
                  "ok. %s() %s @ e=%.3g",FCT,STR(Cb),e);

               if (Cb.cstat!=CGD_STD3_X) {
                  X.checkAdditivityZ(FL,G1.Z,G2.Z,R.Z);

                  X.cstat.ID   =Cb.cstat.ID;
                  X.cstat.ctime=Cb.cstat.ctime;

                  X.save2(Cb);
               }

               if (im>=M) wblog(FL,
                  "ERR %s() index out of bounds (%d/%d)",FCT,im,M);
               R3.cgw_check_std3(FL,d3,im);
            }
            else {
               wblog(FL,"TST %s() 0x%lX  %d: %d %d",FCT,iOM.data,it+1,
                  Cb.getOM(FL), it<iOM.dim1 ? iOM(it,1):-1);
               MXPut(FL).add(J1,"J1").add(J2,"J2")
                  .add(Cb,"Cb").addP(Cb.cgd.toMX(),"C0")
                  .add(X, "X" ).addP(X .cgd.toMX(),"C" )
                  .add(R,"R").add(Sp,"Sp").add(Sz,"Sz").add(U,"U")
                  .add(iOM,"iOM").add(dd,"dd")
                  .add(it+1,"it").add(d3,"d3").add(e,"e");
               wblog(FL,"ERR overwriting cgd: %s <> %dx%dx%d [%s] @ %.3g",
               SSTR(X.cgd),d1,d2,d3,STR(J),e);
            }
         }
      }
      else { n3+=mp; saveC|=2;
         if (im>m || m>=R3.wdim(FL) || (im && Cb.cstat.t!=CGD_STD3))
            wblog(FL,"ERR %s() im=%d/%d (it=%d)\n%s\n%s",
            FCT, im,m,it, SSTR(R3), STR(Cb)
         );

         if (im) { Cb.AddMultiplicity(FL,X); }
         else {
            X.checkAdditivityZ(FL,G1.Z,G2.Z,R.Z);
            X.save2(Cb); 
         }

         if (im>=R3.wdim(FL)) wblog(FL,
            "ERR %s() index out of bounds (%d: %s)",FCT,im,SSTR(R3.cgw));
         m=Cb.getOM(FL);
         if (m!=im+1) wblog(FL,"ERR %s() m=%d/%d",FCT,m,im+1);

         R3.cgw_check_std3(FL,d3,im);
      }

      if (xflag && (Cb.cstat!=CGD_STD3 && Cb.cstat!=CGD_STD3_X)) {
         char s[128];
         snprintf(s,128,"converting %s : %s -> %s",
            STR(Cb), Cb.cstat.tstr(), CGD_TYPE_STR[CGD_STD3_X]);
         wblog(FL,Cb.cstat==CGD_FROM_DEC ? "CHK %s":"WRN %s",s);

         Cb.cstat=CGD_STD3_X; saveC|=4; 
      }

      if (normR3) {
          R3.NormStd(FL,3); 
      }

      if (mp>1) {
         CRef<TQ> &R3p=gCS.map3[q][J21][J];

         if (R3p.cgw) { if (!R3p.cgw.sameSize(R3.cgw)) wblog(FL,
            "ERR %s() CRef size mismatch (%s/%s)",FCT,SSTR(R3p),SSTR(R3)); }
         else {
            R3.permute(p213,R3p);
         }

         R3p.cgw.setCol(im, R3.cgw.col(im), sr<0 ? -1:1);

if (xflag && im>1) { 
   MXPut(FL,"Im","base").add(*R3p.cgb,"Cb").add(R3p,"R3p").add(R3,"R3")
    .add(sr,"sr").add(xflag,"xflag").add(im,"im").add(m,"m");
   wblog(FL,"ERR %s() check this\n%s",FCT,STR(R3));
}
      }

      if (m!=R3.wdim()) { continue; }

      if (saveC){
         ioC=gStore.save_CData(0,0,Cb);

         if (J.allZero()) { 
            QSet<TQ> Q1(Cb);
            if (!J.len || Q1.qdir.len!=3) wblog(FL,"ERR %s()",FCT);
            if (m!=1) wblog(FL,"ERR %s() got 1J symbol @ om=%d",FCT,m);
            Q1.qdir.len=2;
            Q1.qs.len=2*J.len;

            CData<TQ,TD> &S1=gCS.getBUF(0,0,Q1,LB_LOAD);
            if (S1.isEmpty()) {
               Cb.reduceto1JSymbol(FL,S1);
               int io1=gStore.save_CData(0,0,S1); 

               gStore.rclog(q, PFL, CG_VERBOSE>6,
                  "[+]  %c CBUF[%03d] 1J symbol #%05x %s", io1? 'W':'w',
                  gCS.BUF.size(), S1.cstat.ID, STR(Q1)
               );
            }
            else {
               CData<TQ,TD> Z;
               Cb.reduceto1JSymbol(FL,Z); double e=S1.normDiff(FL,Z);

               if (e>CG_EPS2 || CG_VERBOSE>2) { char s[256];
                  snprintf(s,256,"tensorProd() got existing 1J symbol "
                     "@ %.3g\n   %s %s\n<> %s %s",e, STR(Z),
                     STR2(Z.cstat,'v'), STR(S1), STR2(S1.cstat,'v'));
                  if (e>CG_EPS2) { MXPut(FL,"x1j")
                     .add(Cb,"Cb").add(S1,"Z1").add(Z,"Z2").add(e,"e");
                     wblog(FL,"ERR %s",s);
                  }
                  gStore.rclog(Cb.t, FL, CG_VERBOSE>6,"ok. %s",s);
               }
            }
         }
      }

      if (1) { 
         unsigned l=16; char s[l];

         if (it==0) {
            gStore.rclog(q,PFL,CG_VERBOSE>6,
               "CG3 (%s,%s) -> %d ireps @ %d*%d %s",STR2(J1,q),
               STR2(J2,q),RR.len, d1,d2, mp>1 ? "(+reverse)":""
            );
         }

         if (iOM.dim1 && iOM(it,1)>1)
              { snprintf(s,l,"%4d x%d",d3,iOM(it,1)); }
         else { snprintf(s,l,"d=%d",d3); }

         gStore.rclog(q,PFL,CG_VERBOSE>6,
            "[+] %c%c map3[%03d] #%05x %-6s (%s) %s",
            ioR<0? ' ': (ioR?'W':'w'), ioC<0? ' ': (ioC?'W':'w'),
            gCS.map3[q].size(), Cb.cstat.ID, STR(q),
            Cb.QStr().data, s
         );
         doflush();
      }

   } 

   if (nel && nel!=d12*d12) wblog(FL, 
      "ERR %s() state space inconsistency (%d/%dx%d)",FCT,nel,d12,d12);

   io3 =   gStore.save_Std3(0,0,q,J12,&c2eps); if (mp>1) { 
   io3+=10*gStore.save_Std3(0,0,q,J21,&c2eps); }

   if (CG_VERBOSE>8 || isLargeD || q.isLargeD(dd)) {
      if (dd.sum()!=d1*d2) wblog(FL,"ERR %s() "
         "got size inconsistency %d ~= %dx%d",FCT,dd.sum(),d1,d2);

      gStore.rclog(q, PFL, CG_VERBOSE>6 && (isLargeD || q.isLargeD(dd)),
         "%s  %c %s %dx%d = %d, CR @ %.3g", io3%2 ? "WRN":"[+]",
         io3%2 ? 'W':'w', q.type==QTYPE_UNKNOWN ? "general" : STR(q),
         d1, d2, d1*d2, sqrt(r2)
      );
      fflush(0); doflush();
   }

   gCS.reduceMemUsage(FL); 
   gStore.rclog(q,0,0, CG_VERBOSE>6,"\n");

   if ((flag & TP3_ITER)<=1) {
   for (it=0; it<RR.len; ++it) { q.getDual(RR[it].J.data,Jbar.data);
      if (gStore.load_RSet(0,1,q,Jbar)>0) {
         if (RR[it].Z.dim1<(1<<13)) {  
            get1J_gen(FL, RR[it].J,flag); }
      }  
   }}

   return n3;
};

template <class TQ, class TD>
double CG::getSymmetryStates(const char *F, int L, const QType &q, 
   const wbvector< wbSparrayTD > &Sp,
   const wbvector< wbSparrayTD > &Sz,
   wbvector< wbSparrayTD > &UK, 
   wbvector<unsigned> &dd, 
   wbvector<genRG_base<TQ,TD> > &RR,
   wbvector<double> &c2eps,  
   wbMatrix<unsigned> *iOM,  
   char vflag
){
   unsigned i,j,ip,i0=0, D,d0=0, nu=0, m,found,
        it=0, nt, r=q.sub, np=Sp.len, nz=Sz.len;
   double c2=0, r2=0; int largeD;

   wbvector< wbSparrayTD* > uk;
   wbSparrayTD U,x1,x2,v0,vi;
   wbMatrix<TQ> z2,JJ;
   wbvector<widx_t> dJ;
   wbvector<TD> sz;
   wbperm P;

   TD x,x12,vi2;

   TD eps=1E-8, eps2=1E-10;

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
      FCT, q.type==QTYPE_UNKNOWN ? "general" : STR(q), D);

   U.init(D,0); 

   dd.init(D);  
   uk.init(D);

   for (it=0; it<D; ++it) {
      v0.initz(D,1,1); c2=0;
      for (; i0<D; ++i0) {
         v0.setRec(0, i0,0, 1.); if (!i0 && !it) break;

         Wb::MatProd(U,v0,x1,'C'); x=x1.norm();
         if (ABS(x-1)<eps) {
            if (ABS(x-1)>eps2) wblog(FL,"WRN %s() "
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

      wbSparrayTD V(v0); 
      wbSparrayTD Vi;    
      m=0; found=1;

      while (found) { found=0;
      for (ip=0; ip<Sp.len; ++ip) {

          Wb::MatProd(Sp[ip],v0,vi,'C'); 

          vi2=vi.norm2(); x=SQRT(vi2/vi.SIZE[1]);
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
      wbSparrayTD &V=UK[i]; d0=dd[i];

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
            "(%s,%d; %d: %d/%d)", q==QTYPE_UNKNOWN ? 
            "*":q.toStr('t').data, P.data[0]+1,i+1,d0,D);

          MXPut X(FL,"a"); X.add(Sp,"Sp").add(Sz,"Sz").add(UK,"UK")
           .add(dd,"dd"); if (iOM) X.add(*iOM,"M");
          X.add(D,"D").add(V,"V").add(R.Z,"Z").add(P,"P"); X.put();
      }

      if (!P.isIdentityPerm()) {
          R.Z.recPermute(P);
          V.ColPermute(P); 

          for (j=0; j<np; ++j) R.Sp[j].MatPermute(P);
          for (j=0; j<nz; ++j) R.Sz[j].MatPermute(P);

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
         wblog(FL,"NB! %s() "
           "got outer multiplicity (OM<=%d)",FCT,dJ.max());
      }
      return r2;
   }
   else if (iOM) iOM->init();

   return r2;
};

template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::SetupSym(
   const char *F, int L,
   const QType &q0, qset<TQ> *qs_
){
   q0.validType(F_L);

   if (qdef.isEmpty()) { int ep=0; Wb::LogException e;
     #pragma omp critical (within_initBase) 
      if (qdef.isEmpty()) { 
         try {
            gStore.setupDirs(F_L); q=q0;
            switch (q.type) {
               case QTYPE_SUN: Setup_SUN(F_L,&qdef); break;
               case QTYPE_SpN: Setup_SpN(F_L,&qdef); break;
               case QTYPE_SON: Setup_SON(F_L,&qdef); break;
               case QTYPE_SEN: Setup_SEN(F_L,&qdef); break;
               default: wblog(FL,
               "ERR %s() type '%s' not implemented yet",FCT,STR(q));
            }
         }
         catch (Wb::LogException &e_) { e=e_; }
         catch (...) { ++ep; } 
      }

      if (e.type) { ep+=100; e.init(); }
      if (ep) {
         sprintf(str,"%s() for %s [e=%d]",FCT,STR(q),ep);

         if (Wb::SigHandler::check911_()) wblog(FL, 
            "WRN initial RCStore population interrupted\n%s",FCT,str);
         else {
            qdef.init(); q.init(); 
            wblog(FL,"ERR %s",str);
         }
      }
   }

   if (qs_) { (*qs_)=qdef; }
   return *this;
};

template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::Setup_SUN(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i,j,l, nrep=0, r=q.sub, N=r+1, r2=N/2; 
   TD z;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (N<2 || N>10) wblog(FL,"ERR %s() invalid SU(%d)",FCT,N);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (i=1; i<N; ++i) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];

       P.initz(N,N,1); P.setRec(0, i-1,i, 1.);

       Z.initDiag(N); Z[i]=-int(i); for (j=0; j<i; ++j) Z[j]=1;

   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(N,r); X.initIdentity(N);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=N) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d)",FCT,STR(q),R.Sz[i].isDiag(),sz.len,N);
      for (j=0; j<N; ++j) R.Z(j,i)=double(sz[j]);
   }

   CG::findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];
   char s_[64];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() already got existing RSet (%s: %d,%d)",
      FCT,STR(q),R0.Sp.len,R0.Sz.len
   );
   else if (gStore.load_RSet(0,1,q,R.J)>0) {
      double e=R.normDiff(0,0,R0,1E-10);
      if (e>1E-10) wblog(FL,
         "ERR %s() inconsistent %s generators (e=%.3g)",FCT,STR(q),e);
      wbstring fs; gStore.get_file_name(F_L,fs,R,"rep",RC_SAVE);
   }
   else {
      R0=R; 

      R0.err=checkCommRel(FL,R0); 
         snprintf(s_,64,"%s @ e=%.3g",FCT,R0.err);
      R0.istr=s_;

      int ioR=gStore.save_RSet(FL,R0,'q');
      if (!ioR) { ++nrep; }

      gStore.rclog(q,PF_L,CG_VERBOSE>5 || (nrep && CG_VERBOSE>2),
         "%s defining irep (%s) for %s",
         ioR ? (ioR>1? "ok.":"[+] W/") : "[+] w/",
         QSet<TQ>().init1(q,R0.J.data).QStr().data, STR(R0.q)
      ); doflush();

   }

   if (r>1) {
      R.Z*=(-1); R.Z.FlipRecs();

      for (i=0; i<r; ++i) {
         wbsparray<TD> &P=R.Sp[i], &Z=R.Sz[i];
         P.setRec(0, r-i-1,r-i, 1.); 
         for (j=0; j<r2; ++j) { z=Z[j]; Z[j]=-Z[r-j]; Z[r-j]=-z; }
         if (N%2) Z[j]=-Z[j]; 
      }

      CG::findMaxWeight(q,R.Z,&R.J);
      R.err=checkCommRel(FL,R); 

      genRG_base<TQ,TD> &Rd = RSet[R.J];

      if (Rd.Sp.len || Rd.Sz.len) wblog(FL,"WRN %s() "
         "already got existing RSet (%d,%d)",FCT,Rd.Sp.len,Rd.Sz.len);
      else if (gStore.load_RSet(0,1,q,R.J)>0) {
         double e=R.normDiff(0,0,Rd,1E-10); if (e>1E-10) wblog(FL,
        "ERR %s() inconsistent %s generators (e=%.3g)",FCT,STR(q),e);
      }
      else {
         R.save2(Rd);

         snprintf(s_,64,"%s @ e=%.3g",FCT,Rd.err);
         Rd.istr=s_;

         int ioR=gStore.save_RSet(FL,Rd,'q');
         if (!ioR) { ++nrep; }

         gStore.rclog(q,PF_L,CG_VERBOSE>5 || (nrep && CG_VERBOSE>2),
            "%s dual to defining irep for %s is (%s)",
            ioR ? (ioR>1 ? "ok." : "[+] W/") : "[+] w/",
            STR(Rd.q), QSet<TQ>().init1(q,Rd.J.data).QStr().data
         ); doflush();
      }
   }

   unsigned m=0;

   unsigned dmax=pow(10,r);

   j=nrep; 

   if (nrep) { 
      if (N==3) { nrep=3; }
      else if (N<=4) { nrep=3; } 
      else { nrep=3; } 

      wblog_RC_1st_build(FLF,nrep);
   }
   else { nrep=2; }

   try { for (i=0; i<nrep; ++i) {
      if (int(l=genTensorProds(dmax,'i'))>0) m+=l; 
   }}
   catch (...) { Wb::SigHandler::check911_(FL); } 

   if ((m && CG_VERBOSE>7) || (j && CG_VERBOSE>2)) wblog(PFL,
      " *  generated %d CGTs for %s [%d passes]",m,STR(q),nrep);

   return *this;
};

template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::Setup_SpN(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i,j, nrep=0, r=q.sub, D=2*r;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (D<2 || D>10) wblog(FL,"ERR %s() invalid Sp(%d)",FCT,D);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (i=1; i<=r; i++) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];
       P.initz(D,D, i<r ? 2:1); Z.initDiag(D);

       if (i<r) {
          P.setRec(0, i-1, i, 1.);
          P.setRec(1, 2*r-i-1, 2*r-i, +1.);

          Z[i]=-int(i); for (j=0; j<i; ++j) Z[j]=1;
       }
       else {
          P.setRec(0, i-1, i, 1.); 
          for (j=0; j<r; ++j) { Z[j]=1; }  
       }

       for (j=0; j<r; ++j) Z[r+j]=-Z[r-1-j];
   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(D,r); X.initIdentity(D);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=D) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d)",FCT,STR(q),R.Sz[i].isDiag(),sz.len,D);
      for (j=0; j<D; ++j) R.Z(j,i)=double(sz[j]); 
   }

   R.Sort();

   CG::findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];
   char s_[64];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() already got existing RSet (%s: %d,%d)",
      FCT,STR(q),R0.Sp.len,R0.Sz.len
   );
   else if (gStore.load_RSet(0,1,q,R.J)>0) {
      double e=R.normDiff(0,0,R0,1E-10);
      if (e>1E-10) wblog(FL,
         "ERR %s() inconsistent %s generators (e=%.3g)",FCT,STR(q),e);
      wbstring fs; gStore.get_file_name(F_L,fs,R,"rep",RC_SAVE);
   }
   else {
      R.save2(R0);

      R0.err=checkCommRel(FL,R0); 
         snprintf(s_,64,"%s @ e=%.3g",FCT,R0.err);
      R0.istr=s_;

      int ioR=gStore.save_RSet(FL,R0,'q');
      if (!ioR) { ++nrep; }

      gStore.rclog(q,PF_L,CG_VERBOSE>5 || (nrep && CG_VERBOSE>2),
         "%s defining irep for %s is (%s)",
         ioR ? (ioR>1 ? "ok." : "[+]  W") : "[+]  w",
         STR(R0.q), QSet<TQ>().init1(q,R0.J.data).QStr().data
      ); doflush();

   }

   unsigned m=0; int l;

   unsigned dmax=pow(10,r); 

   j=nrep; 

   if (nrep) { 
      if (r==2) { nrep=2; } 
      else      { nrep=3; } 
      wblog_RC_1st_build(FLF,nrep);
   }
   else { nrep=2; }

   try { for (i=0; i<nrep; ++i) {
      l=genTensorProds(dmax,'i');
      if (int(l)>0) { m+=l;
         if (CG_VERBOSE>2) wblog(PFL,"... %4d/%d -> %3d CGTs",i+1,nrep,l);
      }
   }}
   catch (...) { Wb::SigHandler::check911_(FL); } 

   if ((m && CG_VERBOSE>6) || (j && CG_VERBOSE>2)) wblog(PFL,
      " *  generated %d CGTs for %s [%d passes]",m,STR(q),nrep);

   return *this;
};

template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::Setup_SON(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i=1, j, i2=0, nrep=0, r=q.sub, D=2*r+1;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (D<3 || D>12) wblog(FL,"ERR %s() invalid SO(%d)",FCT,D);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (; i<=r; ++i, i2+=2) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];
       P.initz(D,D,2); Z.initDiag(D);

       Z[i2]=1; Z[i2+1]=-1;

       if (i<r) {
          P.setRec(0, i2+2, i2,   1.); 
          P.setRec(1, i2+1, i2+3, 1.); 
       }
       else {
          P.setRec(0, i2+2, 1,  1.);
          P.setRec(1, 0, i2+2,  1.);
       }
   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(D,r); X.initIdentity(D);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=D) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d)",FCT,STR(q),R.Sz[i].isDiag(),sz.len,D);
      for (j=0; j<D; ++j) R.Z(j,i)=double(sz[j]); 
   }

   R.Sort();

   CG::findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];
   char s_[64];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() possibly already called earlier (%d)",
       FCT,R0.Sp.len,R0.Sz.len);
   R.save2(R0);

   R0.err=checkCommRel(FL,R0); 
      snprintf(s_,64,"%s @ e=%.3g",FCT,R0.err);
   R0.istr=s_;

   int ioR=gStore.save_RSet(FL,R0,'q');

   if (!ioR) { ++nrep; }
   i=(ioR ? 0 : 1);
      if (nrep && CG_VERBOSE>2) { i|=2; } else
      if (CG_VERBOSE>5) { i|=4; }

   if (i) {
      unsigned n=64; char s[n];
      snprintf(s,n,"%s defining irep for %s is (%s)",
         ioR ? (ioR>1 ? "ok." : "[+]  W") : "[+]  w",
         STR(R0.q), QSet<TQ>().init1(q,R0.J.data).QStr().data);

      if (i&1) {
         gStore.rclog(q,PF_L,CG_VERBOSE>5 || (nrep && CG_VERBOSE>2),s);
      }  else { wblog(PF_L,s); }
      doflush();
   }

   unsigned m=0; int l;

   unsigned dmax=pow(10,r); 

   j=nrep; 

   if (nrep) { 
      nrep=(r<3 ? 3 : 2);
      wblog_RC_1st_build(FLF,nrep);
   }
   else { nrep=2; }

   try { for (i=0; i<nrep; ++i) {
      if (int(l=genTensorProds(dmax,'i'))>0) m+=l;
   }}
   catch (...) { Wb::SigHandler::check911_(FL); } 

   if ((m && CG_VERBOSE>6) || (j && CG_VERBOSE>2)) wblog(PFL,
      " *  generated %d CGTs for %s [%d passes]",m,STR(q),nrep);

   if (CG_VERBOSE && ioR<2) {
      wblog(PF_L,"<i> CGC spaces for %s: J=[%s]",
         R0.q.toStr().data, STR(R0.J));
      snprintf(s_,64,"ISO%d",D); put(F_L,s_);
   }

   return *this;
};

template <class TQ, class TD>
genRG_struct<TQ,TD>& genRG_struct<TQ,TD>::Setup_SEN(
   const char *F, int L, qset<TQ> *qs
){
   unsigned i=1, j, i2=0, nrep=0, r=q.sub, D=2*r;

   genRG_base<TQ,TD> R;
   cgsparray X;

   if (D<3 || D>12) wblog(FL,"ERR %s() invalid SO(%d)",FCT,D);

   R.q=q; R.Sp.init(r); R.Sz.init(r);

   for (; i<=r; ++i, i2+=2) {
       wbsparray<TD> &P=R.Sp[i-1], &Z=R.Sz[i-1];
       P.initz(D,D,2); Z.initDiag(D);

       Z[i2]=1; Z[i2+1]=-1;

       if (i<r) { 
          P.setRec(0, i2+2, i2,   1.); 
          P.setRec(1, i2+1, i2+3, 1.); 
       }
       else {
          P.setRec(0, 2, 1, 1.);
          P.setRec(1, 0, 3, 1.);
       }
   }

   initCommRel(FL,R);

   R.J.init(r); R.Z.init(D,r); X.initIdentity(D);

   for (i=0; i<r; ++i) {
      const wbvector<TD> &sz=R.Sz[i].D;
      if (!R.Sz[i].isDiag() || sz.len!=D) wblog(FL,"ERR %s(%s) "
         "%d (%d/%d)",FCT,STR(q),R.Sz[i].isDiag(),sz.len,D);
      for (j=0; j<D; ++j) R.Z(j,i)=double(sz[j]); 
   }

   R.Sort();

   CG::findMaxWeight(q,R.Z,&R.J); if (qs) (*qs)=(R.J);

   genRG_base<TQ,TD> &R0 = RSet[R.J];
   char s_[64];

   if (R0.Sp.len || R0.Sz.len) wblog(FL,
      "WRN %s() possibly already called earlier (%d)",
       FCT,R0.Sp.len,R0.Sz.len);
   R.save2(R0);

   R0.err=checkCommRel(FL,R0); 
      snprintf(s_,64,"%s @ e=%.3g",FCT,R0.err);
   R0.istr=s_;

   int ioR=gStore.save_RSet(FL,R0,'q');

   if (!ioR) { ++nrep; }
   i=(ioR ? 0 : 1);
      if (nrep && CG_VERBOSE>2) { i|=2; } else
      if (CG_VERBOSE>5) { i|=4; }

   if (i) {
      unsigned n=64; char s[n];
      snprintf(s,n,"%s defining irep for %s is (%s)",
         ioR ? (ioR>1 ? "ok." : "[+]  W") : "[+]  w",
         STR(R0.q), QSet<TQ>().init1(q,R0.J.data).QStr().data
      );

      if (i&1) {
         gStore.rclog(q,PF_L,CG_VERBOSE>5 || (nrep && CG_VERBOSE>2),s);
      }  else { wblog(PF_L,s); }
      doflush();
   }

   if (CG_VERBOSE>5) {
      MXPut(FL,"I0").add(R0,"R0").add(*this,"all");
   }

   unsigned m=0; int l;

   unsigned dmax=pow(10,r); 

   j=nrep; 

   if (nrep) { 
      nrep=(r<3 ? 3 : 2);
      wblog_RC_1st_build(FLF,nrep);
   }
   else { nrep=2; }

   try { for (i=0; i<nrep; ++i) {
      if (int(l=genTensorProds(dmax,'i'))>0) m+=l;
   }}
   catch (...) { Wb::SigHandler::check911_(FL); } 

   if ((m && CG_VERBOSE>6) || (j && CG_VERBOSE>2)) wblog(PFL,
      " *  generated %d CGTs for %s [%d passes]",m,STR(q),nrep);

   if (CG_VERBOSE>5 && ioR<2) {
      wblog(PF_L,"<i> CGC spaces for %s: J=[%s]",
         R0.q.toStr().data, STR(R0.J));
      snprintf(s_,64,"I_SO%d",D); put(F_L,s_);
   }

   return *this;
};

#endif

