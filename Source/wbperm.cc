/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace permute routines
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

#ifndef __WB_PERM_CC__
#define __WB_PERM_CC__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //
// flattens result -> ensures inv=0

wbperm& wbperm::init(const wbperm &P, char iflag, unsigned r) {

   if (abs(iflag)>1 && iflag!='i') wblog(FL,
      "WRN %s() got iflag=%s",FCT,cSTR(iflag));
   if (abs(P.inv)>1) wblog(FL,"WRN %s() got '%s'",FCT,STR(P));

   if (&P==this) {
      if (iflag) { ++inv; }
      if ((inv%=2)) { flatten(); }
   }
   else if (iflag) {
      if (P.inv)
           { (*this)=P; inv=0; }
      else { P.invert(*this); }  
   }
   else {
      (*this)=P; if (P.inv) {
         if (len) { P.invert_data(data); }
         fac=1./fac; 
         inv=0; 
      }
   }

   if (int(r)>=0) { Extend(r); } 

   return *this;
};

wbperm& wbperm::Invert(char cflag) {

   if ((inv%=2)) { inv=0; }
   else if (!fac) { wblog(FL,"ERR %s() got %s",FCT,STR(*this)); }
   else {
      invert_data(); 
      fac=1./fac;    
   }

   if (cflag) { Wb::conj_add_z2(conj,cflag); }

   return *this;
};

wbperm& wbperm::invert(wbperm &iP) const {

   if (&iP==this) { return iP.Invert(); }
   iP=(*this); 

   if ((iP.inv%=2)) { iP.inv=0; }
   else if (!fac) { wblog(FL,"ERR %s() got %s",FCT,STR(*this)); }
   else {
      this->invert_data(iP.data);
      iP.fac=1./fac;
   }
   return iP;
};

wbperm& wbperm::Extend(wperm_t r) {

   if (r!=len) {
      if (int(r)<0) wblog(FL,"ERR %s() requesting r=%d/%d",FCT,r,len);

      if (r>len) {
         wbperm P(r);  
         if (len) { memcpy(P.data,data,len*sizeof(wperm_t)); }
         P.WBPERM::save2(*this); 
      }
      else {
         for (unsigned i=r; i<len; ++i) {  
            if (data[i]!=i) wblog(FL,
            "ERR %s() invalid r=%d having %s",FCT,r,STR(*this));
         }
         WBPERM::Shorten2(r);
      }
   }
   return *this;
};

wbperm& wbperm::Complete(wperm_t l) {
   if (l<len) {
      if (l) {
         unsigned i=0; wbvector<char> mark(len);
         for (; i<l; ++i) {
            if (data[i]>=len) wblog(FL,"ERR %s() invalid partial perm\n"
               "(index out of bounds i=%d: %d/%d)",FCT,i+1,data[i],len);
            if (++mark[data[i]]!=1) wblog(FL,"ERR %s() "
               "invalid partial perm (non-unique: %s)",FCT,STR(*this));
         }
         for (i=0; i<len; ++i) { if (!mark[i]) { data[l++]=i; }}
      }
      else { init(l=len); }
   }
   else if (l>len) wblog(FL,
     "ERR %s() invalid usage (l=%d/%d)",FCT,l,len);

   return *this;
};

wbperm& wbperm::Complete1() {
   unsigned i,j;
   MVEC mark(len); wbperm P;
   P.init_bare(len); 

   for (i=0; i<len; ++i) {
      if (data[i]) { j=data[i]-1; ++mark[i];
         if (j>=len) wblog(FL,"ERR %s() invalid partial perm\n"
            "(index out of bounds i=%d: %d/%d)",FCT,i+1,data[i],len);
         if (P[j]) wblog(FL,"ERR %s() "
            "invalid partial perm (non-unique: %s)",FCT,STR(*this));
         P[j]=i+1; 
      }
   }

   for (i=-1, j=0; j<len; ++j) {
      if (P[j])
           { --P[j]; } 
      else { for (++i; i<len; ++i) { if (!mark[i]) { P[j]=i; break; }}}
   }
   if (i>=len) wblog(FL,"ERR %s() failed to complete perm "
      "(%s; i=%d/%d, q=%d)",FCT,STR(P),i,len,P.isValidPerm());
   return P.save2(*this);
};

int wbperm::initStr(
   const char *F, int L, const char *s, wperm_t offset) {

   if (!s) { init(); return -1; }     

   while (isspace(*s) && *s) { ++s; } 
   if (!*s){ init(); return 0; }      

   unsigned i=0, n=strlen(s);

   inv=conj=0; fac=1;

   for (; i<n; ++i) { if (s[i]=='x') break; }
   if (i<n) { const char *s1=s+i+1; char *s2;
      double x=strtod(s1,&s2);

      if (s2==s1) wblog(FL,"ERR %s() invalid factor in '%s'",C_FCT,s);
      while (isspace(*s2)) { ++s2; } 
      if (*s2) { wblog(FL,"ERR %s() "
         "invalid trailing string\nafter factor in '%s'",C_FCT,s); }
      fac*=x;

      while (--i<n && isspace(s[i])) {}
      n=i+1;
   }

   if (n) {
      for (i=n-1; i<n; --i) {
         if (CTR_IS_CONJ(s[i])) { ++conj; } else 
         if (!isspace(s[i])) { break; }
      }
      if (conj) { conj%=2; } 
      n=i+1;
   }

   while (n>=3 && !strncmp(s+n-3,"^-1",3)) { ++inv;
      for (n-=3; n && isspace(s[n-1]); --n) { }
   }; inv%=2;

   for (i=0; i<n && isdigit(s[i]); ++i) { }

   if (i==n) { WBPERM::init(n); 
      for (i=0; i<n; ++i) { data[i]=s[i]-'0'; }
      if (int(offset)>=0) {
         for (i=0; i<n; ++i) {
            if (data[i]<offset) { if (F) wblog(FL,"ERR %s() "
               "invalid offset ('%s': %d !?)",FCT,s,offset);
               else return -2;
            }
            data[i]-=offset;
         }
      }
      else {
         if (min()==1) { (*this)-=1; }
      }
   }
   else {
      if (isdigit(s[i]) || isspace(s[i]) || strchr("[,;]",s[i])) {
         i=Wb::Str2Idx(0,0,s,*this,offset,0,n);
      }
      else { i=-(i+1); } 

      if (int(i)<0) { if (F) wblog(FL, 
         "ERR %s() invalid permutation \"%s\" @ %d\n"
         "having extended format (offset=%d)",FCT,s,-i,offset);
         else return -3;
      }
   }

   return ( isValidPerm(FL_)>0 ? len : -4 );
};

wbstring wbperm::toStr() const {
   unsigned i=0, nx=(fac==1 ? 0 : 16);
   wbvec<char> sx;

   if (inv ) { nx+=3; } 
   if (conj) { nx+=1; } 

   if (len<10) { 
      sx.init( (len?len:2) + nx+1);
      if (len) { for (; i<len; ++i) { sx.append(0,0,'1'+data[i]); }}
   }
   else {
      if (nx) { nx+=3; } 
      sx.init(nx + len*(log10(double(len))+2));
      if (nx) { sx.append('['); }
      for (; i<len; ++i) { sx.catf(0,0," %d",(int)(data[i]+1)); }
      if (nx)
           { sx.cat(" ]"); }
      else { sx.append(' '); }
   }

   if (inv   ) { sx.cat("^-1");  }
   if (conj  ) { sx.append('*'); } 
   if (fac!=1) { sx.catf(0,0," x%.4g",fac); }
   sx.check_bounds(FL);

   return sx.data;
};

wbperm& wbperm::init_trafo(const wbperm &p1, const wbperm &p2) {

   if (this==&p1 || this==&p2) {
      wbperm X; X.init_trafo(p1,p2);
      return X.save2(*this);
   }

   if (!p1.len || !p2.len) { 
      if (!p1.len)
           { if (p2.len) { init(p2); } else { init(); }}
      else { init(p1,'i'); }

      conj=p1.conj; Wb::conj_add_z2(conj, p2.conj);
      return *this;
   }

   if (p1.len!=p2.len) wblog(FL,
      "ERR %s() length mismatch (%d/%d)",FCT,p1.len,p2.len);

   init(p1,'i').Permute(p2); 
   return *this;
};

wbperm& wbperm::initMove(unsigned k, unsigned l, unsigned N) {
   size_t i=0;        

   if (int(k)<0) { k+=N; }
   if (int(l)<0) { l+=N; }
   if (k>=N || l>=N) wblog(FL,
      "ERR %s() index out of bounds (%d,%d / %d)",FCT,k,l);

   if (k==l) { return Index(N); }

   RENEW(N,NULL,0,0);

   if (k<l) {
      for (; i<k; ++i) { data[i]=i;   }
      for (; i<l; ++i) { data[i]=i+1; }; data[i++]=k;
      for (; i<N; ++i) { data[i]=i;   }
   }
   else if (k>l) { 
      for (; i<l; ++i) { data[i]=i;   }; data[i++]=k;
      for (;i<=k; ++i) { data[i]=i-1; }
      for (; i<N; ++i) { data[i]=i;   }
   }

   return *this;
};

wbperm& wbperm::initFirstTo(wperm_t k, wperm_t N) {
   size_t i=0; RENEW(N,NULL,0,0); 
   if (k>=N) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,N);
   for (   ; i<k; ++i) { data[i]=i+1; }; data[k]=0;
   for (++i; i<N; ++i) { data[i]=i;   };
   return *this;
};

wbperm& wbperm::initLastTo(wperm_t k, wperm_t N) {
   size_t i=0; RENEW(N,NULL,0,0); 
   if (k>=N) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,N);
   for (   ; i<k; ++i) { data[i]=i;   }; data[k]=N-1;
   for (++i; i<N; ++i) { data[i]=i-1; };
   return *this;
};

wbperm& wbperm::init2Front(wperm_t k, wperm_t N) { 
   if (k>=N) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,N);
   wperm_t i=0; RENEW(N,NULL,0,0); data[0]=k;
   for (   ; i<k; ++i) { data[i+1]=i; }
   for (++i; i<N; ++i) { data[i  ]=i; }
   return *this;
};

wbperm& wbperm::init2End(wperm_t k, wperm_t N) { 
   if (k>=N) wblog(FL,"ERR %s() index out of bounds (%d/%d)",FCT,k,N);
   wperm_t i=0; RENEW(N,NULL,0,0);
   for (   ; i<k; ++i) { data[i  ]=i; }
   for (++i; i<N; ++i) { data[i-1]=i; }; data[N-1]=k;
   return *this;
};

wbperm& wbperm::init2Front(const wbvector<wperm_t> &I, wperm_t N) {

   wperm_t i,j;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   RENEW(N,NULL,0,0);
   for (i=0; i<I.len; ++i) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N);
      data[i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (i<N) { 
      for (j=0; j<N; ++j) { if (!m[j]) data[i++]=j; }
   }
   return *this;
};

wbperm& wbperm::init2End(const wbvector<wperm_t> &I, wperm_t N) {

   wperm_t i,j, l=N-I.len;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   RENEW(N,NULL,0,0);
   for (i=0; i<I.len; i++) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N);
      data[l+i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (l) { 
      for (i=j=0; j<N; j++) { if (!m[j]) data[i++]=j; }
   }
   return *this;
};

wbperm& wbperm::init2FrontB(wperm_t m, wperm_t N) { 

   if (!m || m>N) wblog(FL,
      "ERR %s() m out of bounds (%d/%d)",FCT,m,N);

   wperm_t i,l=N-m; init(N);

   for (i=0; i<m; ++i) { data[i]=i+l; }
   for (   ; i<N; ++i) { data[i]=i-m; }

   return *this;
};

wbperm& wbperm::init2EndB(wperm_t m, wperm_t N) { 

   if (!m || m>N) wblog(FL,
      "ERR %s() m out of bounds (%d/%d)",FCT,m,N);

   wperm_t i,l=N-m; init(N);

   for (i=0; i<l; ++i) { data[i]=i+m; }
   for (   ; i<N; ++i) { data[i]=i-l; }

   return *this;
};

wbperm& wbperm::Cycle(wperm_t k1, wperm_t k2) { 
   if (k1>=len || k2>=len) wblog(FL,
      "ERR %s() index out of bounds [%d %d]/%d !?",FCT,k1,k2,len);
   if (k1==k2) return *this;

   wperm_t x1=data[k1], k=k1;
   if (k1<k2)
        { for (; k<k2; ++k) { data[k]=data[k+1]; }}
   else { for (; k>k2; --k) { data[k]=data[k-1]; }}
   data[k]=x1;

   return *this;
};

wbperm& wbperm::Rotate(wperm_ts k) {

    if (!len) wblog(FL,
      "ERR %s() got empty permutation %d/%d !?",FCT,k,len);
    if ((k%=len)<0) k+=len; 
    if (k) {
       wbperm P0(*this);
       for (wperm_t i=0; i<len; ++i) { data[(i+k)%len]=P0[i]; }
    }
    return *this;
};

wbperm& wbperm::initRotate(wperm_t r, wperm_ts k) {

    fac=1;
    if (!r) {
       if (k) wblog(FL,
          "ERR %s() got empty permutation %d/%d !?",FCT,k,r);
       return Index(r);
    }

    if ((k%=r)<0) k+=r; 
    if (!k) { return Index(r); }

    RENEW(r);
    for (wperm_t i=0; i<len; ++i) { data[(i+k)%len]=i; }
    return *this;
};

bool wbperm::sameAs(const wbperm &b, char lflag) const { 

   bool q=0; if (conj!=b.conj) { return q; }
   unsigned i=0;
   char iflag = (inv+b.inv)%2;    

   if (iflag) 
        { if (!lflag || fabs(fac*b.fac-1)>1e-14) { return q; }}
   else { if (fac!=b.fac) { return q; }}

   if (!len  ) { if ( b.len && !b.isIdentityPerm()) { return q; }} else
   if (!b.len) { if (            !isIdentityPerm()) { return q; }} else
   if (iflag) {
      if (len<=b.len) 
           { for (i=0; i<  len; ++i) { if (b.at(data[i])!=i) { return q; }}}
      else { for (i=0; i<b.len; ++i) { if (at(b.data[i])!=i) { return q; }}}
   }
   else {
      unsigned n=MIN(len,b.len);
      for (i=0; i<n; ++i) { if (data[i]!=b.data[i]) { return q; }}
   }

   if (i && len!=b.len) {
      if (len>b.len)
           { for (; i<  len; ++i) { if (  data[i]!=i) return q; }}
      else { for (; i<b.len; ++i) { if (b.data[i]!=i) return q; }}
   }
   return (q=1);
};

bool wbperm::isReversePerm() const {
    if (!len) { return 0; }
    for (wperm_t m=len-1, i=0; i<len; ++i) {
       if (data[i]!=m-i) { return 0; }
    }
    return 1;
};

bool wbperm::isCyclic2F(wperm_t m, wperm_t n,
    char iflag 
  ) const {

    if (wperm_ts(n)<0) { n = (len>m ? (len-m) : 0); }
    if (!len) { return (m || n ? 0 : 1); }
    if (len!=m+n) wblog(FL,
       "ERR %s() length mismatch (%d / %d + %d)",FCT,len,m,n);

    if (iflag) { SWAP(m,n); }

    wperm_t i=0;
    for (; i<m;   ++i) { if (data[i]!=i+n) { return 0; }}
    for (; i<len; ++i) { if (data[i]!=i-m) { return 0; }}
    return 1;
};

int wbperm::getTranspositionsNN(wbMatrix<unsigned> &T2) const {

   unsigned i,j, l=-1;

   wbperm P(*this); wperm_t p_, *p=P.data;

   T2.init((len*(len-1))/2,2);  
   unsigned *t2=T2.data;

   for (i=0; i<len; ++i) { if (p[i]!=i) {
      for (j=i+1; j<len; ++j) { if (p[j]==i) { break; }}
      if (j==len) wblog(FL,
         "ERR %s() invalid permutation [%d-%d: %s]",FCT,i,j,STR(P));

      for (; j>i; --j, t2+=T2.dim2) {
         if (++l>=T2.dim1) wblog(FL,"ERR %s() "
            "nT2 out of bounds (n=%d -> l=%d/%d",FCT,len,l,T2.dim1);
         if (p[j-1]<p[j]) 
              { t2[0]=p[j-1]; t2[1]=p[j]; }
         else { t2[1]=p[j-1]; t2[0]=p[j]; }

         p_=p[j-1]; p[j-1]=p[j]; p[j]=p_;
      }

      if (p[i]!=i) wblog(FL,"ERR %s() "
      "failed to transform to Id [%s] -> [%s]",FCT,STR(*this),STR(P));
   }}

   if (++l==0)
        { T2.init(); } 
   else { T2.dim1=l; }

   return l;
};

bool wbperm::isOpTranspose() const {

    if (len!=3 && len%2) { return 0; }

    for (wperm_t i=0, r=len/2; i<r; ++i) {
        if (data[i]!=i+r || data[i+r]!=i) return 0;
    }
    return 1;
};

wbperm& wbperm::Permute(wbperm P, unsigned r) {

   char q=(isEmpty() ? 1 : (P.isEmpty() ? 2:0));
   if (q) {
      if (q&1) { init(P); }
      if (int(r)>=0 && len!=r) { Extend(r); } 
      return *this;
   }

   P.flatten();

   if (isIdentityPerm()) {
      double f0=((inv%2) ? 1/fac : fac);
      char c0=conj;

      init(P); 
      if (inv) wblog(FL,"ERR %s() got inv=%d",FCT,inv);

      fac*=f0;
      if (c0) { Wb::conj_add_z2(conj,c0); }

      if (int(r)>=0 && len!=r) { Extend(r); }
      return *this; 
   }

   this->flatten();
   Wb::conj_add_z2(conj,P.conj);
   fac *= P.fac;

   if (!P.isIdentityPerm(FL)) {
      unsigned i, n=P.len;
      wperm_t *p=P.data;

      if ((int(r)<0 || r<P.len) && n) { 
         for (i=n-1; i; --i) { if (p[i]!=i) { n=i+1; break; }}
         if (int(r)>=0 && n<r) { n=r; } 
      }

      if (len>=n && (int(r)<0 || r==len)) {
         WBPERM X(n,data); 
         for (i=0; i<n; ++i) { data[i] = X.data[p[i]]; } 
      }
      else { wbperm X(n);
         for (i=0; i<n; ++i) { X.data[i] = at_(p[i]); }
         if (int(r)>=0 && len!=r) { X.Extend(r); } 
         X.WBPERM::save2(*this); 
      }
   }

   return *this; 
};

#endif

