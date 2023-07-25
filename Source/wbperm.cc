/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace permute routines
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

#ifndef __WB_PERM_CC__
#define __WB_PERM_CC__

// ----------------------------------------------------------------- //
// ----------------------------------------------------------------- //

char wbperm::isValidPerm(const char *F, int L, wperm_t r) const {

   char rval=0, id=0; 

   if (swperm_t(r)>=0 && r!=len) { if (F) wblog(F,L,
      "ERR invalid permutation [%s; %d]",toStr().data,r);
      return rval;
   }

   if (len) { wbvector<char> mark(len);
      for (wperm_t i=0; i<len; ++i) {
         if (i==data[i]) { ++id; } else
         if (data[i]>=len || mark[data[i]]++) { if (F) wblog(F,L,
            "ERR invalid permutation [%s]",toStr().data);
            return rval;
         }
      }
   }
   return (rval=( unsigned(id)==len? 2:1));
};

wbperm& wbperm::init(const wbperm &P, char iflag, unsigned r) {
   if (&P!=this) {
      if (iflag)
           { P.invert(*this); }
      else { RENEW(P.len,P.data); }
   }
   else if (iflag) { Invert(); }

   if (int(r)>=0) { Extend(r); } 

   return *this;
};

wbperm& wbperm::invert(wbperm &iP) const { 
   if (&iP==this) {
      wbperm X; this->invert(X).save2(iP);
      return iP;
   }

 #ifndef SKIP_WB_ASSERT
   isValidPerm(FL);
 #endif

   iP.RENEW(len);

   for (wperm_t *p=iP.data, i=0; i<len; ++i) {
      if (data[i]>=len) wblog(FL,"ERR %s() "
         "permutation out of bounds (%d/%d)",FCT,data[i],len);
      p[data[i]]=i;
   }
   return iP;
};

wbperm& wbperm::Extend(wperm_t r) {

   if (r!=len) {
   if (int(r)<=0) wblog(FL,"ERR %s() requesting r=%d/%d",FCT,r,len);

   if (r>len) {
      wbperm P(r);  
      if (len) { memcpy(P.data,data,len*sizeof(wperm_t)); }
      P.save2(*this);
   }
   else {
      for (unsigned i=r; i<len; ++i) {  
         if (data[i]!=i) wblog(FL,
         "ERR %s() invalid r=%d having %s",FCT,r,STR_(this));
      }
      Shorten2(r);
   }}
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
               "invalid partial perm (non-unique: %s)",FCT,STR_(this));
         }
         for (i=0; i<len; ++i) { if (!mark[i]) { data[l++]=i; }}
      }
      else { init(l=len); }
   }
   else if (l>len) wblog(FL,
     "ERR %s() invalid usage (l=%d/%d)",FCT,l,len);

   return *this;
};

int wbperm::initStr(
   const char *F, int L, const char *s, wperm_t offset
){

   if (!s) { init(); return -1; } 
   while (isspace(*s) && *s) ++s; 
   if (!*s){ init(); return 0; }  

   unsigned i=0, n=0;
   while (isdigit(s[n])) ++n;
   for (i=n; isspace(s[i]); ++i) {}; 

   if (!s[i]) { init(n); 
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
         char mark[n+1]; memset(mark,0,n+1);
         for (i=0; i<n; ++i) { 
            if (data[i]>n || ++mark[data[i]]!=1) { return -4; }
         }
         if (!mark[0]) {
            for (i=0; i<n; ++i) { --data[i]; }
         }
         return len; 
      }
   }
   else {
      if (isdigit(s[i]) || isspace(s[i]) || s[i]==',' || s[i]==';') {
         i=Str2Idx(s,*this,offset);
      }
      else { i=-(i+1); } 

      if (int(i)<0) { if (F) wblog(FL, 
         "ERR %s() invalid permutation (i=%d, offset=%d)"
         "%N%N    s = '%s' %N",FCT,-i,offset,s);
         else return -3;
      }
   }

   i=isValidPerm(F ? F : NULL,F ? L : 0);

   return (i ? len : -4);
};

wbperm& wbperm::init_trafo(const wbperm &p1, const wbperm &p2) {

   if (!p1.len) {
      if (p2.len) return init(p2); 
      else return init();
   }
   if (!p2.len) { return init(p1,'i'); }  

   if (p1.len!=p2.len) wblog(FL,
      "ERR %s() length mismatch (%d/%d)",FCT,p1.len,p2.len);

   if (this==&p1 || this==&p2) {
      wbperm X; X.init_trafo(p1,p2);
      X.save2(*this); return *this;
   }

   init(p1,'i').Permute(p2); 

   return *this;
};

wbperm& wbperm::init2Front(const wbvector<wperm_t> &I, wperm_t N) {

   wperm_t i,j;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   init(N);
   for (i=0; i<I.len; i++) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N);
      data[i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (i<N) 
   for (j=0; j<N; j++) { if (!m[j]) data[i++]=j; }

   return *this;
};

wbperm& wbperm::init2End(const wbvector<wperm_t> &I, wperm_t N) {
   wperm_t i,j, l=N-I.len;
   wbvector<char> mark(N); char *m=mark.data;
   if (I.len>N) wblog(FL,"ERR index too long (%d/%d)",I.len,N);

   init(N);
   for (i=0; i<I.len; i++) { j=I[i];
      if (j>=N) wblog(FL,"ERR index out of bounds (%d/%d)",j,N);
      data[l+i]=j; if ((++m[j])>1) {
      wblog(FL,"ERR %s() index not unique (%d/%d)",FCT,j,N); }
   }

   if (l) 
   for (i=j=0; j<N; j++) { if (!m[j]) data[i++]=j; }

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

wbperm& wbperm::init2Front(wperm_t k, wperm_t N) { 
   if (k>=N) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,N);
   wperm_t i=0; RENEW(N); data[0]=k;
   for (; i<k; ++i) { data[i+1]=i; }; ++i;
   for (; i<N; ++i) { data[i  ]=i; };
   return *this;
};

wbperm& wbperm::init2End(wperm_t k, wperm_t N) { 
   if (k>=N) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,k+1,N);
   wperm_t i=0; RENEW(N); data[N-1]=k;
   for (; i<k; ++i) { data[i  ]=i; }; ++i;
   for (; i<N; ++i) { data[i-1]=i; };
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

wbperm& wbperm::initCycle(wperm_t r, wperm_t k1, wperm_t k2) { 

   if (k1>=r || k2>=r) wblog(FL,
      "ERR %s() index out of bounds [%d %d]/%d !?",FCT,k1,k2,len);
   if (k1==k2) { return Index(r); }

   RENEW(r);

   if (k1<k2) { wperm_t k=0;
      for (; k<k1; ++k) { data[k]=k; }
      for (; k<k2; ++k) { data[k]=k+1; }; data[k++]=k1;
      for (; k<r;  ++k) { data[k]=k; }
   }
   else { wperm_t k=r-1;
      for (; k>k1; --k) { data[k]=k; }
      for (; k>k2; --k) { data[k]=k-1; }; data[k--]=k1;
      for (; k<r;  --k) { data[k]=k; }
   }

   return *this;
};

wbperm& wbperm::Rotate(swperm_t k) {

    if (!len) wblog(FL,
      "ERR %s() got empty permutation %d/%d !?",FCT,k,len);
    if ((k%=len)<0) k+=len; 
    if (k) {
       wbperm P0(*this);
       for (wperm_t i=0; i<len; ++i) { data[(i+k)%len]=P0[i]; }
    }
    return *this;
};

wbperm& wbperm::initRotate(wperm_t r, swperm_t k) {

    if (!r) {
       if (k) wblog(FL,
          "ERR %s() got empty permutation %d/%d !?",FCT,k,r);
       return Index(r);
    }

    if ((k%=r)<0) k+=r; 
    if (!k) return Index(r);

    RENEW(r);
    for (wperm_t i=0; i<len; ++i) { data[(i+k)%len]=i; }
    return *this;
};

bool wbperm::sameAs(const wbperm &b) const { 
   if (!len) {
      if (b.len)
           { return b.isIdentityPerm(); }
      else { return 1; }
   }
   if (!b.len)
        { return isIdentityPerm(); }
   else { return (*this)==b; }
};

bool wbperm::isIdentityPerm() const {
    bool id=1; if (len) { 
       for (unsigned i=0; i<len; ++i) { if (data[i]!=i) { id=0; break; }};
    }
    return id;
};

bool wbperm::isIdentityPerm(const char *F, int L, wperm_t r) const {
    if (int(r)>=0 && len!=r) { if (F) wblog(F,L,
       "ERR %s() invalid permutation (len=%d/%d)",FCT,len,r);
       return 0;
    }
    bool q=1; char mark[len]; memset(mark,0,len);
    wperm_t i=0;

    for (; i<len; ++i) {
       if (data[i]==i) { ++mark[i]; } else { q=0; break; }
    }
    for (; i<len; ++i) {
       if (data[i]>=len) wblog(F_L, "ERR %s() invalid permutation\n"
          "(index out of range %d/%d)",FCT,data[i],len);
       if ((++mark[data[i]])>1) wblog(F_L,"ERR %s() "
          "invalid permutation (non-unique index)",FCT);
    }
    return q;
};

bool wbperm::isReversePerm() const {
    if (len==0) return 0;
    for (wperm_t m=len-1, i=0; i<len; i++) if (data[i]!=m-i) return 0;
    return 1;
};

bool wbperm::isCyclic2F(wperm_t m, wperm_t n,
    char iflag 
  ) const {

    if (swperm_t(n)<0) { n = (len>m ? (len-m) : 0); }
    if (!len) { return (m || n ? 0 : 1); }
    if (len!=m+n) wblog(FL,
       "ERR %s() length mismatch (%d / %d + %d)",FCT,len,m,n);

    if (iflag) { SWAP(m,n); }

    wperm_t i=0;
    for (; i<m;   ++i) { if (data[i]!=i+n) { return 0; }}
    for (; i<len; ++i) { if (data[i]!=i-m) { return 0; }}
    return 1;
};

bool wbperm::isOpTranspose() const {

    if (len!=3 && len%2) { return 0; }

    for (wperm_t i=0, r=len/2; i<r; ++i) {
        if (data[i]!=i+r || data[i+r]!=i) return 0;
    }
    return 1;
};

wbperm& wbperm::Permute(const wbperm &P, char iflag, wperm_t r){

   if (P.len && !P.isIdentityPerm(FL,r)) { 
      if (!len || isIdentityPerm(FL,r)) {
         if (iflag)
              { P.invert(*this); }
         else { init(P); }
      }
      else {
         wbperm X(*this); 
         wperm_t i=0, *x=X.data, *p=P.data;

         if (iflag)
              { for (; i<r; ++i) { data[p[i]] = x[i]; }}
         else { for (; i<r; ++i) { data[i] = x[p[i]]; }}
      }
   }
   return *this;
};

#endif

