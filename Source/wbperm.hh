/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbvector (template vector class)
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

#ifndef __WB_PERM_HH__
#define __WB_PERM_HH__

// ----------------------------------------------------------------- //
// permutation class
// Wb,Apr19,06
// ----------------------------------------------------------------- //

class wbperm : public wbvector<wperm_t> { 

  public:

    wbperm() {}; 

    wbperm(wperm_t n, char reverse=0) { Index(n,reverse); };

    wbperm(const wbperm &P, char iflag=0, unsigned r=-1)
     : WBPERM() { init(P,iflag,r); };

    wbperm(const wbperm &P, size_t r) { init_pad(P,r); };

    template<class T>
    wbperm(const wbvector<T> &b, wperm_t offset=0) {
       init(FL,b,offset);
    };

    wbperm(const wbindex &P); 

    wbperm(const char *s, wperm_t offset=1) {
       initStr(FL,s,offset);
    }

    wperm_t el1(wperm_t i) const { 
       if (len) {
          if (i>=len) wblog(FL,
             "ERR element index out of bounds (%d/%d)",i+1/len);
          return data[i];
       }
       else return i;
    };

    wbperm& operator= (const char *s) {
       initStr(FL,s,0); return *this;
    };

    int initStr( 
       const char *F, int L, const char *s, wperm_t offset=1);

    wbperm& init(const wbperm &P, char iflag=0, unsigned r=-1);
    wbperm& Extend(wperm_t r); 

    wbperm& init_trafo(const wbperm &p1_, const wbperm &p2_);

    wbperm& init_pad(const wbperm &P, wperm_t r) {
       if (r<P.len) wblog(FL,
          "ERR invalid %s contruction (%d/%d)",FCT,r,P.len);
       RENEW(r); memcpy(data,P.data,P.len*sizeof(P.data[0]));
       for (wperm_t i=P.len; i<r; i++) data[i]=i;
       if (!isValidPerm()) dispInvalidPerm(FL);
       return *this;
    };

    wbperm& init(unsigned n=0, char reverse=0) {
       return Index(n,reverse);
    };

    wbperm& init(unsigned n, const wperm_t* d, char ref=0) {
       wbvector<wperm_t>::init(n,d,ref); isValidPerm(FL);
       return *this;
    };

    wbperm& init(const char *F, int L, const mxArray* a){
       return init(F,L,a,mxIsChar(a) ? 1 : 0); 
    };

    wbperm& init(const char *F, int L, 
       const mxArray* a, wperm_t offset
    ){
       if (mxIsNumeric(a)) {
          WBPERM::init(F,L,a); if (offset) { operator-=(offset); }
          isValidPerm(F_L);
       }
       else if (mxIsChar(a)) {
          int l=mxGetNumberOfElements(a)+1;
          char s[l]; mxGetString(a,s,l);
          initStr(F_L,s,offset);
       }
       else wblog(FL,
          "ERR %s() invalid input (%s) !?",FCT,mxGetClassName(a));
       return *this;
    };

    template<class T>
    wbperm& init(
       const char *F, int L, const wbvector<T> &b, wperm_t offset=0
    ){
       initT(b.len,b.data);
       if (offset) { (*this)-=offset; }; isValidPerm(F_L);
       return *this;
    };

    wbperm& initTranspose(wperm_t r) {
       WBPERM::init(r); wperm_t r2=r/2;
       if (r%2) wblog(FL,"ERR %s() even input rank required (%d)",FCT,r);
       for (wperm_t i=0; i<r2; ++i) { data[i]=r2+i; data[r2+i]=i; }
       return *this;
    };

    wbperm& initPerm2(wperm_t r) {
       WBPERM::init(r); wperm_t r2=r/2;
       if (r%2) wblog(FL,"ERR %s() even input rank required (%d)",FCT,r);
       for (wperm_t l=0,i=0; i<r2; ++i) {
          data[l++]=i;
          data[l++]=r2+i;
       }
       return *this;
    };

    wbperm& initFirstTo (wperm_t k, wperm_t N) { 
       Index(N); mvFirstTo(k); 
       return *this;
    };

    wbperm& initLastTo (wperm_t k, wperm_t N) {
       Index(N); mvLastTo(k); 
       return *this;
    };

    wbperm& init2Front(const WBPERM &I, wperm_t N);
    wbperm& init2End  (const WBPERM &I, wperm_t N);

    wbperm& init2FrontB(wperm_t m, wperm_t N); 
    wbperm& init2EndB  (wperm_t m, wperm_t N);

    wbperm& init2Front(wperm_t k, wperm_t N);  
    wbperm& init2End  (wperm_t k, wperm_t N);  

    wbperm& Cycle(wperm_t k1, wperm_t k2); 
    wbperm& initCycle(wperm_t n, wperm_t k1, wperm_t k2);

    wbperm& Index(wperm_t n, char reverse=0) {
       RENEW(n); if (n==0) return *this;
       if (reverse)
            for (wperm_t l=n-1, i=0; i<n; ++i) data[i]=l-i;
       else for (wperm_t i=0; i<n; ++i) data[i]=i;
       return *this;
    };

    wbperm& Index(wperm_t i1, wperm_t i2) {
       if (i1>i2) { RENEW(0); return *this; }

       RENEW(i2-i1+1);
       for (wperm_t i=0; i<len; ++i) data[i]=i+i1;

       return *this;
    };

    char isValidPerm(const char *F, int L, wperm_t r=-1) const;
    char isValidPerm(wperm_t r=-1) const { return isValidPerm(0,0,r); };

    wbperm& Complete(wperm_t l);

    wbperm& Invert() { 
       wbperm iP; this->invert(iP).save2(*this);
       return *this;
    };

    wbperm& invert(wbperm &iP) const;

    wbperm& Rotate(swperm_t l); 
    wbperm& initRotate(wperm_t n, swperm_t k);

    wbperm& Permute(const wbperm &P, char iflag, wperm_t r);
    wbperm& Permute(const wbperm &P, char iflag=0) {
       return Permute(P, iflag, len>P.len ? len : P.len);
    };

    bool sameAs(const wbperm &b) const; 
    bool isIdentityPerm() const;
    bool isIdentityPerm(const char *F, int L, wperm_t r=-1) const;

    bool isReversePerm() const;
    bool isCyclic2F(wperm_t m, wperm_t n=-1,  char iflag=0) const;

    bool isOpTranspose() const;

    explicit operator bool() const {
         return (isEmpty() || isIdentityPerm() ? 0 : 1); };
    bool operator! () const { return !bool(*this); }

    wbperm& times(const wbperm &p2, wbperm &Pout) const {
       if (!isValidPerm()) wblog(FL,
          "ERR %s() invalid permutation (%s)",__FUNCTION__,toStr().data);
       if (!p2.isValidPerm()) wblog(FL,
          "ERR %s() invalid permutation (%s)",__FUNCTION__,STR(p2));
       if (len!=p2.len) wblog(FL,
          "ERR %s() incompatible permutations (%d/%d)",
         __FUNCTION__,len,p2.len);

       wbperm P(len); 
       for (wperm_t i=0; i<len; i++) P.data[i]=data[p2.data[i]];

       P.save2(Pout); return Pout;
    };

    wbperm& flipIdx(wbperm &P) const;
    wbperm& FlipIdx();
    wbperm  flipIdx() const {
       wbperm P; flipIdx(P); 
       return P;
    };

    void ToRaw() { ToRaw(len); };
    void ToRaw(wperm_t n) {
        if (n<len) wblog(FL,"WRN ToRaw() with n=%d/%d !?",n,len);
        if (n==0) return; else n--; 

        Flip();

        for (wperm_t i=0; i<len; i++) {
            if (i>n) wblog(FL,"ERR wbperm - invalid n=%d (%d)",n,i);
            data[i]=n-data[i];
        }
    }

    wbstring toStr() const { 

       wbstring s_; 
       unsigned i=0; char *s;

       if (len<10) { 
          s_.init(len+1); s=s_.data;
          for (; i<len; ++i) { s[i]='1'+data[i]; }
          s[i]=0;
       }
       else {
          unsigned l=0, n=len*(log10(double(len))+2);
          s_.init(n); s=s_.data;
          for (; i<len && l<n; ++i) {
             l+=snprintf(s+l,n-l," %d",(int)(data[i]+1));
          }
          if (l+3>=n) wblog(FL,"WRN %s() "
             "string out of bounds (%d/%d; %d/%d)",FCT,l,n,i+1,len);
          else { s[l]=' '; s[l+1]=0; }
       }
       return s_;
    };

    mxArray* toMx(const char tflag=0) const {
       return toMx_offset(1,tflag); 
    };

  protected:
  private:

    void dispInvalidPerm(const char* F, int L) const {
    wblog(F,L,"WRN invalid permutation [%s]",toStr().data); };

};

int isValidPerm(const char *s, wperm_t r=-1) {
   int q=0; wbperm P; int i;
   if ((i=P.initStr(0,0,s))<0) { q=-1; } else
   if (P.isValidPerm(r)) { q=P.len; } else { q=-2; }
   return q;
};

int isValidPerm(const mxArray *a, wperm_t r=-1) {
   int q=0; 

   if (!a) { q=-11; } else
   if (mxIsChar(a)) {
      wbstring S(FL,a); if (S) {
         q=isValidPerm(S.data,r); 
      } 
   }
   else if (!mxIsNumeric(a)) { q=-12; }
   else if (!Mx::IsVector(a)) {
      if (mxGetNumberOfElements(a)) q=-13; 
   }
   else {
      wbvector<double> X(FL,a);
      if (int(r)>=0 && X.len!=r) { q=-3; } else
      if (X) {
         unsigned i=0, j, n=X.len; double *x=X.data;
         wbvector<char> m(n+1);
         for (; i<n; ++i) { j=x[i];
            if (x[i]!=j || j>n || ++m[j]!=1) { break; }
         }
         q=(i<n || (m[0] && m[n]) ? -4 : n);
      }
   }
   return q;
};

wbperm& wbperm::flipIdx(wbperm &P) const {

    if (&P==this) wblog(FL,
       "ERR flipIdx(P) calls itself\nused flipIdx() instead!");

    P.RENEW(len); if (!len) return P;

    for (wperm_t m=len-1, i=0; i<len; i++)
    P.data[i]=m-data[i];

    return P;
};

inline wbperm& wbperm::FlipIdx() {

    for (wperm_t m=len-1, i=0; i<len; ++i) {
       if (data[i]>m) wblog(FL,"ERR index out of bounds (%d/%d)",data[i],m);
       data[i]=m-data[i];
    }

    return *this;
};

#endif

