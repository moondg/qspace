/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbvector (template vector class)
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

#ifndef __WB_PERM_HH__
#define __WB_PERM_HH__

// ----------------------------------------------------------------- //
// permutation class
// Wb,Apr19,06
// ----------------------------------------------------------------- //

class wbperm : public wbvector<wperm_t> { 

  public:

    wbperm() : inv(0), conj(0), fac(1) {};

    wbperm(wperm_t n, char reverse=0)
    : inv(0), conj(0), fac(1) { Index(n,reverse); };

    wbperm(const wbperm &P, char iflag=0, unsigned r=-1)
     : inv(0), conj(0), fac(P.fac) { init(P,iflag,r); };

    wbperm(const char *F, int L,
       const mxArray *a, unsigned offset=-1, unsigned r=-1)
     : inv(0), conj(0), fac(1) { init(F,L,a,offset,r); };

    template<class T>
    explicit wbperm(const wbvector<T> &b, wperm_t offset=0, double fac_=1)
     : inv(0), conj(0), fac(1) { init(FL,b,offset,fac_); };

    explicit wbperm(const wbindex &P, double fac_=1); 

    wbperm(const char *s, wperm_t offset=1)
     : inv(0), conj(0), fac(1) { initStr(FL,s,offset); };

    wbperm& init_bare(unsigned n) { WBPERM::init(n); return *this; };

    wbperm& init(unsigned n=0, char reverse=0) {
       inv=conj=0; fac=1;
       return Index(n,reverse);
    };

    wbperm& init(unsigned n, const wperm_t* d, char ref=0, double fac_=1) {
       WBPERM::init(n,d,ref); inv=conj=0; fac=fac_; isValidPerm(FL);
       return *this;
    };

    wbperm& init(const wbperm &P, char iflag=0, unsigned r=-1);

    wbperm& init_trafo(const wbperm &p1, const wbperm &p2);

    int init(const char *s, wperm_t offset=1) { 
       return initStr(FL,s,offset);
    };

    int initStr( 
       const char *F, int L, const char *s, wperm_t offset=1);

    wperm_t at_(wperm_t i) const { return (i<len ? data[i] : i); };

    wbperm& operator= (const wbperm &P) {
       WBPERM::init(P); inv=P.inv; conj=P.conj; fac=P.fac;
       return *this;
    };

    wbperm& operator= (const char *s) {
       initStr(FL,s,1); 
       return *this;
    };

    void swap(wbperm &P) { WBPERM::swap(P);
       SWAP(inv ,P.inv );
       SWAP(conj,P.conj);
       SWAP(fac ,P.fac );
    };

    wbperm& save2(wbperm &P) {
       WBPERM::save2(P); P.inv=inv; P.conj=conj; P.fac=fac;
       inv=conj=0; fac=1;
       return P;
    };

    wbperm& Extend(wperm_t r);

    wbperm& init(const char *F, int L, 
       const mxArray* a, unsigned offset=-1, unsigned r=-1) {

       inv=conj=0; fac=1;

       if (mxIsNumeric(a)) {
          WBPERM::init(F,L,a); if (len) {
             if (int(offset)<0) { offset=min(); }
             if (offset) {
                if (offset!=1) wblog(FL,
                   "WRN %s() using offset=%d",C_FCT,offset);
                operator-=(offset);
             }
             isValidPerm(F_L);
          }
       }
       else if (mxIsChar(a)) {
          wbvec<char> s(mxGetNumberOfElements(a)+1);
          mxGetString(a,s.data,s.len);
          if (int(offset)<0) { offset=1; } 
          initStr(F_L,s.data,offset); 
       }
       else wblog(FL,
          "ERR %s() invalid input (%s) !?",FCT,mxGetClassName(a));

       if (int(r)>=0) { Extend(r); }

       return *this;
    };

    template<class T>
    wbperm& init(const char *F, int L,
       const wbvector<T> &b, wperm_t offset=0, double fac_=1){
       initT(b.len,b.data); inv=conj=0; fac=fac_;
       if (offset) { (*this)-=offset; }; isValidPerm(F_L);
       return *this;
    };

    wbperm& initTranspose(wperm_t r) {
       RENEW(r,NULL,0,0); wperm_t r2=r/2;
       if (r%2) wblog(FL,"ERR %s() even input rank required (%d)",FCT,r);
       for (wperm_t i=0; i<r2; ++i) { data[i]=r2+i; data[r2+i]=i; }
       inv=conj=0; fac=1;
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

    wbperm& initMove(unsigned k, unsigned l, unsigned N); 

    wbperm& initFirstTo(wperm_t k, wperm_t N); 
    wbperm& initLastTo (wperm_t k, wperm_t N); 

    wbperm& init2Front(const WBPERM &I, wperm_t N);
    wbperm& init2End  (const WBPERM &I, wperm_t N);

    wbperm& init2FrontB(wperm_t m, wperm_t N); 
    wbperm& init2EndB  (wperm_t m, wperm_t N);

    wbperm& init2Front(wperm_t k, wperm_t N);  
    wbperm& init2End  (wperm_t k, wperm_t N);  

    wbperm& Cycle(wperm_t k1, wperm_t k2); 

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

    inline wbperm& Conj(char cflag=1) {
       Wb::conj_add_z2(conj,cflag);
       return *this;
    };

    inline wbperm& flatten(unsigned r=-1) {
       if ((inv%=2)) {
          if (len) { invert_data(); } 
          fac=1./fac; 
          inv=0; 
       }
       if (int(r)>=0) { Extend(r); } 
       return *this;
    };

    wbperm& compact() { 
       if (len) { unsigned i=len-1;
          for (; i<len; --i) { if (data[i]!=i) break; }
          if (++i<len) {
             if (i) { len=i; }
             else { WBPERM::init(); inv=0; }
          }
       }
       return flatten();
    };

    inline wbperm& adapt(char iflag, char cflag=0, double fac_=1) {
       if (iflag) { ++inv; } 
       if (cflag) { Wb::conj_add_z2(conj,cflag); }

       flatten(); fac*=fac_;
       return *this;
    };

    inline char isValidPerm(wperm_t r=-1, char f=1) const;

    inline char isValidPerm(
       const char *F, int L, wperm_t r=-1, char f=1) const {

       char q=isValidPerm(r,f);
       if (q<=0 && F) { const unsigned n=32; char sx[n];
          snprintf(sx,n,"(len=%ld/%ld, q=%d)",len,r,q);
          if (len<10)
               { wblog(FL,"ERR %s() '%s' %s",FCT,STR(*this),sx); }
          else { wblog(FL,"ERR %s() %s",FCT,sx); }
       }
       return q;
    };

    wbperm& Complete (wperm_t l);

    wbperm& Complete1();

    wbperm& invert(wbperm &iP) const;
    wbperm& Invert(char cflag=0);
    wbperm  inverse() const { wbperm iP; return invert(iP); }

    wbperm& Rotate(wperm_ts l); 
    wbperm& initRotate(wperm_t n, wperm_ts k);

    wbperm& Permute(wbperm P, unsigned r); 
    wbperm& Permute(const wbperm &P) {
       return Permute(P, MAX(len,P.len));
    };

    bool sameAs(const wbperm &b, char lflag=1) const; 

    bool operator==(const wbperm &b) const { return  sameAs(b,0); };
    bool operator!=(const wbperm &b) const { return !sameAs(b,0); };

    bool isIdentityPerm(const char *F, int L, wperm_t r=-1) const;

    bool isIdentityPerm() const { 
       for (wperm_t i=0; i<len; ++i) { if (data[i]!=i) { return 0; }}
       return 1;
    };

    bool isReversePerm() const;
    bool isCyclic2F(wperm_t m, wperm_t n=-1, char iflag=0) const;

    bool isOpTranspose() const;

    int getTranspositionsNN(wbMatrix<unsigned> &T2) const;

    explicit operator bool() const {
       return (!isIdentityPerm() || conj || fac!=1 ? 1 : 0); };

    bool operator! () const {
       return (!isIdentityPerm() || conj || fac!=1 ? 0 : 1); };

    wbperm& Strip() { 
       fac=1; conj=0; return flatten();
    };

    inline wbperm operator() (char ref=0) const {
       wbperm B; 
       B.WBPERM::init(len,data,ref); B.inv=inv;
       B.conj=0; B.fac=1; 
       return B;
    };

    bool isEmpty() const {
       return (len || (conj%2) || fac!=1 ? 0 : 1);
    };

    char extras(char qref=0) const { char q=0; 
       if (inv %2) { q|=1; }
       if (conj%2) { q|=2; }
       if (fac!=1) { q|=4; }
       if (abs(inv)>1 || conj<0 || conj>1) { q|=64; } 
       return (q^=(q&qref)); 
    };

    char relevant(const char *F=NULL, int L=0) const { 
       char q=0;             
          if (!isIdentityPerm(F_L)) { q|=1; } 
          if (fac!=1) { q|=2; }
          if (conj) { q|=4; }
       return q;
    };

    char relevant(const char *F, int L, char qm) const {
       char q=relevant(F_L);
       if (qm>=0) {
        if ( q!=(q&(qm|1)) || (!(qm&1) && inv) ) { 
           if (F) wblog(F,L,
              "ERR %s() unexpected %s (qm=%d)",FCT,STR(*this),qm);
           q=-q;
       }}
       return q;
    };

    wbperm& times(const wbperm &p2, wbperm &Pout) const {
       if (   isValidPerm()<=0) wblog(FL,
          "ERR %s() invalid permutation (this: %s)",__FUNCTION__,STR(*this));
       if (p2.isValidPerm()<=0) wblog(FL,
          "ERR %s() invalid permutation (p2: %s)",__FUNCTION__,STR(p2));
       if (len!=p2.len) wblog(FL,
          "ERR %s() incompatible permutations (len=%d/%d)",
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

    wbstring toStr() const;

    void wberr_invalid_perm(const char *F, int L, unsigned l=-1) const {
       wbvec<char> s(64);
          s.catf(0,0,"%s() invalid perm ",PROG); if (len<16) {
          s.catf(0,0," %s",STR(*this)); };
          if (int(l)>=0) { s.catf(0,0," (%d/%d)",len,l); } else
          if (    len>4) { s.catf(0,0," (len=%d)", len); }
       wblog(FL,"ERR %s",s.data);
    };

    mxArray* toMx(const char tflag=0) const {
       if (fac==1 && !conj) {
          if (inv)
               { wbperm iP; invert(iP); return iP.toMx_offset(1,tflag); }
          else { return toMx_offset(1,tflag); } 
       }
       else { return toStr().toMx(); }
    };

    char inv;   

    char conj;  
    double fac; 

  protected:
  private:

   void invert_data(wperm_t *p) const {
      wbvector<char> mark(len); char *m=mark.data;

      for (unsigned i=0; i<len; ++i) {
         if (data[i]<len && ++m[data[i]]==1) { p[data[i]]=i; }
         else wblog(FL,"ERR %s() invalid permutation %s",FCT,STR(*this));
      }
   };

   wbperm& invert_data() { 
      wbvector<char> mark(len); char *m=mark.data;
      WBPERM P(len,data); wperm_t *p=P.data;

      for (unsigned i=0; i<len; ++i) {
         if (p[i]<len && ++m[p[i]]==1) { data[p[i]]=i; }
         else wblog(FL,"ERR %s() invalid permutation [%s]",FCT,STR(P));
      }
      return *this;
   };

};

namespace Wb {

template <class T, ENABLE_IF_isINT(T)>
char is_valid_perm( 
    const T *data, size_t len, size_t n, char f=1) {

    char q=(len ? 1 : 2); 
    size_t i;

    if (f>3) {
       if (f=='l') { f=1; } else
       if (f=='L') { f=3; } 
    }
    if (f<0 || f>3) wblog(FL,"ERR %s() invalid f=%s",FCT,cSTR(f));

    if (long(n)<0) { n=len; } else
    if (len<n) { q|=4; if (!(f&1)) { return -q; } else n=len; } else
    if (len>n) { q|=8; if (!(f&2)) { return -q; } else {
       for (i=n; i<len; ++i) { if (data[i]!=i) { return -q; }}
    }}

    if (n) { 
       if (n==1) { if (data[0]) { return -q; } else { q|=2; }}
       else {
          size_t id=0; 
          wbvector<char> mark(n);

          for (i=0; i<n; ++i) {
             if (data[i]<0 || data[i]>=n || mark[data[i]]++) { return -q; }
             if (data[i]==i) { ++id; }
          }
          if (id==n) { q|=2; } 
       }
    }
    return q;
};

char isValidPerm( 
   const char *s, char f=1, wperm_t *r_=NULL) {

   wbperm P; int l=P.initStr(0,0,s);
   wperm_t r=-1; if (r_) { r=*r_; *r_=P.len; }
   if (l>=0) 
        { return is_valid_perm(P.data,P.len,r,f); }
   else { return -32; } 
};

char isValidPerm( 
  const mxArray *a, char f=1, wperm_t *r_=NULL) {

   char q=0; 
   int i;
   wperm_t r=(r_ ? *r_ : -1);

   if (!a) { q=-11; } else
   if (mxIsChar(a)) { wbstring S(FL,a);
      if (S) { q=isValidPerm(S.data,f,&r); } 
   }
   else if (!mxIsNumeric(a)) { q=-12; }
   else if ((i=Mx::IsVector(a))<=0) {
      q=(i && !(r=mxGetNumberOfElements(a)) ? 2 : -13);
   }
   else {
      wbvector<double> P(FL,a);
      if (P) {
         if (P.min()==1) { P-=1; }
         q=is_valid_perm(P.data,P.len,r,f);
         r=P.len;
      }
   }
   if (r_) { *r_=r; }
   return q;
};

template <class T, class TI, ENABLE_IF_isINT(TI)>
inline void perm_data(const char *F, int L, 
   T* b, const T *a, size_t N, const TI *p, size_t n, char iflag=0) {

   size_t i;

   if (!a || !b || !p || a==b) wblog(F_L,
      "ERR %s() invalid input (%x / %x / %x)",FCT,a,b,p);
   if (n>N) {
      for (i=N; i<n; ++i) { if (p[i]!=i) wblog(F_L,
         "ERR %s() permutation out of bounds (len=%d/%d) ",FCT,n,N); }
      n=N;
   }

   if (!iflag)
		{ for (i=0; i<n; ++i) { b[i]=a[p[i]]; }}
   else { for (i=0; i<n; ++i) { b[p[i]]=a[i]; }}

   for (; i<N; ++i) { b[i]=a[i]; }
};

} 

char wbperm::isValidPerm(wperm_t r, char f) const {
   if (!isfinite(fac)) { return -32; } 
   return Wb::is_valid_perm(data,len,r,f);
};

bool wbperm::isIdentityPerm(const char *F, int L, wperm_t r) const {

   char q=Wb::is_valid_perm(data,len,r,1); 

   if (q<=0) { if (F) wblog(FL,
      "ERR %s() invalid permutation (len=%d/%d)",FCT,len,r);
      return 0;
   }

   return (q&2 ? 1 : 0); 
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

