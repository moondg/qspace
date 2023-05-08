/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace index routines
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

#ifndef __WB_INDEX_CC__
#define __WB_INDEX_CC__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */
// if toend==0 : add towards the end, i.e. P=[I,...]  (default)
// else P=[...,I].

int wbindex::extend2Perm(widx_t N, wbperm &P, char toend) const {

    widx_t i,k;
    wbvector<char> mark(N);

    for (i=0; i<len; ++i) {
       if (data[i]>=N) { sprintf(str,
          "%s:%d %s() index out of bounds (%ld/%ld)",FLF,data[i],N);
          return 1;
       }
       if (mark[data[i]]++) { sprintf(str,
          "%s:%d %s() index not unique",FLF);
          return 1;
       }
    }

    P.init(N);

    if (!toend) { 
       memcpy(P.data,data,len*sizeof(widx_t));

       for (k=len, i=0; i<N; ++i) 
       if (!mark[i]) P[k++]=i;
    }
    else {
       memcpy(P.data+N-len,data,len*sizeof(widx_t));

       for (k=i=0; i<N; ++i)
       if (!mark[i]) P[k++]=i; 
    }

    return 0;
};

wbindex& wbindex::invert(widx_t n, wbindex &I, char uflag) const {

    wbvector<char> M(n); char *m=M.data;
    widx_t i,l;

    for (i=0; i<len; i++) { l=data[i];
        if (l<n) m[l]=1; else wblog(FL,
        "ERR %s() index out of bounds (%d: %d/%d)",FCT,i,l,n);
    }

    l=M.nnz();
    if (uflag && l!=len) wblog(FL,
       "ERR %s() index not unique (%d/%d)",FCT,l,len);

    I.init(n-l);
    if (I.len) {
       for (l=i=0; i<n; ++i) if (!m[i]) I.data[l++]=i;
    }

    return I;
};

wbindex wbindex::flipIdx(widx_t ndim) const {
    wbindex a(len);

    if (ndim==0) {
       if (len==0) return *this; else { dbstop(FL);
       wblog(FL,"ERR wbindex - invalid ndim=%d (%d)", ndim, len); }
    }

    for (widx_t m=ndim-1, i=0; i<len; i++) {
        if (i>m) wblog(FL,"ERR wbindex - invalid ndim=%d (%d)",ndim,i);
        a.data[i]=m-data[i];
    }
    return a;
};

inline char wbindex::isUnique(widx_t imax, widx_t offset) const {
    wbvector<char> mark(imax);

    if (offset==0) {
       for (widx_t i=0; i<len; i++) {
          if (data[i]<imax) { if ((++mark[data[i]])>1) return 0; }
          else return 0;
       }
    }
    else {
       widx_t i=0, i0=offset; imax+=offset;
       for (; i<len; i++) {
          if (data[i]>=i0 && data[i]<imax) {
             if ((++mark[data[i]-i0])>1) return 0;
          }
          else return 0;
       }
    }

    return 1;
};

inline wbindex& wbindex::init(const ctrIdx &ic) {
   initT(ic.len,ic.data); 
   return *this;
};

inline wbindex& wbindex::init( 
   const char *F, int L, const mxArray *a, widx_t offset
){
   if (!a) { init(); return *this; }

   const size_t *sz=mxGetDimensions(a);

   Mx::Array<double> A(a);
   wbvector<double> x;

   if (A.rank!=2 || (sz[0]!=1 && sz[1]!=1)) wblog(F_L,
      "ERR %s() numeric vector required\n%s",FCT,STR(A));

   init(A.len); if (!A.len) { return *this; }
   if (A.data)
        { x.init(A.len,A.data,'r'); }
   else { x.init(A.len); A.rcopy_to(x.data,0,'!'); }

   for (widx_t i=0; i<A.len; ++i) { data[i]=widx_t(x.data[i]);

      if (double(data[i])!=x.data[i]) wblog(F_L,
         "ERR %s() type mismatch (%d: %d %g)",i+1,data[i],x.data[i]);

      if (offset) {
         if (data[i]>=offset) { data[i]-=offset; }
         else if (offset==1) wblog(F_L,"ERR %s() "
            "1-based index required (%d: %d)",FCT,i+1,data[i]);
         else wblog(F_L,"ERR %s() "
            "index out of bounds (%d: %d/%d)",FCT,i+1,data[i],offset
         );
      }
   }

   return *this;
};

template <class T>
template <class T2>
void groupIndex<T>::getRecs0(
   const wbMatrix<T2> &AB, wbMatrix<T2> &X, T m
 ) const {

   widx_t i0=0; if (swidx_t(m)<0) { m=-m; i0=AB.dim2-m;  }

   if (swidx_t(i0)<0 || i0+m>AB.dim2) wblog(FL,
      "ERR %s() size out of bounds (%d+%d/%d)",FCT,i0+m,AB.dim2);

   X.init(D.len,m); if (!D.len || !m) return;

   if (I0.len!=D.len+1) wblog(FL,
      "ERR %s() unexpected size %d/%d",FCT,I0.len-1,D.len);
   if (AB.dim1!=P.len) wblog(FL,"ERR %s() size mismatch "
      "(%dx%d <> %dx%d)",FCT,P.len, m, AB.dim1, AB.dim2);

   for (widx_t i=0; i<D.len; ++i) {
      X.recSetP(i,AB.ref(I0.data[i],i0));
   }
};

int ctrIdx::init(const char *F, int L, const char *s, unsigned offset) {

   int i, n=Str2Idx(FL,s,(wbvector<unsigned>&)(*this),offset,'x');

   if (!n) {
      if (F) wblog(F,L,"ERR %s() got empty/null string !?",FCT);
      return -1;
   }

   if (n<0) { int m=0; i=-n-1; 

      for (; s[i]; ++i) { if (!isspace(s[i])) {
         if (CTR_IS_SEP_ANY(s[i])) { ++m; } else break;
      }}
      if (m<=1 && CTR_IS_CONJ(s[i])) { conj=1;
         for (++i; s[i]; ++i) { if (!isspace(s[i])) break; }
      }

      if (s[i]) { if (F) wblog(F,L,"ERR %s() "
         "invalid ctr-index set '%s' (%d: %c<%d>)",FCT,s,i+1,s[i],s[i]);
         return -2;
      }
   }

   if (!len) { if (F) wblog(F,L,
      "ERR %s() got empty ctr-index set '%s' !?",FCT,s);
   }
   else if (!wbvector<unsigned>::isUnique()) { if (F) wblog(F,L,
       "ERR %s() invalid index set '%s' [%s]",
        FCT,STR_(this),wbvector<unsigned>::toStr().data);
      return -3;
   }
   else if (F) { 
      if (len>32 || anyGT(62)) wblog(F,L,"WRN "
      "got ctr-index [%s] (len=%d)",wbvector<unsigned>::toStr().data,len);
   }

   return len;
};

ctrIdx& ctrIdx::init(const char *F, int L, const mxArray *a) {

   if (!a) wblog(FL,"ERR %s() got null mxArray",FCT);
   if (mxIsChar(a)) {
      wbstring s(a); int e=init(F,L,s.data);
      if (e<=0) wblog(FL,
         "ERR %s() '%s' => [%s] @ e=%d !?",FCT,s.data,STR_(this),e);
   }
   else {
      wbvector<unsigned>::init(F,L,a);
      if (anyEqual(0)) wblog(FL,
         "ERR %s() contraction index must be 1-based",FCT);
      (*this)-=1; 
      conj=0;
   }
   return *this;
};

template<class TQ, class TD>
ctrIdx& ctrIdx::initOp(const char *F, int L,
   itag_ t1, const QSpace<TQ,TD> &B, const char *sopt) {

   unsigned l=0, r=-1;

   if (!sopt || (l=strlen(sopt))>8 || strncmp(sopt,"-op",3)) wblog(F_L,
      "ERR %s() invalid ctr-string '%s'",FCT,sopt?sopt:"");

   if (B.isOperator(&r)<=0 || r>3 || r!=B.itags.len) wblog(F_L,
	  "ERR %s() got '%s' with QSpace %s !?",FCT,sopt,STR(B.itags));

   init(1);
   if (sopt[l-1]=='*') { data[0]=0; conj=1; } else { data[0]=1; }

   newtags.initOp(t1,B.itags,sopt+1);

   return *this;
};

int ctrIdx::extend2Perm(widx_t N, wbperm &P, char toend) const {

    widx_t i,k;
    wbvector<char> mark(N);

    for (i=0; i<len; ++i) {
       if (data[i]>=N) { sprintf(str,
          "%s:%d %s() index out of bounds (%d/%ld)",FLF,data[i],N);
          return 1;
       }
       if (mark[data[i]]++) { sprintf(str,
          "%s:%d %s() index not unique",FLF);
          return 1;
       }
    }

    P.init(N);

    if (!toend) { 
       for (k=0; k<len; ++k) { P.data[k]=data[k]; }
       for (i=0; i<N; ++i) { if (!mark[i]) P[k++]=i; }
    }
    else {
       for (k=i=0; i<N; ++i) { if (!mark[i]) P[k++]=i; }
       for (i=0; i<len; ++i) { P.data[i+k]=data[i]; }
    }

    return 0;
};

wbperm&
ctrIdx::extend2Perm(    unsigned ra, unsigned ma,
   const ctrIdx &icb__, unsigned rb, unsigned mb,
   const wbperm &pab, char conj, unsigned mc,
   wbperm &Pbc, 
   ctrIdx &kb,  
   ctrIdx &kcb, 
   char Bflag   
) const {

   ctrIdx ka,   
          kca;  

   unsigned i=0, l=0, m=0, rca=ra-len, rcb=rb-icb__.len, rc=rca+rcb;

   ctrIdx ica,icb; wbperm p,ipab;
   if (Bflag)
        { icb__.sort(icb,p); this->wbvector<unsigned>::permute(ica,p); }
   else { this->sort(ica,p); icb__.wbvector<unsigned>::permute(icb,p); }

   if (pab.len)
        ipab.init(pab,'i');
   else ipab.init(rc);

   if (len!=icb.len || ra<len || rb<icb.len) wblog(FL,"ERR %s() "
      "invalid ctrIdx data (%d/%d, %d/%d) !?",FCT,ra,len,rb,icb.len);
   if (ra>127 || rb>127 || rc>127) wblog(FL, 
      "ERR %s() got rank %d/%d/%d !?",FCT,ra,rb,rc);
   if (ipab.len!=rc) wblog(FL,"ERR %s() "
      "invalid permutation pab=[%s] (%d)",FCT,STR(pab),rc);

   ica.invert(ra,ka); 
   icb.invert(rb,kb); 

   kca.initT(rca, ipab.data,     !conj);
   kcb.initT(rcb, ipab.data+rca, !conj);

   if (Bflag) {
      wbperm Pa(ra+(mb?1:0)+(mc?1:0)); m=ra;

      for (; l<len; ++l) { Pa[l]=ica[l]; } 
      if (mb) Pa[l++]=(m++); 

      kca.Sort(p); ka.Permute(p);
      for (i=0; i<ka.len; ++i, ++l) { Pa[l]=ka[i]; }
      if (mc) Pa[l++]=(m++); 

      if (l!=Pa.len || !Pa.isValidPerm()) wblog(FL, 
         "ERR %s() got invalid Pa=[%s] (%d/%d)",FCT,STR(Pa),l,Pa.len);
      Pa.invert(Pbc);
   }
   else {
      wbperm Pb(rb+(ma?1:0)+(mc?1:0)); m=rb;

      for (; l<len; ++l) { Pb[l]=icb[l]; } 
      if (ma) Pb[l++]=(m++); 

      kcb.Sort(p); kb.Permute(p);
      for (i=0; i<kb.len; ++i, ++l) { Pb[l]=kb[i]; }
      if (mc) Pb[l++]=(m++); 

      if (l!=Pb.len || !Pb.isValidPerm()) wblog(FL, 
         "ERR %s() got invalid Pb=[%s] (%d/%d)",FCT,STR(Pb),l,Pb.len);

      Pb.invert(Pbc); 
      ka.save2(kb);
      kca.save2(kcb);
   }

   return Pbc;
};

int ctrIdx::checkUniqueS(const char *F, int L, const unsigned r) const {

   if (len) {
      unsigned i=1;
      for (; i<len; ++i) {
         if (data[i-1]>=data[i]) {
             if (data[i-1]>data[i]) { if (F) wblog(F,L,
                "ERR %s() ctrIdx not sorted [%s]",FCT,STR_(this));
                return 1;
             }
             if (data[i-1]==data[i]) wblog(F_L,
             "ERR %s() ctrIdx not sorted [%s]",FCT,STR_(this));
         }
      }

      if (int(r)>0) { 
         if (data[i-1]>=r) wblog(F_L,
         "ERR %s() ctrIdx out of bounds (%s; %d) !?",FCT,STR_(this),r);
      }
   }

   return 0;
};

ctrIdx& ctrIdx::invert(unsigned r, ctrIdx &I) const {

    wbvector<char> M(r); char *m=M.data;
    unsigned i=0, l;

    for (; i<len; ++i) { l=data[i];
       if (l>=r) wblog(FL,
          "ERR %s() index out of bounds (%d: %d/%d)",FCT,i,l,r);
       if ((++m[l])>1) wblog(FL,
          "ERR %s() index not unique: %s",FCT,toStr().data
       );
    }

    I.init(r-len);
    if (I.len) {
       for (l=i=0; i<r; ++i) { if (!m[i]) I.data[l++]=i; }
    }

    I.conj=conj;

    return I;
};

wbstring ctrIdx::toStr(char xflag) const {
   wbstring s; 

   if (!xflag && this->anyGT(60)) { 
      unsigned l, n=4*len+4; s.init(n); 

      l=snprintf(s.data,n,"%s%s",
         ((wbvector<unsigned>&)(*this)+1).toStrf("%d",",").data,
         conj? ";*":"");
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)",FCT,l,n
      );
   }
   else {
      unsigned i=0, n=len+(conj ? 2 : 1); s.init(n);
      for (; i<len; ++i) {
         if (data[i]< 9) { s.data[i]='1'+ data[i];     } else
         if (data[i]<35) { s.data[i]='A'+(data[i]- 9); } 
         else            { s.data[i]='a'+(data[i]-35); } 
      }
      if (conj) { s.data[i++]='*'; }
      s.data[i]=0;
   }

   return s;
};

itag_& itag_::init(const char* F, int L, const char *s, unsigned n) {

   if (!s || !s[0] || !n) { t=0; return *this; }

   unsigned i=0, l=0, n_=ITAG_LEN;
   char *q=(char*)(&t); t=0;

   if (int(n)<0 || n>n_) n=n_;

   while (s[i]==' ') { ++i; };
   if (i) wblog(FL,"WRN %s() ignoring leading blanks (%d)",FCT,i);

   for (; i<n && s[i]; ++i) {
      if (IS_CHAR_ITAG(s[i])) { q[l++]=s[i]; }
      else break;
   }

   while (s[i]==CC_ITAG) { ++i;
      q[0] ^= char(128); 
   }

   if (i<n && s[i]) {
      if (i>=n_)
           sprintf(str,"invalid itag `%s' (len=%li/%d)",s,strlen(s),n);
      else sprintf(str,"invalid itag `%s' (len=%li/%d; i=%d)",s,strlen(s),n,i);
      if (F) wblog(FL,"ERR %s",str); else wblog(F_L,"WRN %s",str);
      t=0; return *this;
   }
   if (!t && F) wblog(F,L,"WRN %s() invalid itag '%s'",FCT,s);

   return *this;
};

template<class TQ, class TD>
itag_& itag_::init(const char *F, int L,
   const wbvector< QSpace<TQ,TD> > &A,
   wbindex ia 
){
   t=0; if (!A.len) return *this;

   unsigned k=0, l=0; itag_ x;

   if (!ia.len) {
      unsigned r=-1, rx;
      for (; k<A.len; ++k) { rx=r;
         if (int(r=A[k].isOperator(&rx,'x'))<=0) wblog(FL,
            "ERR %s() invalid rank-%d operator \n(%d: %d/%d, `%s')",
            FCT,A[k].itags.len,k+1,r,rx,A[k].otype2Str().data
         );
      }
      ia.Index(r); k=0;
   }
   if (!ia.len) wblog(FL,"ERR %s() got empty index set !?",FCT);

   for (; k<A.len; ++k) {
      const iTags &it=A[k].itags;
      for (l=0; l<ia.len; ++l) {
         if (ia[l]>=it.len) wblog(F_L,
            "ERR %s() index out of bounds (%d/%d)",FCT,ia[l],it.len);
         (x=it[ia[l]]).deConj();

         if (!t) t=x.t; else
         if (x.t && x.t!=t) wblog(F_L,"ERR %s() got itags "
            "mismatch (%s [%s]; %s)",FCT, A[k].itags.toStr().data,
            (ia+1).toStr().data, this->toStr().data
         );
      }
   }
   return *this;
};

itag_& itag_::SetSI(const char *F, int L,
   const char *tag, unsigned k, char conj, char nfmt) {

   unsigned l, n=16; char s[n];
   if (!tag) wblog(F_L,"ERR %s() got null tag",FCT);

   if (nfmt<=1) {
      l=snprintf(s,--n,"%s%d",tag,k); 
   }
   else { l=9; memset(s,' ',l); 
      l+=snprintf(s+l,n-l,"%%s%%0%d",nfmt); if (l<n) { n=9;
      l=snprintf(s,n,s+n,tag,k); } 
   }
   if (l>=n) wblog(F_L,
      "ERR %s() string out of bounds (%s: %d/%d)",FCT,s,l,n);

   if (conj=='*' || (conj<0 && isConj())) {  
      s[l]=CC_ITAG; s[l+1]=0;
   }
   else if (conj>0) wblog(F_L,"ERR %s() got conj=%d",FCT,conj);

   init(s);

   return *this;
};

#ifndef NOMEX

itag_& itag_::init(const char* F, int L, const mxArray *a) {

   if (!a) { t=0; return *this; }
   if (mxIsChar(a)) { wbstring s(a); init(F,L,s.data); }
   else wblog(F_L,
      "ERR %s() invalid QSpace.itag upon input (%s)",
       FCT,mxGetClassName(a)
   );
   return *this;
};

#endif

char itag_::sameAs(const itag_ &B, char lflag) const {

   char q=0;  
   IDT c128=128, a=t, b=B.t;
   bool cflag = (a & c128) != (b & c128); 

   a &= ~c128; 
   b &= ~c128;

   if (a==b) { q=(a ? 3 : 4); } else
   if (lflag) {
      if (!b) { q|=1; } else
      if (!a) { q|=2; }
   }

   if (cflag) { q=-q; }
   return q;
};

unsigned itag_::to_str(char *s, unsigned len,
   char cflag 
 ) const {

   unsigned e=0, l=0, i=0, m=0, n=ITAG_LEN, nx=2;
   const char *st=(const char *)(&t), c128=char(128);

   if (int(len)<int(n)) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)",FCT,len,ITAG_LEN);

   for (; i<n && st[i]; ++i) { s[i]=st[i];
      if (s[i]<0) { s[i]+=c128;
         if (i) { ++m; }; if (i>=nx) { ++e; }
         if (!s[i]) { break; }
      }
      if (!IS_CHAR_ITAG(s[i])) { e+=10; }
   }; s[i]=0;

   for (l=i; i<n; ++i) {
      if (st[i]) {
         if (st[i]==c128)
              { if (i>=nx) { e+=100; } else if (i) { ++m; }}
         else { e+=100; }
      }
   }

   if (e) wblog(FL,
      "WRN %s() got non-char itag (%s; e=%d) !?",FCT,s,e);

   if (l+m+3>=len) wblog(FL,"ERR %s() "
      "string out of bounds (%d+%d+3=%d/%d; %s)",FCT,l+m,l+m+3,len,s);

   if (m) { s[l++]='\''; 
      for (i=1; i<n; ++i) { if (st[i]<0) { s[l++]='0'+i; }}
   }
   if (st[0]<0 && cflag) { s[l++]=CC_ITAG; }

   s[l]=0;

   return l;
};

int itag_::isValid() const {

   unsigned l=0, i=0, n=ITAG_LEN, nx=2;
   const char *s=(const char *)(&t), x128=~char(128); 
   char x;

   for (; i<n; ++i) { x=(s[i] & x128);
      if (x) { if (!IS_CHAR_ITAG(x)) return 0; } else break;
      if (s[i]<0 && i>=nx) return 0;
   }; l=i; 

   for (; i<n; ++i) { if (s[i]) {
      if (i<nx) { if (s[i] & x128) return 0; }
      else return 0;
   }}

   return (l ? 2: 1);
};

bool itag_::isConj_x(const itag_ &x, char lflag) const {

   if (lflag<=0) { return isConj(x); }
   if (lflag>3) {
      if (lflag=='l') lflag=1; else
      if (lflag=='L') lflag=3; else
      wblog(FL,"ERR %s() invalid lflag=%d !?",FCT,lflag);
   }

   unsigned i=1;
   const char *a=(const char*)(&t), *b=(const char*)(&x.t);

   if (!(lflag & 2)) { 
      if ((a[0]^char(128)) != b[0]) return 0;
      for (; i<ITAG_LEN; ++i) { if (a[i]!=b[i]) break; }
   }
   else { 
      if (tolower(a[0]^char(128)) != tolower(b[0])) return 0;
      for (; i<ITAG_LEN; ++i) { if (tolower(a[i])!=tolower(b[i])) break; }
   }

   if (i<ITAG_LEN) {
      if (!(lflag & 1) ||
         (a[i] && a[i]!='~') || (b[i] && b[i]!='~')) return 0;
      for (++i; i<ITAG_LEN; ++i) { if (a[i] || b[i]) return 0; }
   }

   return 1;
};

itag_& itag_::AppendChar(char q, const char *qs) {

   unsigned i=0, n=ITAG_LEN;
   char *s=(char*)(&t);

   if (!IS_CHAR_ITAG(q)) wblog(FL,
      "ERR %s() invalid char %c<%d>",FCT,q,q);

   for (; i<n; ++i) {
      if (!( s[i] & ~char(128) )) break; 
   }

   if (!i) { s[i]+=q; } 
   else {
      char x=(s[i-1] & ~char(128));
      if (qs && strchr(qs,x)) { --i; } 
      else if (i==n) 
      wblog(FL,"ERR %s() itag out of bounds (%s)",FCT,STR(*this));

      s[i] = (s[i] & char(128)) + q;
   }

   return *this;
};

itag_& itag_::PrependChar(char q, const char *qs) {

   unsigned i=0, n=ITAG_LEN;
   char *s=(char*)(&t), x=(s[0] & ~char(128));

   if (!IS_CHAR_ITAG(q)) wblog(FL,
      "ERR %s() invalid char %c<%d>",FCT,q,q);

   if (!qs || !strchr(qs,x)) { 
      for (; i<n; ++i) {
         if (!( s[i] & ~char(128) )) break; 
      }
      if (i==n) wblog(FL,
         "ERR %s() itag out of bounds (%s)",FCT,toStr().data);
      for (--i; i>0; --i) { s[i+1]=s[i]; }
      s[1] = (s[0] & ~char(128));
   }

   s[0] = (s[0] & char(128)) + q;

   return *this;
};

int itag_::CheckFirstChar(char q, const char *qs) {

   char *s=(char*)(&t), x=(s[0] & ~char(128));

   if (!IS_CHAR_ITAG(q)) wblog(FL,
      "ERR %s() invalid char %c<%d>",FCT,q,q);

   if (qs && strchr(qs,x)) {
      s[0] = (s[0] & char(128)) + q;
      return 1;
   }

   return 0; 
};

unsigned iTags::init(const char *F, int L, const char *s) {

   if (s && s[0]) {
      unsigned i,l=0, i0=0, n=0;
      for (i=0; s[i]; ++i) {
         if (!IS_CHAR_ITAGS(s[i])) wblog(F_L,
            "ERR %s() invalid itag >%s<\n"
            "[required: printable string with no space]",FCT,s);
         if (IS_ITAGS_SEP(s[i])) ++l;
      }
      wbvector<itag_>::init(++l);

      for (l=0, n=i, i=0; i<=n; ++i) {
         if (IS_ITAGS_SEP(s[i]) || i==n) {
            data[l].init(F_L,s+i0,i-i0);
            ++l; i0=i+1;
         }
      }
   }
   else wbvector<itag_>::init();

   return len;
};

iTags& iTags::init_alpha(
   unsigned r, unsigned s, 
   const char *tag,  
   char t,           
   unsigned k        
){
   unsigned i=0, l, n=10; char st[n];

   init(r+s);
   for (; i<len; ++i) {          
      l=snprintf(st,n,"%s%c%d",tag ? tag : QS_ITAG_,t+i,k);
      if (i>=r && ++l<n) { st[l-1]=CC_ITAG; st[l]=0; }
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%s: %d/%d)",FCT,st,l,n);
      data[i].init(st);
   }

   return *this;
};

#ifdef LOAD_CGC_QSPACE

unsigned iTags::init_qdir(const char *s) {

   QDir qdir(FL,s); init(qdir.len);

   for (unsigned i=0; i<qdir.len; ++i) {
      if (qdir[i]<0) data[i].Conj(); 
   }

   return qdir.len;
};

#endif

#ifndef NOMEX
unsigned iTags::init(const char *F, int L, const mxArray *a) {

   if (!a || mxIsEmpty(a)) { init(); return len; }

   if (mxIsChar(a)) { wbstring sz(a); init(sz.data); }
   else if (mxIsCell(a)) {
      unsigned i=0, n=mxGetNumberOfElements(a);
      wbvector<itag_>::init(n);
      for (; i<n; ++i) data[i].init(F_L,mxGetCell(a,i));
   }
   else wblog(FL,"ERR invalid QSpace.itag upon input");

   return len;
};
#endif

iTags& iTags::initOp(
    itag_ t, const iTags &b, const char *ot) {

   if (!b.isOp()) wblog(FL,"ERR %s() invalid usage (%s)",FCT,STR(b));
   init(b);

   if (!t.isEmpty()) { t.deConj();
      for (unsigned i=0; i<2; ++i) {
          if (data[i].isConj()) { (data[i]=t).Conj(); }
          else { data[i]=t; }
      }
   }

   if (len>2 && data[2].isEmpty()) {
      unsigned i=0;
      char s[10]; strncpy(s, ot?ot:"op", 8); s[8]=0;

      for (; s[i]; ++i) { if (s[i]==CC_ITAG) { s[i]=0; break; }}

      if (data[2].isConj()) { s[i]=CC_ITAG; s[i+1]=0; }
      data[2].init(s);
   }

   return *this;
};

unsigned iTags::Init( 
   const char *F, int L, const char *tag, unsigned k, char lflag
){
   size_t l=0, n=64; char s[n]; s[0]=0;
   iTags Ix(*this);  

   if (!tag) wblog(F_L,"ERR %s() invalid tag (empty)",FCT);
   if (!strncmp(tag,"nrg:",4)) {
      const char *t=tag+4;
      if (strlen(t)!=2) wblog(FL,"ERR %s() invalid tag `%s'",FCT,tag);

      if (int(k)<0) wblog(FL, 
         "ERR %s() invalid k=%d for '%s'",FCT,k,tag);

      if (t[0]=='H' && strchr("KD",t[1])) { 
         l=snprintf(s,n,"%c%02d;%c%02d*",t[1],k,t[1],k);
      }
      else if (t[0]=='A' && strchr("KD",t[1])) { 
         if (k) 
              { l=snprintf(s,n,"K%02d;%c%02d*;s%02d",k-1,t[1],k,k); }
         else { l=snprintf(s,n,  "L00;%c%02d*;s%02d",    t[1],k,k); }
      }
      else wblog(FL,"ERR %s() invalid tag `%s'",FCT,tag);
   }
   else if (!strncasecmp(tag,"op:",3) && int(k)<0) {
      int i=0, n1=0, n2=-1; const char *t=tag+3;
      unsigned r=-k; 

      if (lflag<0) {
         lflag = (tag[0]=='O'? 1:0) + (tag[1]=='P'? 2:0);
      }

      if (r<2 || (lflag ? r>4 : r>3)) wblog(FL,
         "ERR %s() invalid rank-%d operator",FCT,r);

      for (; t[i]; ++i) {
         if (t[i]==':') {
            n1=i++; while (t[i]) { ++i; }; n2=i-n1-1;
            if (n1>8 || n2>8) wblog(FL,"ERR %s() "
               "itag string out of bounds ('%s' %d,%d/8)",FCT,t,n1,n2);
            break;
         }
      }
      if (n2>=0) { char fmt[20];
         if (r>3) wblog(FL, 
            "ERR %s() got rank-%d QSpace (l=%d)",FCT,r,lflag);

         if (r==2) {
            snprintf(fmt,20,"%%%d.%ds;%%%d.%ds*",n1,n1,n1,n1);
            l=snprintf(s,n,fmt,t,t);
         }
         else {
            snprintf(fmt,20,"%%%d.%ds;%%%d.%ds*;%%s*",n1,n1,n1,n1);
            l=snprintf(s,n,fmt,t,t,t+n1+1);
         }
      }
      else if (r==3)
           { l=snprintf(s,n,"%s;%s*;op*",t,t); }
      else { l=snprintf(s,n,"%s;%s*",t,t); } 

      if (r==3 && !data[2].isConj() && l>1 && l<n) { s[--l]=0; }
   }
   else wblog(F_L,"ERR %s() unexpected tag '%s'",FCT,tag);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%s/%d)",FCT,l,n);
   n=0; 

   init(F,L,s);

   if (len<Ix.len) {

      for (unsigned i=0; i<len; ++i) {
         if (data[i].isConj(Ix.data[i],0)) { ++n; }
         else { Ix[i]=data[i]; }
      }
      if (!n) { Ix.save2(*this); }
   }
   else if (Ix.len && !Ix.sameConj(*this)) {
      n=99; 
   }

   if (n) wblog(FL,"ERR inconsistent conj labels "
      "(%s => %s; e=%d)",STR(Ix), STR_(this), n);
   return len; 
};

iTags& iTags::Init(unsigned r, int k, const char **ss) {

   unsigned i=0, l, n=9; char st[n+1]; 
   init(r);

   for (; i<r; ++i) { l=snprintf(st,n,"%s", ss[i]);
      if (k>=0 && l<n) {
         if (l>1 && st[l-2]==CC_ITAG) { --l;
              l+=snprintf(st+l,n-l,"%d%c",k,CC_ITAG)-1; }
         else l+=snprintf(st+l,n-l,"%d",k);
      }
      if (l>=n) { st[n]=0; wblog(FL,
         "ERR %s() string out of bounds (%s /%d)",FCT,st,n); }
      data[i].init(st);
   }

   return *this;
};

unsigned iTags::Set(const char* F, int L,
   const char *tag, unsigned r,  
   const char *tom, 
   unsigned m,      
   bool conjOM      
){
   unsigned rval=0, i,l;
   if (!len && int(r)<=0 && int(m)<=0) { return rval; }

   if (int(r)<0) {
      if (int(m)<=0) { r=len; } else
      if (m+2<len) { r=len-m; } else  
      wblog(FL,"ERR %s() requesting %d OM legs for rank-%d",FCT,m,len-m);
   } else
   if (r>len) wblog(FL, 
      "ERR %s() rank out of bounds (%d/%d)",FCT,r,len);

   if (int(m)<0) { if (tom) { m=len-r; }}
   else {
      if (m && r< 2) wblog(FL,"ERR %s() "
         "requesting %d OM leg%s for rank-%d",FCT,m,m==1?"":"s",r);

      if ((l=r+m)<len) wblog(FL,
         "ERR %s() asking for lower rank (%d+%d / %d)",FCT,r,m,len);
      if (l>len) { this->Resize(l); }
   }

   if (!tag || (r<len && int(m)>0 && !tom)) wblog(FL,
      "ERR %s() got null strings (%p / %p)",FCT,tag,tom);
   if (int(m)<0) { m=0; }

   for (i=0; i<r; ++i) {
      if (data[i].isEmpty()) { ++rval;
         data[i].SetSI(FL,tag,i+1); 
      }
   }

   for (i=0; i<m; ++i) {
      if (!data[r+i]) { ++rval;
         data[r+i].SetSI(FL,tom,i+1, conjOM ? '*' : 0); 
      }
      else wblog(FL,
      "ERR %s() non-empty OM itags (%s; %d+%d)",FCT,STR_(this),r,m);
   }

   return rval;
};

wbstring iTags::toStr(char vflag) const {

   unsigned k, n=ITAG_LEN;
   size_t l=0, N=(n+2)*len + (vflag ? 16 : 0); char s[N+1];

   for (k=0; k<len; ++k) { if (k) s[l++]=';';
      l+=data[k].to_str(s+l,N);
      if (l>=N) wblog(FL,"ERR %s() string out of bounds (%s)",FCT,s);
   }

   if (vflag) {
      l+=snprintf(s+l,N-l," (rank-%ld)",len);
      if (l>=N) wblog(FL,"ERR %s() string out of bounds (%s)",FCT,s);
   }
   else s[l]=0;

   return s;
};

int iTags::getCtr( 
   const char *F, int L,  const iTags &B,
   ctrIdx &ia, ctrIdx &ib,
   iTags *C, 
   wbvector<char> *ma_, wbvector<char> *mb_
 ) const {

   int nc=0; 
   wbvector<char> ma(len), mb(B.len); 
   unsigned i,j, it=0, l=0, cflag=(ia.conj ^ ib.conj);
   int q, w=0;

   for (i=0; i<  len; ++i) { if (  data[i].GotFlag(1)) { ma[i]=-1; }}
   for (i=0; i<B.len; ++i) { if (B.data[i].GotFlag(1)) { mb[i]=-1; }}

   for (; it<2; ++it) {
      q=(cflag ? SGN(+1) : SGN(-1));

      for (i=0; i<  len; ++i) { if (ma[i]>=0) {
      for (j=0; j<B.len; ++j) { if (mb[j]>=0) {
         if (SGN(data[i].sameAs(B.data[j]))==q) {
            if (!ma[i] && !mb[j])  
                 { ma[i]=j+1; mb[j]=i+1; ++nc; } 
            else { ++w; } 
         }
      }}}}
      if (nc) { break; } else { cflag=!cflag; } 
   }

   ia.wbvector<unsigned>::init(nc);
   ib.wbvector<unsigned>::init(nc);

   if (nc) {
      for (i=0; i<len; ++i) { if (ma[i]>0) {
         ia[l]=i;
         ib[l]=ma[i]-1; ++l;
      }}

      if (it) { wblog(F_L,"WRN contract::match() "
         "check missing conj flag\nA: %-15s -> %s\nB: %-15s -> %s",
          STR_(this), STR(ia), STR(B), STR(ib));
         mexIssueWRN("applying conj(A) flag");
         ia.Conj(); 
      }
      if (w) { wblog(F_L,
         "WRN %s() non-unique itags `%s'/`%s'",FCT,STR(B),STR_(this));
         sprintf(str,"using %s / %s",STR(ia),STR(ib));
         mexIssueWRN(str);
      }
   }
   else if (F || !L) wblog(FL, 
     "ERR %s() got itag mismatch (conj=%d,%d)\n(%s <> %s)",
      FCT, ia.conj, ib.conj, STR_(this), STR(B)
   ); 

   if (C) { C->init( len + B.len - 2*nc );
   if (C->len) { l=0; itag_ *c=C->data;

      for (i=0; i<len; ++i) { if (ma[i]<=0) {
         c[l]=  data[i]; if (ia.conj) { c[l].Conj(); }; ++l;
      }}
      for (i=0; i<B.len; ++i) { if (mb[i]<=0) {
         c[l]=B.data[i]; if (ib.conj) { c[l].Conj(); }; ++l;
      }}
      if (l!=C->len) wblog(FL,"ERR %s() %d/%d",FCT,l,C->len);
   }}

   if (ma_) { 
      for (i=0; i<ma.len; ++i) { if (ma[i]<0) { ma[i]=0; }}
      ma.save2(*ma_);
   }
   if (mb_) { 
      for (i=0; i<mb.len; ++i) { if (mb[i]<0) { mb[i]=0; }}
      mb.save2(*mb_);
   }

   return nc;
};

int iTags::UpdateItagsCtr(const char *F, int L,
   const ctrIdx &ia, iTags &B,
   const ctrIdx &ib,
   const char *tag
){
   int q, s=0, e=0; unsigned i=0, k=0, l;

   if (!ia || ia.len!=ib.len) wblog(F_L,
      "ERR %s() invalid usage (ia.len=%d/%d)",FCT,ia.len,ib.len);
   if (ia.len>15) wblog(F_L, 
      "ERR %s() contracting %d>15 indices not supported",FCT,ia.len);

   if (ia.conj) { this->Conj(); }
   if (ib.conj) {     B.Conj(); }

   for (; i<ia.len; ++i) {
      if (ia[i]>=  len) wblog(F_L,
         "ERR %s() index out of bounds (A: %s / %d)",FCT,STR(ia),  len);
      if (ib[i]>=B.len) wblog(F_L,
         "ERR %s() index out of bounds (B: %s / %d)",FCT,STR(ib),B.len);

      itag_ &a=this->data[ia[i]], &b=B.data[ib[i]];
      q=a.sameAs(b,'l');

      if (!q && (l=strlen(QS_ITAG_))) {
         if (!strncmp((const char*)a,QS_ITAG_,l)) { a.init_tag(); } else
         if (!strncmp((const char*)b,QS_ITAG_,l)) { b.init_tag(); }
         else { l=0; }

         if (!l || !(q=a.sameAs(b,'l'))) { 
            e+=1; continue;
         }
      }

      if (i) { if (SGN(q)!=s) {
         e+=(1<<4); continue;
      }}
      else { s=SGN(q); }

      if ((q=abs(q))==3) { continue; } 

      if (q==1) { b.tSet(a); } else
      if (q==2) { a.tSet(b); } else
      if (q==4) { 
         unsigned n=9; char st[n];
         if (!tag) wblog(F_L,"ERR %s() got null tag (%p)",FCT,tag);
         l=snprintf(st,n,"%s%d",tag,++k); 
         if (l>=n) { st[n-1]=0; wblog(FL,
            "ERR %s() string out of bounds (%s: %d/%d)",FCT,st,l,n); }
         a.Set(st); b.Set(st);
      }
      else wblog(F_L,"ERR %s() q=%d",FCT,q);
   }

   if (s>=0) { e|=(1<<8); }

   if (e && (F || L)) {
      sprintf(str,"'%s' @ %s / '%s' @ %s (e=%d)",
          STR_(this), STR(ia), STR(B), STR(ib), e);
      if (e&15) wblog(F_L,"ERR %s() itag mismatch %s",FCT,str);
      else if (unsigned((e>>=4)&15) + 1 >= ia.len)
           wblog(F_L,"ERR %s() missing conj-flag %s", FCT,str);
      else wblog(F_L,"ERR %s() itag conj mismatch %s",FCT,str);
   }

   if (!e) {
      FlagItagsCtr(F_L,ia,B,ib,0,-1);
   }

   return e;
};

int iTags::FlagItagsCtr(
   const char *F, int L, const ctrIdx &ia,
   iTags &B, 
   const ctrIdx &ib,
   char all, 
   char sgn  
) const {    

   int rval=0;
   unsigned i,j;
   wbvector<char> ma(len), mb(B.len);

   if (sgn)
        { sgn=SGN(sgn); } 
   else { sgn=((!ia.conj) ^ (!ib.conj) ? SGN(+1) : SGN(-1)); }

   for (i=0; i<ia.len; ++i) {
      if (ia[i]>=ma.len || ++ma[ia[i]]!=1) wblog(F_L,
      "ERR %s() invalid ctr-index A: %s /%d",FCT,STR(ia),ma.len);
   }
   for (i=0; i<ib.len; ++i) {
      if (ib[i]>=mb.len || ++mb[ib[i]]!=1) wblog(F_L,
      "ERR %s() non-unique ctr-index B: %s /%d",FCT,STR(ib),mb.len);
   }

   if (!all) { 
      for (i=0; i<ma.len; ++i) {
         if (!ma[i] &&   data[i].isEmpty()) { ma[i]=-9; }}
      for (i=0; i<mb.len; ++i) {
         if (!mb[i] && B.data[i].isEmpty()) { mb[i]=-9; }}
   }

   for (j=0; j<B.len; ++j) { if (!mb[j]) {
   for (i=0; i<  len; ++i) { if (!ma[i] &&
       SGN(data[i].sameAs(B.data[j]))==sgn) {
          B.data[j].SetFlag(1); ++rval; break;
       }}
   }}

   if (all) { ctrIdx ia_,ib_; 
      this->getCtr(F,L,B,ia_,ib_);
      if (ia_.len!=ib.len) { wblog(F_L,"ERR %s() ctr-index mismatch\n"
         "   A: %s | %s -> %s\n"
         "   B: %s | %s -> %s", FCT,
         STR_(this),STR(ia),STR(ia_), STR(B),STR(ib),STR(ib_));
      }
   }
   return rval;
};

int iTags::findRegEx(const char *pat) const {

   regex_t r; int e, k=-1;

   unsigned i, m=0, l=ITAG_LEN+2; 
   char s[l+1]; s[l]=0; 
   wbstring ptx;

   if (!pat || !pat[0]) wblog(FL,
      "ERR %s() got invalid regexp (empty)",FCT);

   if (strstr(pat,"\\d") && strlen(pat)<32) { unsigned j=0; str[0]=0;
      for (i=0; pat[i]; ++i) {
         if (pat[i]=='\\' && pat[i+1]=='d') { strcpy(str+j,"[0-9]"); j+=5; ++i; }
         else if (pat[i]=='[') {
            if (strstr(pat+i,"\\d")) wblog(FL,
               "WRN failed to translate non-posix \\d -> [0-9]\n"
               "in itag regexp '%s'",FCT,pat);
            strcpy(str+j,pat+i); break;
         }
         else { str[j++]=pat[i]; }
      }
      str[j]=0; ptx=str;
      pat=ptx.data;
   }

   if ((e=regcomp(&r,pat, REG_EXTENDED | REG_NOSUB))) {
      return -1;
   }

   for (i=0; i<len; ++i) {
       data[i].to_str(s,l);
       if (!regexec(&r,s,0,NULL,0)) { 
          if (!m || !data[k].sameAs(data[i])) { k=i; ++m; }
       }
   }
   regfree(&r); 

   if (m>1) {
      return -m;
   }

   return (k+1);
};

template <class T>
sparseIndex2D<T>& sparseIndex2D<T>::Setup(char use_rmaj) {

   if (IJ.dim2!=2) wblog(FL,
      "ERR %s() got invalid IJ (%s)",FCT,IJ.sizeStr().data);
   if (S.len!=2) wblog(FL,
      "ERR %s() got invalid S (len=%d)",FCT,S.len);

   const T *s=S.data;

   if (IJ.dim1<=1) { P.init(IJ.dim1);
      cIdx.init(s[1]+1);
      if (IJ.dim1) { T *idx=IJ.data; 
         if (idx[0]>=s[0] || idx[1]>=s[1]) wblog(FL,
            "ERR %s() index out of bounds (%ld,%ld) having (%dx%d) !?",
            FCT,idx[0],idx[1],s[0],s[1]);
         for (T i=idx[1]+1; i<=s[1]; ++i) cIdx[i]=1;
      }
      else P.init();
      return *this;
   }

   T j=0, l=0, lmax=IJ.dim1-1, *idx=IJ.data;

   if (!use_rmaj || s[0]>128 || (s[0]>1 && s[1]<(1<<24))) {
      cIdx.init(s[1]+1); rmaj=0;
      IJ.SortRecs(P,+1, -1);

      while (1) {
         for (; j<=idx[1]; ++j) { cIdx[j]=l; }
         for (; l<lmax && idx[1]==idx[3]; ++l, idx+=2) {
            if (idx[0]>=idx[2]) wblog(FL,
               "ERR %s() IJ not sorted or unique!",FCT);
         }; if (l>=lmax) { ++l; break; }

         if (idx[1]>idx[3]) wblog(FL,"ERR %s() IJ(:,2) not sorted!",FCT);
         if (idx[0]>=s[0]) wblog(FL,
            "ERR %s() index out of bounds (%ld/%ld) !?",FCT,idx[0],s[0]);
         ++l; idx+=2;
      }
   }
   else {
      cIdx.init(s[0]+1); rmaj=1;
      IJ.SortRecs(P,+1, +1);

      while (1) {
         for (; j<=idx[0]; ++j) { cIdx[j]=l; }
         for (; l<lmax && idx[0]==idx[2]; ++l, idx+=2) {
            if (idx[1]>=idx[3]) wblog(FL,
               "ERR %s() IJ not sorted or unique!",FCT);
         }; if (l>=lmax) { ++l; break; }

         if (idx[0]>idx[2]) wblog(FL,"ERR %s() IJ(:,2) not sorted!",FCT);
         if (idx[1]>=s[1]) wblog(FL,
            "ERR %s() index out of bounds (%ld/%ld) !?",FCT,idx[1],s[1]);
         ++l; idx+=2;
      }
   }

   for (; j<cIdx.len; ++j) { cIdx[j]=l; }

   return *this;
};

template <class T>
inline T sparseIndex2D<T>::cidx(T j) const {
#ifndef WB_SKIP_ASSERT
   if (j>=cIdx.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,j,cIdx.len);
#endif
   return cIdx.data[j];
};

template <class T>
inline T sparseIndex2D<T>::ridx(T ic) const {
#ifndef WB_SKIP_ASSERT
   if (ic>=IJ.dim1) wblog(FL,
      "ERR %s() index out of bounds (%d/%d)",FCT,ic,IJ.dim1);
#endif
   return IJ.data[2*ic]; 
};

#endif

