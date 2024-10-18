/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace MEX routines
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

#ifndef __WB_MEXLIB_CC__
#define __WB_MEXLIB_CC__

// also using mxPutAndDestroy() outside MEX file
// by writing data to file to capture debug data // Wb,Oct25,16

void mxPutAndDestroy(
    const char *F, int L, mxArray *a, const char *vn,
    const char* ws 
){
    int i=0;
    if (!ws ||!ws[0]) wblog(FL,"ERR %s() invalid workspace `%s'",FCT,ws);
    if (!strcmp(ws,"caller") || !strcmp(ws,"base")) { i=1; }

 #ifndef MATLAB_MEX_FILE
    if (i>0) { i=-1; }
 #endif
    if (i>0 && Wb::my_caller_tid!=omp_get_thread_num()) { i=-2;
    }

    if (i<0) {
       char tmp[64]; 
       tmpmatFL(F_L,vn && vn[0] && !strcmp(vn,"ans")? vn:ws, tmp,64);
       wblog(F_L,"WRN thread %d/%d cannot put `%s' to %s (i=%d)",
          omp_get_thread_num(), omp_get_num_threads(), vn?vn:"",ws,i);

       if (!strcmp(tmp,"caller") || !strcmp(tmp,"base")) wblog(FL,
          "ERR %s() invalid tmp=`%s' !?",FCT,tmp);
       mxPutAndDestroy(F,L,a,vn,tmp); 
       return;
    }

    #ifdef QS_USING_OMP
       Wb::ompGuard myLK(mxapi_lock); 
    #endif

    if (i) {
       i=mexPutVariable(ws,vn,a); 
       if (i) wblog(F_L,
          "ERR mexPutVariable('%s','%s',..) returned e=%d !?",ws,vn,i);
       if (F) wblog(F_L,"I/O putting `%s' to %s",vn,ws);
    }
    else {
       wbstring tmp;
       if (!ws || !ws[0] || !strcmp(ws,"tmpfile")) { 
          tmp.init(64); tmpmatFL(F_L,
             vn && vn[0] && strcmp(vn,"ans")? vn:"", tmp.data,tmp.len);
          ws=tmp.data; i=1;
       }
       else {
          unsigned j=0, k=0; 
          for (; ws[j]; ++j) { if (ws[j]=='.') { k=j+1; }}
          if (!k || strcasecmp(ws+k,"mat")) {
             tmp.init(k+6); snprintf(tmp.data,tmp.len,"%s.mat",ws);
             ws=tmp.data;
          }
       }
       if (F) wblog(F_L,"%s saving '%s' to '%s'",i?"TMP":"I/O",vn,ws);

       char x[2]; strcpy(x, Wb::isFile(ws)? "u":"w");
       Wb::matFile f(F_L,ws,x); 

       matPutVariable(f.mfp,vn,a); 
    }

    mxDestroyArray(a); 
};

int Mx::IsNumArray( 
   const char *F, int L, const mxArray* a,
   unsigned r, 
   char cflag, 
   char type   
){
   int rval=0;

   if (!a) {
      if (L) {
         strcpy(str,"got a=null");
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return (rval=-128);
   }

   if (mxIsDouble(a) ) { rval|= 2; } else 
   if (mxIsNumeric(a)) { rval|= 8; } else 
   if (mxIsChar(a)   ) { rval|=16; } else 
   if (mxIsLogical(a)) { rval|=32; }      

   if (rval && mxIsComplex(a)) { rval|=4; }

   if (type=='d') {
      if (!(rval & 2)) { 
         if (L) {
            sprintf_str("got data of type %s",mxGetClassName(a));
            if (F) wblog(F,L,"ERR %s() %s",FCT,str);
         }
         return (rval=-129);
      }
   }
   else if (type=='*') { 
      if (!(rval & 58)) { 
         if (L) {
            sprintf_str("got data of type %s",mxGetClassName(a));
            if (F) wblog(F,L,"ERR %s() %s",FCT,str);
         }
         return (rval=-130);
      }
   }
   else wblog(F_L,"ERR %s() invalid type=%d<%c>",FCT,type,type);

   if (mxIsSparse(a)) {
      if (L) {
         strcpy(str,"got sparse input");
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return (rval=-131);
   }

   if (!cflag && (rval & 4)) { 
      if (L) {
         strcpy(str,"got complex input");
         if (F) wblog(F,L,"ERR %s() %s",FCT,str);
      }
      return (rval=-132);
   }

   if (int(r)>0) {
      unsigned ra=mxGetNumberOfDimensions(a);
      if (ra>r) { 
         if (ra==2) {
            const size_t *D=mxGetDimensions(a);
            if (D[0]==1 || D[1]==1) { ra=1; }
         }
         if (ra>r) { 
            if (L) {
               sprintf_str("got invalid rank (r=%d/%d)",ra,r);
               if (F) wblog(F,L,"ERR %s() %s",FCT,str);
            }
            return (rval=-ra); 
         }
      }
   }
   rval|=1; 

   return rval;
};

template <class T>
Mx::Array<T>::Array(const wbvector<size_t> &S, char t_, char c_)
 : ax(0), data(0), len(0), trans(t_), cmplx(c_), mxref(0) {

   if (cmplx<0) { cmplx=ISCOMPLX_(T); }
   init(S.len,S.data);
};

template <class T>
int Mx::Array<T>::init_dptr() {

   data=0; if (!ax) { return -1; }

   if (mxref<=0) wblog(FL,"ERR %s() got ref=%d !?",FCT,mxref);
   if (mxIsSparse(ax)) wblog(FL,
      "ERR Mx::Array got invalid sparse `%s'",mxGetClassName(ax));

   mxClassID t=mxGetClassID(ax); int e=0;

   if (!mxIsComplex(ax)) {
      switch (t) { 
         case mxDOUBLE_CLASS  : safe_set_dptr( mxGetDoubles (ax));  break;
         case mxCHAR_CLASS    : safe_set_dptr( mxGetChars   (ax));  break;
         case mxINT8_CLASS    : safe_set_dptr( mxGetInt8s   (ax));  break;
         case mxUINT8_CLASS   : safe_set_dptr( mxGetUint8s  (ax));  break;
         case mxINT16_CLASS   : safe_set_dptr( mxGetInt16s  (ax));  break;
         case mxUINT16_CLASS  : safe_set_dptr( mxGetUint16s (ax));  break;
         case mxINT32_CLASS   : safe_set_dptr( mxGetInt32s  (ax));  break;
         case mxUINT32_CLASS  : safe_set_dptr( mxGetUint32s (ax));  break;
         case mxINT64_CLASS   : safe_set_dptr( mxGetInt64s  (ax));  break;
         case mxUINT64_CLASS  : safe_set_dptr( mxGetUint64s (ax));  break;
         case mxLOGICAL_CLASS : safe_set_dptr( mxGetLogicals(ax));  break;
         default: e=1; break; 
      }
   }
   else {
      if (t==mxDOUBLE_CLASS) safe_set_dptr( mxGetComplexDoubles(ax) );
      else { e=2; }
   }

   if (e && !mxIsNumChar(ax))
      wblog(FL,"ERR %s got invalid %s `%s' (e=%d)",xStr(),
      mxIsSparse(ax) ? "sparse":"type", mxGetClassName(ax),e);
   if (e==2 && !data) 
      wblog(FL,"ERR %s got complex `%s'",xStr(),mxGetClassName(ax));
   if (e) wblog(FL,"WRN %s not match found for %s",xStr(),TSTR(T));

   if (!e) { mxref+=8; } 

   return (e<0 ? e : (data ? 0 : 1));
};

template<class T>
Mx::Array<T>& Mx::Array<T>::init(
   unsigned ndim, size_t *dims, char t_, char c_) {

   if (ax) { mxDestroyArray(ax); ax=NULL; data=NULL; }; mxref=0;

   if (t_>=0) { trans=t_; }
   if (c_>=0) { cmplx=c_; }
   if ((rank=ndim)) {
      len=dims[0]; for (unsigned i=1; i<ndim; ++i) { len*=dims[i]; }
   } else len=0;

   MX_INIT_DATA(ndim,dims);
   init_xstr();

   if (ax && (cmplx ^ mxIsComplex(ax))) {
      if (cmplx) { cmplx=0; wblog(FL,
      "WRN %s ignoring cmplx=%d for %s",xStr(),mxGetClassName(ax)); }
      else wblog(FL,
      "ERR %s complex mismatch (%d/%d)",xStr(),cmplx,mxIsComplex(ax));
   }

   return *this;
};

template<class T>
Mx::Array<T>& Mx::Array<T>::init(const mxArray *a, char t_, char c_) {

   int e=0;
   if (ax) { if (mxref<=0) mxDestroyArray(ax); data=NULL; }
   mxref=1; len=rank=0;

   if (t_>=0) { trans=t_; }
   if (c_>=0) { cmplx=c_; }

   ax=(mxArray*)a; if (!a) { return *this; }

   e=init_dptr();
   if (e<0) wblog(FL,
      "ERR %s() %s <> %s => e=%d !?",FCT,mxGetClassName(ax),TSTR(T),e);

   len=mxGetNumberOfElements(a);
   rank=mxGetNumberOfDimensions(a); init_xstr();

   if (ax) { bool x=mxIsComplex(ax);
      if (cmplx<0) { cmplx=x; } else
      if (cmplx ^ x) { const char *s=xStr();
         if (cmplx) {
            if (WBLOG_TCAST) wblog(FL,
               "WRN %s unsetting cmplx->0 for %s%s",s,
               mxIsComplex(ax)? "complex ":"", mxGetClassName(ax));
            cmplx=0;
         }
         else wblog(FL,"ERR %s complex type mismatch (%d/%d)",s,cmplx,x);
      }
   }

   return *this;
};

template <class Ta, class Tb>
size_t copy_row2col_major(
   size_t r, const size_t *S, Ta *a, const Tb *b,
   char from 
){
   if (!r || !S) {
      wblog(FL,"WRN %s() r=%d / sz=%p !?",FCT,r,S);
      return 0;
   }
   if (!a || !b) wblog(FL,
      "ERR %s() got null pointers (%p/%p) !?",FCT,a,b);

   unsigned i=0, l=r-1;
   size_t j, k=0, len=S[0], I[r]; for (; i<r; ++i) { I[i]=0; }

   for (i=1; i<r; ++i) len*=S[i];

   for (k=0; k<len; ++k) {
       for (j=I[0], i=1; i<r; ++i) j = j*S[i] + I[i];

       if (from)
            { a[k]=Ta(b[j]); } 
       else { a[j]=Ta(b[k]); } 

       i=0; ++I[0]; 
       while (I[i]>=S[i] && i<l) { I[i]=0; ++I[++i]; }
   }

   return len;
};

template<class T>
template <class Tb>
int Mx::Array<T>::copy_from(const Tb* b, char tcheck) {

   if (!data || !ax) { if (data || b) wblog(FL,
      "ERR %s() got null data (%p, %p, %p)",FCT,ax,data,b);
      return 0;
   }

   if (!trans) {
      Wb::cpyRange(data,b,len,tcheck);
   }
   else {
      size_t l=
      copy_row2col_major(
         rank, mxGetDimensions(ax), data, b,
         'b' 
      );
      if (l!=len) wblog(FL,"ERR %s() %d/%d",FCT,l,len);
   }

   return (len ? 1 : 0);
};

template<class T>
size_t Mx::Array<T>::rcopy_to( 
   T* b, const wbperm *P, char ck) const { 

   size_t l=0;

   if (!ax) { if (data || b) wblog(FL,
      "ERR %s() got ax=%p (%p, %p)",FCT,ax,data,b);
      return l;
   }
   if (mxIsComplex(ax)) wblog(FL,
      "ERR %s got complex `%s'",xStr(),mxGetClassName(ax));

   if (data) { 
      return l=cpy_to_(b,P,ck,data);
   }

   mxClassID t=mxGetClassID(ax); 
   switch (t) {
      case mxDOUBLE_CLASS : l=cpy_to_(b,P,ck,mxGetDoubles (ax)); break;
      case mxCHAR_CLASS   : l=cpy_to_(b,P,ck,mxGetChars   (ax)); break;
      case mxINT8_CLASS   : l=cpy_to_(b,P,ck,mxGetInt8s   (ax)); break;
      case mxUINT8_CLASS  : l=cpy_to_(b,P,ck,mxGetUint8s  (ax)); break;
      case mxINT16_CLASS  : l=cpy_to_(b,P,ck,mxGetInt16s  (ax)); break;
      case mxUINT16_CLASS : l=cpy_to_(b,P,ck,mxGetUint16s (ax)); break;
      case mxINT32_CLASS  : l=cpy_to_(b,P,ck,mxGetInt32s  (ax)); break;
      case mxUINT32_CLASS : l=cpy_to_(b,P,ck,mxGetUint32s (ax)); break;
      case mxINT64_CLASS  : l=cpy_to_(b,P,ck,mxGetInt64s  (ax)); break;
      case mxUINT64_CLASS : l=cpy_to_(b,P,ck,mxGetUint64s (ax)); break;
      case mxLOGICAL_CLASS: l=cpy_to_(b,P,ck,mxGetLogicals(ax)); break;
      default: wblog(FL,"ERR %s unexpected type `%s'",xStr(),mxGetClassName(ax));
   }

   if (l!=len) wblog(FL,"ERR %s len=%d/%d !?",xStr(),l,len);
   return l;
};

namespace Mx { 

template<>
size_t Array<wbcomplex>::zcopy_to(
   wbcomplex *b, const wbperm *P, char ck) const { 

   size_t l=0;

   if (!ax) { if (data || b) wblog(FL,
      "ERR %s() got null data (%p -> %p, %p)",FCT,ax,data,b);
      return l;
   }

   if (data) { return l=cpy_to_(b,P,ck,data); }

   mxClassID t=mxGetClassID(ax); 
   if (t==mxDOUBLE_CLASS) {
      if (mxIsComplex(ax)) 
           l=cpy_to_(b,P,ck,(wbcomplex*)mxGetComplexDoubles(ax));
      else l=cpy_to_(b,P,ck,mxGetDoubles(ax));
   }
   else wblog(FL,"ERR %s() %s %s",FCT,
      mxIsComplex(ax) ? "complex":"type",mxGetClassName(ax));

   if (l!=len) wblog(FL,"ERR %s len=%d/%d !?",xStr(),l,len);
   return l;
}

template <>
size_t Mx::Array<wbcomplex>::copy_to(wbcomplex *b, const wbperm *P, char tcheck)
const { return zcopy_to(b,P,tcheck); };

};

template<class T>
template<class Tx>
size_t Mx::Array<T>::cpy_to_( 
   T* b, const wbperm *P, char ck, const Tx *xd) const {

   if (P && trans) wblog(FL,
      "ERR %s() got trans=%d together with permutation",FCT,trans);
   if (data && (void*)data!=(void*)xd) wblog(FL,
      "ERR %s got pointer mismatch %p / %p\n%s",xStr(),data,xd,STR_(this));

   if (!trans && !P) {
      Wb::cpyRange(b,xd,len,ck);
      return len;
   }

   if (!ax) wblog(FL,"ERR %s() %s",FCT,STR_(this));

   if (trans) {
      size_t l=
      copy_row2col_major(
         rank, mxGetDimensions(ax), b, xd,
         0 
      );
      if (l!=len) wblog(FL,"ERR %s() %ld/%ld",FCT,l,len);
   }
   else {
      size_t i=-1,j;
      unsigned k, l=rank-1;
      wbIndex I(rank);

      const size_t *s0=mxGetDimensions(ax), *p=P->data;

      if (!rank || P->len!=rank || len!=(i=mxGetNumberOfElements(ax)))
         wblog(FL,"ERR %s() size mismatch (r=%ld/%ld, len=%ld/%ld)",
         FCT, rank, P->len, len, i);

	  for (i=0; i<len; ++i) {
		  for (j=I[p[l]], k=l-1; k<l; --k) j = j*s0[p[k]] + I[p[k]];

          b[j]=T(xd[i]); 

          if (ck && Tx(b[j])!=xd[i]) wblog(FL,
          "ERR %s() rounding error (%g/%g)",FCT,double(b[j]),double(xd[i]));

		  k=0; ++I[0];
		  while (I[k]>=s0[k] && k<l) { I[k]=0; ++I[++k]; }
	  }
   }

   return len; 
};

bool Mx::IsIndex(const mxArray *a, int base) { 
   if (!a || !mxIsNumeric(a) || mxIsComplex(a)) { return 0; }

   const size_t *s=mxGetDimensions(a);
   int r=mxGetNumberOfDimensions(a);

   if (r!=2 || (s[0]!=1 && s[1]!=1)) { return 0; }

   const double *d=mxGetDoubles(a);
   if (d) { return IsIndex(d,s[0]*s[1],base); }
   else {
      wbvector<int> x(s[0]*s[1]);
      Mx::Array<int>(a).rcopy_to(x.data,0,'!');
      return IsIndex(x.data,x.len,base);
   }
};

template<class T>
wbstring Mx::Array<T>::toStr(char vflag, unsigned l) const {

   wbstring s(l); 

   if (ax) {
      unsigned i=0;
      const size_t *S=mxGetDimensions(ax);

      l=snprintf(s.data,s.len,"%s%s ",
         mxIsComplex(ax)? "complex ":"",mxGetClassName(ax));
      for (; i<rank && l<s.len; ++i) {
         l+=snprintf(s.data+l,s.len-l,"%ldx",S[i]); }

      if (l<s.len) { if (i) { --l; }
         l+=snprintf(s.data+l,s.len-l," [%s",TSTR(T)); }
      if (l<s.len) {
         l+=snprintf(s.data+l,s.len-l," @ ");
         if (!data) { l+=snprintf(s.data+l,s.len-l,"d=null; "); } else
         if (vflag) { l+=snprintf(s.data+l,s.len-l,"d=%p; ",data); }
      }
      if (l<s.len) {
         l+=snprintf(s.data+l,s.len-l,"%d,%d,%d; l=%ld]",
            cmplx,trans,mxref,len);
      }

      if (l>=s.len) wblog(FL,"ERR %s() string out of bounds "
         "(%d/%d)%N'%s'",FCT,l,s.len,s.data);
   }
   else {
      if (vflag) snprintf(s.data,s.len,
         "ax=%p -> %p (%s; %d,%d,%d; len=%ld)",ax,data,TSTR(T),
         cmplx,trans,mxref,len);
      else snprintf(s.data,s.len,
         "null%s (%s; %d,%d,%d; n=%ld)", data ? " -> d!=0":"",
         TSTR(T), cmplx,trans,mxref,len
      );
   }

   return s;
};

template<class T>
wbstring Mx::Array<T>::toStrT(char vflag) const {

   wbstring s(24); 
   unsigned l=0;

   if (!data || vflag) l=snprintf(s.data,s.len,"%s%s -> %s",
      mxIsComplex(ax)? "complex ":"",mxGetClassName(ax),TSTR(T));
   else l=snprintf(s.data,s.len,"%s",TSTR(T));

   if (l>=s.len) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)%N`%s'",FCT,l,s.len,s.data);
   return s;
};

template<class T>
void Mx::Array<T>::issue_type_error(
   const char *F, int L, unsigned ndim, const size_t *dims) {

   if (!ax || !mxGetNumberOfElements(ax)) {
      if (ndim) {
         size_t l=dims[0]; for (unsigned i=0; i<ndim; ++i) { l*=dims[i]; }
         if (l!=len) wblog(F_L,
            "ERR %s length mismatch (%ld/%ld)",xStr(),l,len);
         if (l) wblog(FL,
            "ERR %s failed to allocate memory\n%s (%ld)",xStr(),
            SSTR(wbvector<size_t>(ndim,dims,'r')), len
         );
      }
      if (len || data) wblog(F_L,"ERR %s got null (%d;%p)",xStr(),len,data);
   }
   else {
      wblog(F_L, "ERR %s type mismatch (%s%s/%s)",xStr(),
      mxIsComplex(ax) ? "complex ":"",mxGetClassName(ax),TSTR(T));
   }
};

template<class T> inline
int mxGetNumber(const mxArray *a, T &d, const char qflag) {

    int n=0; d=T(0);
    if (!a) { Wb::quietErrLog(FL,"mxArray is NULL",qflag); return 1; }

    try {
       Mx::Array<T> A(a); n=A.len;
       A.ncopy_to(&d,1,1); 
    }
    catch (...) { n=-1; }

    if (n==1) { return 0; }
    else {
       unsigned l,m=64; char msg[m];
       l=snprintf(msg,m,"%s %s() invalid input",SHORT_FL,FCT);
       if (n>1) { l+=snprintf(msg+l,m-l," (got array, n=%d)",n); } else
       if (!n ) { l+=snprintf(msg+l,m-l," (empty)"); }
       Wb::quietErrLog(FL,msg,qflag);
       return 1;
    }
};

inline int mxGetString(const mxArray *a, char *s) {
   if (mxGetString(a,s,127)) {
      sprintf_str("%s:%d ERR could not read string (%s) !?", FL, s);
      return 1;
   }
   return 0;
}

inline int mxGetString(const mxArray *a, wbstring &s) {
   char istr[128]; istr[0]=0;
   if (mxGetString(a,istr,127)) {
      sprintf_str("%s:%d ERR could not read string (%s) !?", FL, istr);
      return 1;
   }
   s=istr; return 0;
}

wbstring mxSize2Str(const mxArray *a) {
   if (!a) { return "(null)"; }

   size_t l=0, n=32; char s[n];
   int i=0, r=mxGetNumberOfDimensions(a);
   const size_t *S=mxGetDimensions(a);

   for (; i<r; ++i) {
       l+=snprintf(s+l,n-l,"%s%d", i ? "x":"", unsigned(S[i]));
       if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   }
   return s;
};

wbstring mxTypeSize2Str(const mxArray *a) {
   size_t n=64; char s[n];
   if (!a) { strcpy(s,"(null)"); return s; }

   size_t l=snprintf(s,n,"%s: ",mxGetClassName(a));
   int r=mxGetNumberOfDimensions(a);
   const size_t *S=mxGetDimensions(a);

   for (int i=0; i<r; ++i) {
      l+=snprintf(s+l,n-l,"%s%d",i?"x":"", (unsigned)S[i]);
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);
   }
   return s;
};

template<class T>
void mxGetVector(const char* F, int L,
   const mxArray *S, const char *vn, 
   wbvector<T> &v
){
   if (!mxIsStruct(S) || !Mx::IsScalar(S)) wblog(F_L,
      "ERR input must be scalar structure");

   mxArray *a=mxGetField(S,0,vn);
   if (!a) wblog(F_L,"ERR field %s does not exist",vn);

   mxGetVector(F_L,a,v);
};

double mxGetNumber( 
   const char *F, int L,
   const mxArray *S, const char *vn, const char qflag
){
   mxArray *a;
   wbcomplex z;

   if (!mxIsStruct(S)) wbdie(F_L,
      "need scalar structure on input");
   if (mxGetNumberOfElements(S)>1) wbdie(F_L,
      "need SINGLE structure on input");

   a=mxGetField(S,0,vn);

   if (!a) wblog(F_L,
      "ERR field %s does not exist in structure", vn);

   if (mxGetNumber(a,z,qflag)) wbdie(FL,str);
   if (z.i!=0.) wblog(FL,"WRN Input number is complex! (%g)", z.i);

   return z.r;
}

template<class T>
mxArray* cpyRange2Mx(const T* d, const wbvector<size_t> &S0){

   mxArray *a;
   unsigned i=0, n=S0.prod(0), r=S0.len, len=MAX(2U,r);
   size_t s[len]; 

   for (; i<r; ++i) s[i]=S0[i];
   if (i<len) {
      if (n) for (; i<len; ++i) s[i]=1; 
      else   for (; i<len; ++i) s[i]=0; 
   }

   a=Mx::Array<double>(len,s).copyFromTR(d);
   if (!a && d) wblog(FL,
     "ERR %s() got %p -> %p @ len=%d !?",FCT,d,a,len);

   return a;
};

template<>
mxArray* cpyRange2Mx(const char* d0, const wbvector<size_t> &S0){

   unsigned i, r=S0.len, n=S0.prod(0), len=MAX(2U,r);
   mxArray *a;

   size_t s[len];
   for (i=0; i<r; i++) s[i]=S0[i];
   if (i<2) { if (i) s[1]=1; else s[0]=s[1]=0; }

   a=Mx::Array<mxChar>(len,s).Return();

   if (a==NULL) {
      if (!d0) wblog(FL,
         "ERR %s() invalid usage (%lX, %lX)",__FUNCTION__,d0,a);
      else return a;
   }
   if (!n) return a;

   mxChar *d = (mxChar*)mxGetData(a);
   if (!d) wblog(FL,"ERR invalid mxCharArray");

   for (i=0; i<n; ++i) { d[i]=mxChar(d0[i]);
      if (char(d[i])!=d0[i]) wblog(FL,
         "ERR %s() got data conversion error (%d/%d)",
         FCT,int(d[i]),int(d0[i])
      );
   }

   return a;
};

template<>
mxArray* cpyRange2Mx(const wbcomplex* z, const wbvector<size_t> &S){

   size_t i=0, n=S.prod(0); int isz=0; 
   mxArray *a;

   for (; i<n; ++i) { if (z[i].i!=0) { isz=1; break; }}

   if (isz) {
      a=Mx::Array<wbcomplex>(S).copyFromTR(z);
   }
   else {
      Mx::Array<double> A(S); if (A.data) {
      for (i=0; i<n; ++i) { A.data[i]=z[i].r; }}
      a=A.Return(); 
   }

   if (!a && z) wblog(FL,"ERR %s() failed to allocate data !?",FCT);
   return a;
};

mxArray* matGetVariable(
    const char *F, int L,
    const char *fname, const char *vn, char force
){
    Wb::matFile f;
    int i=f.open(F,L,fname,"r", force); if (!i) return NULL;
    mxArray *a=matGetVariable(f.mfp,vn);
    return a;
};

mxArray* matGetVariableInfo(
    const char *F, int L,
    const char *fname, const char *vn, char force
){
    Wb::matFile f;
    int i=f.open(F,L,fname,"r", force); if (!i) return NULL;
    mxArray *a=matGetVariableInfo(f.mfp,vn);
    return a;
};

double* mexSetCell2Matrix(mxArray *C, int fid, unsigned d1, unsigned d2) {

    double *data;
    Mx::Array<double> A(d1,d2); data=A.data;

    if (A.data) { mxSetCell(C,fid,A.Return()); } else
    if (A.ax || (d1 && d2)) { wblog(FL,"ERR %s() %s",FCT,STR(A)); }

    return data;
};

namespace Wb {

template <class TQ, class TD>
mxArray* mxCreateSparse(
   const char *F, int L, unsigned d1, unsigned d2,
   const wbMatrix<TQ> &IJ, 
   const wbvector<TD> &D,  
   const wperm_t *p,       
   int *info
){
   const char lex=-1;
   widx_t n=D.len;
   bool gotp=(p!=NULL);

   if (D.data && ISCOMPLX_(TD)) wblog(F,L,
      "ERR %s() got complex sparse data set (%dx%d; %d)",FCT,d1,d2,D.len);
   if (IJ.dim1!=n || (n && IJ.dim2!=2)) wblog(F,L,
      "ERR %s() size mismatch (%dx%d : ij=%dx%d @ %d)",
       FCT,d1,d2,IJ.dim1,IJ.dim2,D.len);

   if (!gotp && !IJ.isSorted(+1,lex)) {
      wbMatrix<TQ> IX(IJ);
      wbperm P;

      IX.SortRecs(P,+1,lex);

      return mxCreateSparse(F_L,d1,d2,IX,D,P.data);
   }

   mxArray *a=mxCreateSparse(d1,d2,D.len,mxREAL);
   if (!a) wblog(F,L,"ERR %s() "
      "failed to allocate sparse matrix (%dx%d; %d)",FCT,d1,d2,D.len);
   if (!n) { return a; }

   double *sr=mxGetDoubles(a);
   mwIndex *is=mxGetIr(a), *js=mxGetJc(a),k=0;
   unsigned i,j=0, nz=0, m=0; TQ *ij=IJ.data;

   if (!WbUtil<TQ>::isInt()) {
      wblog(FL,"ERR %s() invalid non-integer data type '%s'",FCT,
      TSTR(TQ));
   }

   for (i=0; i<n; ++i, ij+=2) {
      const TD& x=(gotp ? D.data[p[i]] : D.data[i]);
      if (int(ij[0])<0 || ij[0]>=d1 || int(ij[1])<0 || ij[1]>=d2) wblog(FL,
         "ERR %s() index out of bounds (%d,%d) [%dx%d]",
         FCT,ij[0],ij[1],d1,d2);
      if (x!=0) {
         if (i) {
            if (ij[-2]==ij[0] && ij[-1]==ij[1]) {
               sr[k]+=double(x); ++m; continue; 
            }
            else if (ij[-1]>ij[1] || (ij[-1]==ij[1] && ij[-2]>ij[0])) {
               wblog(FL,"ERR %s() data not sorted (%d/%d) !?",FCT,i+1,n);
            }
            else ++k;
         }
         sr[k]=double(x); 
         is[k]=ij[0];

         for (; j<=ij[1]; ++j) js[j]=k;
      }
      else if (ij[0]!=ij[1]) ++nz; 
   }

   for (++k; j<=d2; ++j) js[j]=k;

   if (F) {
      if (nz) { wblog(F,L, 
         "WRN %s() got off-diag zero data (%d/%d) !?",FCT,nz,D.len);
      }
      if (m) { wblog(F,L,
         "WRN %s() got non-unique sparse data (%d/%d)",FCT,m,D.len);
      }
   }
   if (info) (*info)=nz+m;

   return a;
};

} 

int mxAddField2Scalar(
   const char *F, int L,
   mxArray* S, const char *vn, mxArray *a
){
   int fid;

   if (!mxIsStruct(S)) wblog(F_L,
      "ERR scalar input structure required");
   if (mxGetNumberOfElements(S)>1) wblog(F_L,
      "ERR single input structure required");

   fid=mxAddField(S,vn);

   if (fid<0) wblog(F_L,
      "ERR failed to add field '%s' !? (%d)", vn,fid);

   if (a) mxSetFieldByNumber(S,0,fid,a);
   return fid;
};

int mxAddField(
    const char *F, int L,
    mxArray* S, const char *name
){
    int fid = mxAddField(S, name);

    if (fid<0) wblog(F,L,
    "ERR Failed to add field %s to structure (%d) ???", name, fid);

    return fid;
}

void mxAppendStructToStruct(
   const char *F, int L,
   mxArray *&S0, mxArray *S 
){
   unsigned i,n;
   const char* name;

   if (!mxIsStruct(S0) || mxGetNumberOfElements(S0)!=1 ||
       !mxIsStruct(S ) || mxGetNumberOfElements(S )!=1) wblog(F,L,
      "ERR %s requires two scalar input structures",__FUNCTION__);

   n=mxGetNumberOfFields(S0);

   for (i=0; i<n; i++) {
      name=mxGetFieldNameByNumber(S0,i);

      if (mxGetFieldNumber(S,name)>=0) wblog(F,L,
         "WRN overwriting field `%s' in structure", name);

      mxAddField2Scalar(F,L,
        S, name, mxGetFieldByNumber(S0,0,i)
      );

      mxSetFieldByNumber(S0,0,i,NULL);
   }

   mxDestroyArray(S0); S0=NULL;
}

int mxUpdateField(
    const char *F, int L,
    mxArray* S, const char *name, unsigned k, mxArray *a
){
    int fid=mxGetFieldNumber(S,name);
    unsigned n;

    if (fid>=0) {
       mxArray *x=mxGetFieldByNumber(S,k,fid);
       if (x) mxDestroyArray(x);
    }
    else fid=mxAddField(F,L, S, name);

    n=(unsigned)mxGetNumberOfElements(S);

    if (k>=n) wblog(FL,
    "ERR %s: index out of bounds (%d/%d)",FCT,k,n);

    mxSetFieldByNumber(S,k,fid,a);

    return fid;
}

void matPutVariable(const char *F, int L,
   const char *file, const char *name, mxArray *a, unsigned keep
){
   Wb::matFile f;

   f.Open(F,L,file,"u"); 

   if (matPutVariable(f.mfp, name, a)) wblog(F,L,
      "ERR Failed to write parameter`%s' to file `%s'", name, file);

   f.close();
   if (!keep) mxDestroyArray(a);
}

void mxReplaceField(
    const char *F, int L,
    mxArray* S, unsigned k, int fid, mxArray *a 
){
    if (!S || !mxIsStruct(S)) wblog(F_L,
       "ERR %s() got %s",FCT, S ? mxGetClassName(S):"(empty)");

    mxArray *x=mxGetFieldByNumber(S,k,fid);
    if (x) mxDestroyArray(x);

    mxSetFieldByNumber(S,k,fid,a); 
};

void mxSetFieldToNumber(
    mxArray *S, unsigned i, unsigned fid, wbcomplex z, char check) {

    if (check) { unsigned n;
       if (!mxIsStruct(S)) wblog(FL,"ERR invalid input structure");
       if (i>=(n=mxGetNumberOfElements(S))) wblog(FL,
          "ERR index out of bounds (i=%d/%d)",i,n);
       if (fid>=(n=mxGetNumberOfFields(S))) wblog(FL,
          "ERR index out of bounds (fid=%d/%d)",fid,n);
    }

    mxArray *a; size_t sz[]={1,1};
    if (z.i)
         { a=Mx::Array<wbcomplex>(2,sz).copyFromTR(&z  ); }
    else { a=Mx::Array<double   >(2,sz).copyFromTR(&z.r); }

    mxSetFieldByNumber(S,i,fid,a);
};

int isFlagGlobal(const char *vn) {
    double flag=0;
    if (getNumGlobal(vn,flag)) return 0;
    return (flag!=0);
}

template<class T>
mxArray* vecVec2Mx(
    const wbvector< wbvector<T> > &V,
    const char tflag
){
    mxArray *C=mxCreateCellMatrix(1,V.len);

    for (unsigned i=0; i<V.len; i++)
    mxSetCell(C, i, V[i].toMx(tflag));

    return C;
}

template<class T>
mxArray* vecMat2Mx(
    const wbvector< wbMatrix<T> > &M,
    const char rawflag
){
    mxArray *C=mxCreateCellMatrix(1,M.len);

    for (unsigned i=0; i<M.len; i++)
    mxSetCell(C, i, M[i].toMx(rawflag));

    return C;
}

template<class T>
mxArray* vecVec2Mx(
    const wbvector< wbvector<T>* > &V,
    const char tflag
){
    mxArray *C=mxCreateCellMatrix(1,V.len);
    wbvector<T> emptyVec;

    for (unsigned i=0; i<V.len; i++)
    if (V[i]) mxSetCell(C,i,V[i]->toMx(tflag));
    else mxSetCell(C, i, emptyVec.toMx(tflag));

    return C;
}

#endif

