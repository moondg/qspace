/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace memory routines
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

#ifndef __WB_MEM_LIB_CC__
#define __WB_MEM_LIB_CC__

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //

template <class T>
wbstring Wb::sptr<T>::toStr() const {

   wbvec<char> s(128);

   if (mtype!=Wb::MEM_DEF) { s.catf(FL,"%s ",Wb::MTYPE_STR[mtype]); }

   if (data || len)
        { s.catf(FL,"%p -> %p: %s[%3g]",this,data,sTSTR(T),double(len)); }
   else { s.catf(FL,"%p -> %p",this,data); }

   if (nref || mtype) s.catf(FL,
     " @ %d ref%s / %s",nref, nref>1?"s":"", MTYPE_STR[mtype]);

   return s.data;
};

template <class T> inline
T* Wb::sptr<T>::malloc_base(size_t n, char init, Wb::MTYPE mt) {

   char do_alloc=(n>0);
   if (data) {
      if (check_consistency(0)) wblog(FL,"ERR Wb::sptr() %s",STR(*this));
      if (n!=len || nref) rm_dref();
      else {
         do_alloc=0;
      }
   }

   if (do_alloc) {
      if (data || len) wblog(FL,"ERR %s() %s !?",FCT,STR(*this));
      len=n; mtype=mt;

      if (mtype==Wb::MEM_DEF) {
         WB_NEW(data,n, init ? 0 : 1); 
         init=0; 
      }
      else if (mtype==Wb::MEX_RETURN) {
         if (typeid(T)==typeid(double)) {
            mxArray *a=Mx::Array<double>(1,&n).Return();
            init=0; 

            data=(T*)mxGetDoubles(a); 
            gP2X.add(data,a); 

            nref=0; 
         }
         else if (typeid(T)==typeid(wbcomplex)) {
            mxArray *a=Mx::Array<wbcomplex>(1,&n).Return();
            init=0; 

            data=(T*)mxGetComplexDoubles(a); 
            gP2X.add(data,a); 
         }
         else wblog(FL,"ERR %s() got %s",FCT,TSTR(T));

         if (WBLOG_MMEX) wblog(FL,"new %s",STR(*this));
      }
   }
   if (n && init) { T z=T(); 
      for (size_t i=0; i<n; ++i) { data[i]=z; }
   }

   if ((!len ^ !data) || (nref && !len)) wblog(FL,
      "ERR Wb::sptr<%s>=%p -> %s",sTSTR(T),this,STR(*this));

   return data;
};

template <class T> inline
T* Wb::sptr<T>::malloc_data(size_t n, const T *d0, Wb::MTYPE mt) {

   malloc_base(n,0,mt); 
   if (n) {
      if (d0)
           MEM_CPY<T>(data,n,d0);
      else MEM_SET<T>(data,n);
   }

   return data;
};

template <class T> inline
T* Wb::sptr<T>::init2ref(size_t n, const T* d0) {

   if (data==d0 && len==n && nref) { return data; }
   if (data || len || nref) wblog(FL,"ERR %s() got %s",FCT,STR(*this));
   if (!d0 ^ !n) wblog(FL,"ERR %s() got n=%ld with %p",FCT,n,d0);

  #pragma omp critical (manage_SPTR) 
   {  data=(T*)d0;  
      len=n; mtype=Wb::MEM_REF;
      nref=1;
   }
   return data;
};

template <class T> inline
unsigned Wb::sptr<T>::add_dref() {

   if (!data) { check_consistency(); return 0; }

  #pragma omp critical (manage_SPTR) 
   { ++nref; }

   if (nref>128) wblog(FL,"WRN %s() %p @ nref=%d",FCT,data,nref);

   return nref;
};

template <class T> inline
int Wb::sptr<T>::rm_dref(char mflag) {

   int done=0, e=0;

   #pragma omp critical (manage_SPTR)
   { try {

      if (!len ^ !data) { e|=1; }  
      if (nref && !len) { e|=2; }

      if (!e && data) {
         if (mtype==Wb::MEM_DEF) {
            if (nref) { --nref; } 
            else {
               WB_DELETE(data); 
               done=1; 
            }
         }
         else if (mtype==Wb::MEM_REF) {
            if (nref==1) { done=2; } 
            else if (nref>1) { --nref; }
            else { e|=16; }
         }
         else if (mtype==Wb::MEX_RETURN) {
            if (!nref) { done=3; } 
            else { e|=32; }

            if (mflag || (wbl::status & WBL_ERR__)) { 
               mxArray *a=gP2X.Return(0,0,data);
               if (a) { mxDestroyArray(a); data=NULL; }
            }
            else if (gP2X.erase(data)) { e|=64; } 
         }
         else e|=128;
      }
   } catch (...) { e|=256; }}

   if (e) {
      if (e==64)
           { wblog(FL,"WRN %s() got %s",FCT,STR(*this)); }
      else { wblog(FL,"ERR %s() got e=%d\n%s",FCT,e,STR(*this)); }
   }

   if (done) { init_def(); } 
   else if (!data) { mtype=Wb::MEM_DEF; } 

   return done;
};

template <class T> inline
T* Wb::sptr<T>::save_dref(const char *F, int L) {

  T* dref=data; int e=0;

  #pragma omp critical (manage_SPTR)
   { if (!len ^ !data) { e|=1; }  
     if (nref || mtype) { e|=2; }
     else if (len) { len=0; data=NULL; }
   }

   if (e) wblog(F_L,"ERR %s() %s got e=%d !?",FCT,STR(*this),e);
   return dref;
};

template <class T1, class T2>
bool gotMemOverlap(
   const T1* src, size_t len0,
   const T2* dest, size_t len2
){
   if (len0 && len2) {
      const char
        *L=(char*)src,  *R=(char*)(src+len0)-1,
        *a=(char*)dest, *b=(char*)(dest+(long(len2)<0 ? len0 : len2))-1;

      if (!src || !dest) wblog(FL,
         "ERR %s() received null pointer (%p, %p)",FCT,src,dest);
      return ((a<=L && b>=L) || (a<=R && b>=R) || (a>=L && b<=R) );
   }
   return 0;
};

#endif

