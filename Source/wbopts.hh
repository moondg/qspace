/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace MEX routine options class
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

#ifndef __WB_OPTS_HCC__
#define __WB_OPTS_HCC__

/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* class to manage optional input arguments */

class OPTS {

  public:

    OPTS () : mto(0) {};

    OPTS (const mxArray** ain, const unsigned n)
    : mto(0) { init(ain,n); };

    OPTS (const char *fname)
    : mto(0) { init(fname); };

    void init (const mxArray** ain, const unsigned n);
    void init (const char *fname);

   ~OPTS () { init(); }
    void init() {
       name.init();

       if (MAT.mfp && aa.len) 
       for (unsigned i=0; i<aa.len; ++i) mxDestroyArray(aa[i]);

       aa.init();
    };

    int len() { return ( MAT.mfp ? -1 : (int)aa.len ); };

    template <class T>
    bool getOpt(const char*, int,
         const char* vn, T &x, char force=0, const char *istr=0);

    template <class T>
    bool getOpt(const char*, int,
         const char* vn, wbvector<T> &x, char force=0, const char *istr=0);

    template <class T>
    bool getOpt(const char*, int,
        const char* vn, wbMatrix<T> &x, char force=0, const char *istr=0);

    bool getOpt(const char*, int, const char* vn, const char *istr=0);

    template <class T>
    bool getOpt(const char* vn, T &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0, vn,x,force,istr);
    };

    template <class T>
    bool getOpt(const char* vn,
       wbvector<T> &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0,vn,x,force,istr);
    };

    template <class T>
    bool getOpt(const char* vn,
       wbMatrix<T> &x, char force=0, const char *istr=0) {
       return getOpt(NULL,0,vn,x,force,istr);
    };

    bool getOpt(const char* vn, const char *istr=0) {
       return getOpt(NULL,0,vn,istr);
    };

    void checkAnyLeft(const char *MAT=NULL, int L=0);

  protected:
  private:

    wbvector<mxArray*> aa;
    wbvector<wbstring> name;

    unsigned mto; 

    Wb::matFile MAT;

    int findOpt(const char* vn);
};

int OPTS::findOpt(const char* vn) {

   unsigned i, n=name.len;
   int k=-1;

   for (i=0; i<n; ++i)
   if (name[i].len) if (!strcmp(name[i].data,vn)) {
      if (k<0) k=i; else ++mto;
   };

   return k;
};

void OPTS::init (const char *fname) {

   if (aa.len) wblog(FL,
      "ERR OPTS already set by mxArray[] (n=%d)\ngot %s",aa.len,fname);

   MAT.Open(FL, fname, "r"); 
};

bool OPTS::getOpt(
   const char *F, int L, const char *vn, const char *istr
){
   unsigned r=0;

   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, vn); r=(a!=NULL);
      if (a) {
         double dbl=0.;
         if (mxGetNumber(a,dbl)) { wblog(FL,"ERR %s",str); }
         if (dbl==0) r=0;
      }
      mxDestroyArray(a);
   }
   else {
      int i=findOpt(vn); r=(i>=0);
      if (r) name[i].init(); 
   }

   if (r && F) {
      if (istr && istr[0])
           { wblog(F,L," *  %s",istr); } 
      else { wblog(F,L," *  %s",vn); }
   }

   return r;
};

template <class T>
bool OPTS::getOpt(
   const char *F, int L, const char* vn, T &x0,
   char force, const char *istr
){
   unsigned k,e=0;
   mxArray *a; T x=0;

   if (MAT.mfp) {
      a = matGetVariable(MAT.mfp, vn);
      if (!a) {
         if (force) wblog(F_L,
            "ERR %s() option '%s' required (%s)",FCT,vn, istr?istr:"");
         return 0;
      }
   }
   else {
      int i=findOpt(vn);
      if (i<0) {
         if (force) wblog(F_L,
            "ERR %s() option '%s' required (%s)",FCT,vn, istr?istr:"");
         return 0;
      }

      k=(unsigned)(i+1);

      if (k>=name.len) { ++e; } else
      if (aa[k]==NULL) { ++e; }
      if (e) wblog(F_L,"ERR %s() missing value for '%s'",FCT,vn);

      a=aa[k];

      name[k-1].init(); aa[k]=NULL;
   }

   if (mxGetNumber(a,x)) wblog(F_L, 
      "ERR %s() failed to read value for '%s'\n%s",FCT,vn,str);

   if (F && x0!=x) {
      if (istr && istr[0])
           { wblog(F,L," *  %-8s: %s [%s]", vn, NSTR(x), istr); }
      else { wblog(F,L," *  %-8s: %s",      vn, NSTR(x)); }
   }

   x0=x;

   if (MAT.mfp) mxDestroyArray(a);

   return 1;
};

template<>
bool OPTS::getOpt(const char *F, int L,
   const char* vn, 
   wbstring &x, char force, const char *istr
){
   int lflag=0;
   if (MAT.mfp) {
      wbvec<char> s(128);
      mxArray *a = matGetVariable(MAT.mfp, vn);
      if (!a) {
         if (force) wblog(F_L,"ERR %s() missing input for '%s'",FCT,vn);
         return 0;
      }

      if (!mxIsChar(a)) wblog(F_L,
         "ERR %s() invalid type `%s' for '%s'",FCT,mxGetClassName(a),vn);
      if (mxGetString(a,s.data,s.len)) wblog(F_L,
         "ERR %s() failed to read string for '%s' (len<%d?)",FCT,vn,s.len);

      if (x.data && !strcmp(x.data,s.data)) { ++lflag; }
      x=s.data; mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(vn);
      if (i<0) {
         if (force) wblog(F_L,"ERR %s() option '%s' required",FCT,vn);
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) { ++e; } else
      if (aa[k]!=NULL) { ++e; }
      if (e) wblog(F_L,"ERR %s() missing value for option '%s'",FCT,vn);

      if (x.data && !strcmp(x.data,name[k].data)) ++lflag;
      x=name[k];

      name[k-1].init(); name[k].init();
   }

   if (F && lflag) {
      if (istr && istr[0])
           { wblog(F,L," *  %-8s: %s [%s]",vn, x.data, istr); }
      else { wblog(F,L," *  %-8s: %s",     vn, x.data); }
   }

   return 1;
};

template<>
bool OPTS::getOpt(const char *F, int L,
   const char* vn, mxArray* &a, char force, const char *istr
){
   if (MAT.mfp) {
      a=matGetVariable(MAT.mfp, vn);
      if (a) aa.Append((mxArray*)a); 
   }
   else {
      int i=findOpt(vn);

      if (i<0) a=NULL;
      else {
         unsigned e=0, k=unsigned(i+1);

         if (k>=name.len) { ++e; } else
         if (aa[k]==NULL) { ++e; }
         if (e) wblog(F_L,"ERR %s() value missing for '%s'",FCT,vn);

         a=aa[k];

         name[k-1].init(); aa[k]=NULL;
      }
   }

   if (!a) {
      if (force) wblog(F_L,"ERR %s() missing value for '%s'",FCT,vn);
      return 0;
   }

   if (F) {
      if (istr && istr[0])
           { wblog(F,L," *  %-8s: (%s) [%s]",vn, mxGetClassName(a), istr); }
      else { wblog(F,L," *  %-8s: (%s)",     vn, mxGetClassName(a)); }
   }

   return 1;
};

template<>
bool OPTS::getOpt(
   const char *F  QS_UNUSED_VAR,
   int L          QS_UNUSED_VAR,
   const char* vn, itag_ &t, char force, const char *istr 
){
   wbstring q;
   getOpt(vn,q,force,istr);
   if (q) { t.init(0,0,q.data); return 1; }
   return 0;
};

template<class T>
bool OPTS::getOpt(const char *F, int L,
   const char* vn, wbvector<T>&x_, char force, const char *istr
){
   int lflag=0;
   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, vn);
      if (!a) {
         if (force) wblog(F_L,"ERR %s() missing input for '%s'",FCT,vn);
         return 0;
      }

      if (!mxIsDouble(a)) wblog(F_L,
         "ERR %s() invalid type `%s' for '%s'",FCT,mxGetClassName(a),vn);

      wbvector<T> x; x.init(FL,a);
      if (x!=x_) { x.save2(x_); ++lflag; }
      mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(vn);
      if (i<0) {
         if (force) { wblog(F_L,"ERR %s() missing input %s",FCT,vn); }
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) { ++e; } else
      if (aa[k]==NULL) { ++e; }
      if (e) wblog(F_L,"ERR %s() missing value for '%s' (vector)",FCT,vn);

      wbvector<T> x; x.init(FL,aa[k]);
      if (x!=x_) { x.save2(x_); ++lflag; }

      name[k-1].init(); aa[k]=NULL;
   }

   if (F && lflag) {
      if (istr && istr[0])
           { wblog(F,L," *  %-8s: [%s] (%s)", vn, x_.toStr().data, istr); }
      else { wblog(F,L," *  %-8s: [%s]",      vn, x_.toStr().data); }
   }

   return 1;
}

template<class T>
bool OPTS::getOpt(const char *F, int L,
   const char* vn, wbMatrix<T>&x_, char force, const char *istr
){
   int lflag=0;
   if (MAT.mfp) {
      mxArray *a = matGetVariable(MAT.mfp, vn);
      if (!a) {
         if (force) wblog(F_L,
            "ERR %s() missing value for %s (matrix)",FCT,vn);
         return 0;
      }

      if (!mxIsDouble(a)) wblog(F_L,
         "ERR %s() invalid type `%s' for '%s'",FCT,mxGetClassName(a),vn);

      wbMatrix<T> x; x.init(FL,a);
      if (x!=x_) { x.save2(x_); ++lflag; }
      mxDestroyArray(a);
   }
   else {
      int i; unsigned k,e=0;

      i=findOpt(vn);
      if (i<0) {
         if (force) wblog(F_L,
            "ERR %s() missing input %s (matrix)",FCT,vn);
         return 0;
      }
      k=unsigned(i+1);

      if (k>=name.len) { ++e; } else
      if (aa[k]==NULL) { ++e; }
      if (e) wblog(F_L,
         "ERR %s() missing input data for '%s' (matrix)",FCT,vn);

      wbMatrix<T> x; x.init(FL,aa[k]);
      if (x!=x_) { x.save2(x_); ++lflag; }

      name[k-1].init(); aa[k]=NULL;
   }

   if (F && lflag) {
      if (istr && istr[0])
           { wblog(F,L," *  %-8s: [%s] (%s)", vn, x_.toStr().data, istr); }
      else { wblog(F,L," *  %-8s: [%s]",      vn, x_.toStr().data); }
   }

   return 1;
}

void OPTS::checkAnyLeft(const char *F, int L) {

   unsigned i,e=0;

   if (MAT.mfp) {
      MAT.init();
      for (i=0; i<aa.len; ++i) mxDestroyArray(aa[i]);
      return;
   }

   for (i=0; i<name.len; ++i) { if (name[i].len || aa[i]) {
       if ((++e)==1) printf("\n");

       if (name[i].len) {
          printf("%5d: '%s'\n", i+1, name[i].data);
       }
       else if (aa[i]) {
          int m=mxGetM(aa[i]), n=mxGetN(aa[i]);
          if (mxIsCell(aa[i]))
               printf("%5d: {%dx%d %s}\n",i+1,m,n,mxGetClassName(aa[i]));
          else printf("%5d: [%dx%d %s]\n",i+1,m,n,mxGetClassName(aa[i]));
       }
   }}

   if (e) { PRINTF("\n"); if (mto) {
      wblog(F_L,"WRN some options appear more than once (m=%d)",mto); }
      wblog(F_L,"ERR unused/unrecognized options");
   }
};

#endif

