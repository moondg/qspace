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

#ifndef __WB_MEXLIB_HH__
#define __WB_MEXLIB_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

   void mexDisp(const mxArray *a, const char *vname=0);

   void mexWRN(const char* s); 

   wbstring sprint_info(const mxArray *a); 

   double wbtoc(
      const char *F=NULL, int L=0,
      const char *istr="", char restart=0
   );

   void wbtic() { wbtoc(0,0); }

   template<class T>
   int mxGetNumber(const mxArray *a, T &d, const char qflag=0);

   double mxGetNumber(
      const char *F, int L,
      const mxArray *S, const char *vname, const char qflag=0
   );

   int getFlag(const mxArray *a) {
      double q=0;
      Mx::Array<double> A(a);
      if (A.len>1) wblog(FL,
         "ERR %s() invalid flag (got array, n=%d)",FCT,A.len);
      A.ncopy_to(&q,1); 
      return int(q);
   };

   int mxGetString(const mxArray *a, char *s);
   int mxGetString(const mxArray *a, wbstring &s);

   wbstring mxTypeSize2Str(const mxArray *a);
   wbstring mxSize2Str(const mxArray *a);

   template<class T>
   void mxGetVector(const char* F, int L, const mxArray *a, wbvector<T> &v) {
      size_t n=mxGetNumberOfElements(a);
      v.init(n);
      try {
         Mx::Array<T>(a).copyTo(v.data);
      }
      catch (...) { wblog(F_L,"ERR %s()",FCT); }
   };

   template<class T>
   void mxGetVector( const char* F, int L,
      const mxArray *S, const char *vname, 
      wbvector<T> &v);

   mxArray* matGetVariable(const char *F, int L,
      const char *fname, const char *vname, char force=0
   );
   mxArray* matGetVariableInfo(const char *F, int L,
      const char *fname, const char *vname, char force=0
   );

   double* mexSetCell2Matrix(mxArray *a, int fid, unsigned d1, unsigned d2);
   void mxSetFieldToNumber(mxArray *S, unsigned i, unsigned fid,
        wbcomplex z, char dblcheck=1);

   template<class T>
   int getNumGlobal(const char *vname, T &x, const char *ws="global");
   int isFlagGlobal(const char *vname);

   template<class T>
   mxArray* vecVec2Mx(
   const wbvector< wbvector<T> > &V, const char tflag=0);

   template<class T>
   mxArray* vecVec2Mx(
   const wbvector< wbvector<T>* > &V, const char tflag=0);

   template<class T>
   mxArray* vecMat2Mx(
   const wbvector< wbMatrix<T> > &M, const char rawflag=0);

   template<class T>
   mxArray* cpyRange2Mx(const T* d, const wbvector<size_t> &S0);

   int mxAddField2Scalar(
       const char *F, int L,
       mxArray* S, const char *vname, mxArray *a=NULL
   );

   void mxAppendStructToStruct(
      const char *F, int L,
      mxArray *&S0, mxArray *S 
   );

   int mxAddField(
       const char *F, int L,
       mxArray* S, const char *vname
   );

   void mxReplaceField( 
       const char *F, int L,
       mxArray* S, unsigned k, int fid, mxArray *a
   );

   int mxUpdateField(   
       const char *F, int L,
       mxArray* S, const char *name, unsigned k, mxArray *a
   );

   void matPutVariable(const char *F, int L,
      const char *file, const char *name, mxArray *a, unsigned keep=0
   );

namespace Wb {

   template <class TQ, class TD>
   mxArray* mxCreateSparse(
      const char *F, int L, unsigned d1, unsigned d2,
      const wbMatrix<TQ> &IJ, const wbvector<TD> &D,
      const wperm_t *p=NULL, int *info=NULL
   );

}

template<class T>
mxArray* numtoMx(const T &x) { 
    mxArray *a=mxCreateDoubleMatrix(1,1,mxREAL);
    mxGetDoubles(a)[0]=double(x); return a;
};

template<>
mxArray* numtoMx(const wbcomplex &z) {
    mxArray *a;
    if (z.i!=0.) {
       a=mxCreateDoubleMatrix(1,1,mxCOMPLEX);
       z.explicit_copy2(mxGetComplexDoubles(a)[0]);
    }
    else a=numtoMx(z.r);
    return a;
};

template<>
mxArray* numtoMx(const char &c) {
    if (isalpha(c)) {
       char s[2]={c,0};
       return mxCreateString(s);
    }
    else return numtoMx(double(c));
};

class mxStruct {
 public:

    mxStruct(const char *file, int line) : S(NULL) { init(file,line); };
   ~mxStruct() {
       if (S) {
          wblog(FL,"WRN mxStruct() destroyed without further ado");
          mxDestroyArray(S);
       }
   };

    void init(const char *file, int line) {
       F=file; L=line;
       if (S) mxDestroyArray(S);
       S=mxCreateStructMatrix(1,1,0,NULL);
    };

    void put(const char *name) {
       if (S) { mxPutAndDestroy(F.data,L,S,name); S=NULL; } else
       wblog(FL,"WRN mxStruct::put() not initialized yet - return");
       F.init(); L=0;
    };

    void save2Struct(mxArray *S0, const char *name) {
       if (S) {
          mxAddField2Scalar(F.data, L, S0, name, S);
          S=NULL;
       } else
       wblog(FL,"WRN mxStruct::save2Struct() not initialized yet");
       F.init(); L=0;
    };

    void checkNewField(const char *f=NULL, int l=0);

    mxStruct& operator()(const char *fn) { fld=fn; return *this; };

    template<class T>
    mxStruct& operator=(const T&x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, numtoMx(x));
       return *this;
    };

    mxStruct& operator=(mxArray *a) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, a);
       return *this;
    };

    template<class T>
    mxStruct& operator=(const wbvector<T>& x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, x.toMx());
       return *this;
    };

    template<class T>
    mxStruct& operator=(const wbMatrix<T>& x) {
       checkNewField(FL);
       mxAddField2Scalar(F.data, L, S, fld.data, x.toMx());
       return *this;
    };

    mxArray *S;
    wbstring F, fld;
    unsigned L;

 protected:
 private:

};

template<>
mxStruct& mxStruct::operator=(const wbperm &P) {
   checkNewField(FL);
   mxAddField2Scalar(F.data, L, S, fld.data, P.toMx());
   return *this;
}

template<>
mxStruct& mxStruct::operator=(const wbindex &I) {
   checkNewField(FL);
   mxAddField2Scalar(F.data, L, S, fld.data, I.toMx());
   return *this;
}

inline void mxStruct::checkNewField(const char *f, int l) {
   int i; if (!f) { f=F.data; l=L; }

   if (!S) wblog(f,l,"ERR mxStruct() S not set yet");
   if (fld.isEmpty()) wblog(f,l,"ERR mxStruct() field not specified");

   i=mxGetFieldNumber(S,fld.data);
   if (i>=0) wblog(f,l,"ERR Field `%s' already exists!");
};

class MXPut {
 public:

    MXPut(const char *vn=0, const char *ws=0)
     : S(0), vname(0),F(0), L(0) {
       strncpy(wsp, ws && ws[0]? ws:"caller",wlen-1); wsp[wlen-1]=0;
       init(0,0,vn);
    };

    MXPut(const char *file, int line, const char *vn="ans", const char *ws=0)
     : S(0), vname(0),F(0), L(0) {
       strncpy(wsp, ws && ws[0]? ws:"caller",wlen-1); wsp[wlen-1]=0;
       init(file,line,vn);
    };

   ~MXPut() { put(); }

    void init() {
       if (F) { WB_DELETE(F); }; L=0;
       if (vname) { WB_DELETE(vname); }
       if (S) {
          if (mxGetNumberOfFields(S)) { const char *v="I_check_this";
             wblog(F_L,"WRN %s() clearing non-empty S (see %s)",FCT,v);
             mxPutAndDestroy(F,L,S,v,"caller");
          }
          else { mxDestroyArray(S); }
          S=NULL;
       }
    };

    void init(const char *file, int line, const char *vn=NULL) {
       if (S || F || vname) init();
       if (file && file[0]) { L=line;
          WB_NEW(F,strlen(file)+1); strcpy(F,file);
       }

       if (vn && vn[0])
            { WB_NEW(vname,strlen(vn)+1); strcpy(vname,vn); }
       else { WB_NEW(vname,4); strcpy(vname,"ans"); }

      #pragma omp critical (using_MEX_API)
       S=mxCreateStructMatrix(1,1,0,NULL);

       if (L) addFL();
    };

    void put(const char *ws=0, const char *vn=0) {
       if (S) {
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name",FCT,vname);
          const char s0[]="ans", *s=((vn && vn[0]) ? vn : vname);
          if (!s || !s[0]) s=s0;

          if (!ws || !ws[0]) {
             if (!wsp[0]) wblog(FL,"ERR %s() workspace not set",FCT);
             ws=wsp;
          }

          mxPutAndDestroy(F,L,S,s,ws);
          S=NULL;
       }
       init();
    };

    void save2(mxArray* &S0) { S0=S; S=NULL; init(); };

    mxArray* toMx() { mxArray* S0=S; S=NULL; init(); return S0; };

    void add2Struct(mxArray *S0) {
       if (S) {
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name !?",FCT,vname);
          mxAddField2Scalar(F_L,S0,vname,S); S=NULL;
       }
       init();
    };

    MXPut& addFL() {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       else {
          if (!vname) wblog(F_L,
             "ERR %s() got empty variable name !?",FCT,vname);
          if (mxGetFieldNumber(S,"FL")<0)
          mxAddField2Scalar(F_L,S,"FL",wbstring(shortFL(F_L)).toMx());
          else wblog(F_L,"WRN %s() field FL already set",FCT);
       }
       return *this;
    };

    template<class T>
    MXPut& add(const T&x, const char *fn) {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       mxAddField2Scalar(F_L,S,fn,x.toMx());
       return *this;
    };

    template<class T>
    MXPut& addP(
       const T* x     QS_UNUSED_VAR,
       const char *fn QS_UNUSED_VAR
     ) {
       if (!S) wblog(F_L,"ERR MXPut::%s() S not initialized yet",FCT);
       else wblog(F_L,
          "ERR %s() not defined yet for pointer to %s",FCT,TSTR(T));
       return *this;
    };

    mxArray *S;
    char wsp[16]; 
    char *vname, *F;
    unsigned L;

    static const unsigned wlen;

 protected:
 private:
};

   const unsigned MXPut::wlen=16;

template<>
MXPut& MXPut::addP(const mxArray* a, const char *fn) {
   if (!S    ) wblog(F_L,"ERR MXPut::%s() S not set yet",FCT);
   if (!vname) wblog(F_L,"ERR %s() got empty variable name",FCT,vname);
   mxAddField2Scalar(F_L,S, fn, (mxArray*)a);
   return *this;  
};

#ifndef NOMEX
template<>
MXPut& MXPut::addP(const char* s, const char *fn) {
   return addP(mxCreateString(s),fn);
};
#endif

template<>
MXPut& MXPut::add(const double &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const float &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const long double &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const wbcomplex &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const char &x, const char *fn) {
   return addP(numtoMx(x),fn);
};

template<>
MXPut& MXPut::add(char* const &s, const char *fn) {
   return add(wbstring(s),fn);
};
template<>
MXPut& MXPut::add(const char* const &s, const char *fn) {
   return add(wbstring(s),fn);
};

template<>
MXPut& MXPut::add(const unsigned &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const int &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const bool &x, const char *fn) { 
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const long &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
template<>
MXPut& MXPut::add(const unsigned long &x, const char *fn) {
   return addP(numtoMx(x),fn);
};

#ifdef __WB_MPFR_HH__
template<>
MXPut& MXPut::add(const Wb::quad &x, const char *fn) {
   return addP(numtoMx(x),fn);
};
#endif

class mxArray_buf { 

 public:
    mxArray_buf(const char *n=NULL) {
       BUF.max_load_factor(0.2);
       if (n && n[0])
            { name = new char[strlen(n)+1]; strcpy(name,n); }
       else { name = new char[4]; strcpy(name,"MXB"); }
    };

   ~mxArray_buf() {
       size_t l=BUF.size();
       if (l) wblog(FL,"WRN %s() %s got %d entries still",PROG,name,l);
       BUF.clear();
       if (name) { delete [] name; name=NULL; }
    };

    void add(void *p, mxArray *a) {
       if (p) {
          #pragma omp critical (using_mxArray_P2X)
          { auto im=BUF.insert({p,a});
            if (!im.second) wblog(FL,
               "ERR %s() %s already got entry for %p",PROG,name,p);
          }
       }
       else wblog(FL,"WRN %s() got null pointer",FCT,p);
    };

    mxArray* find(const char *F, int L, void *p) {
       mxArray *a=NULL;

       #pragma omp critical (using_mxArray_P2X)
       { auto im=BUF.find(p);
         if (im!=BUF.end()) { a=im->second; }
       }
       if (!a && F) wblog(F,L,
          "WRN %s() %s missing entry for %p",PROG,name,p);
       return a;
    };

    mxArray* Return(const char *F, int L, void *p) {
       mxArray *a=NULL;

       #pragma omp critical (using_mxArray_P2X)
       { auto im=BUF.find(p);
         if (im!=BUF.end()) {
            a=im->second; BUF.erase(im->first); 
         }
       }
       if (!a && F) wblog(F,L,
          "WRN %s() %s missing entry for %p",PROG,name,p);
       return a;
    };

    int erase(void *p) {
       int e=0;
       #pragma omp critical (using_mxArray_P2X)
       { e=BUF.erase(p); } 
       return e;
    };

    unordered_map<void*, mxArray*> BUF;

    char* name;

 protected:
 private:

};

   mxArray_buf P2X("P2X");

#endif

