/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wblog (logging routines)
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

#ifndef __WB_LOG_H__
#define __WB_LOG_H__

// e.g. set in matlab via
// >> setenv('WB_VERBOSE',num2str(1<<10 + 1<<21));
// >> setenvb WB_VERBOSE  10 21
// read in wblib.h -> envVRB=gathered.cc -> get_WB_VERBOSE()
// default 0xF1; // 0xF1 = 241: bit 0 and 4..7

#define WBL_MSEC__   (1U<<10) 
#define WBL_TCAST__  (1U<<12) 

#define WBL_CLCK__   (1U<<20) 
#define WBL_CTR__    (1U<<21) 
#define WBL_MMEX__   (1U<<22) 

#define WBL_IO__     (1U<<23) 
#define WBL_RESIZE__ (1U<<24) 
#define WBL_RLARGE__ (1U<<25) 

#define WBLOG_MSEC   (Wb::envVRB & WBL_MSEC__  )
#define WBLOG_TCAST  (Wb::envVRB & WBL_TCAST__ )
#define WBLOG_CLCK   (Wb::envVRB & WBL_CLCK__  ) 
#define WBLOG_CTR    (Wb::envVRB & WBL_CTR__   )
#define WBLOG_MMEX   (Wb::envVRB & WBL_MMEX__  )
#define WBLOG_IO     (Wb::envVRB & WBL_IO__    )
#define WBLOG_RESIZE (Wb::envVRB & WBL_RESIZE__)
#define WBLOG_RLARGE (Wb::envVRB & WBL_RLARGE__)

#define WBL_HLEN 20
#define WBL_GOT_SHORTFL -99

const char* shortFL(
   const char *F, int L=-1, unsigned n=WBL_HLEN,
   const char *P=NULL, char sep=':');

const char* tmpmatFL(
   const char *F, int L, const char *vn, char *s, unsigned n);

   enum WBLOG_TYPE { OK, WRN, ERR, NUM_TYPE }; 

   const char* WBLOG_TYPE_STR[NUM_TYPE] = { "OK","WRN","ERR" };

namespace Wb {

   struct LogException : public exception { 
     public:
       LogException (WBLOG_TYPE l=OK, const char *s=NULL);

       LogException (const LogException &e)
        : type(e.type), nrefs(e.nrefs), count(e.count), ith(e.ith), nth(e.nth) {
          if (nrefs) { ++(*nrefs); ++count;
          }
          istr[63]=0;
          strncpy(istr,e.istr,64);
          if (istr[63]) { istr[61]=istr[62]='.'; istr[63]=0; }
       };

      ~LogException() noexcept {
          if (nrefs) {
          if (--(*nrefs)<1) {
             delete nrefs; nrefs=NULL;
          }
          else {
          }
       }};

       LogException& operator=(const LogException &e) {
          type=e.type; count=e.count; ith=e.ith; nth=e.nth;
          nrefs=e.nrefs; if (nrefs) { ++(*nrefs);
          }
          istr[63]=0;
          strncpy(istr,e.istr,64);
          if (istr[63]) { istr[61]=istr[62]='.'; istr[63]=0; }
          return *this;
       };

       LogException& operator++() { 
          ++count;

          int n=omp_get_num_threads();
          if (nth<n) { nth=n; ith=omp_get_thread_num(); }

          return *this;
       };

       LogException& operator+=(const LogException &e) { 

          if (!type) 
               { (*this)=e; count=1; }
          else { ++count; }

          int n=omp_get_num_threads();
          if (nth<n) { nth=n; ith=omp_get_thread_num(); }

          return *this;
       };

       bool operator!() const { return (!type && !count); } 
       explicit operator bool() const { return (type || count); } 

       LogException& init() {
          if (nrefs) { --(*nrefs);
             if (*nrefs<1) { delete nrefs; nrefs=NULL; }
          }
          type=OK; count=ith=nth=0; istr[0]=0;
          return *this;
       };

       const char* what() const throw () { return istr; };

       wbstring toStr() const;
       void report(const char *F=0, int L=0, const char *func=0);

       WBLOG_TYPE type;
       int *nrefs;
       char istr[64];
       unsigned count; 
       int ith, nth;   

     private:
   };

   char* surdStrf(wbstring &s, const char *fmt, ...);
   char* surdStrf(char *s, unsigned n, const char *fmt, va_list args);

enum class TCOLS {
   BLACK, RED,  GREEN, YLW,   BLUE, PINK, CYAN, LGRAY, DGRAY,
   LRED, LGREEN, LYLW, LBLUE, LPNK, LCYAN, WHITE
};

class termcolor { 
  public:

   termcolor(const char *tag)
    : em(NULL), e1(NULL) { int k=ID(tag); init(k); };

   termcolor(Wb::TCOLS k)
    : em(NULL), e1(NULL) { init(unsigned(k)); };

   termcolor(unsigned k)
    : em(NULL), e1(NULL) { init(k); };

  ~termcolor() {
     if (em) { delete [] em; }
     em=e1=NULL;
   }

   termcolor& init(unsigned k);

   static int ID(const char *tag);

   char *em, *e1;

  private:
};

}; 

namespace wblog { 

class SBUF {      
  public:
    SBUF(unsigned l0=wblog::SLEN) : l(0), slen(0), sbuf(NULL) {
       init(l0);
    };

   ~SBUF() { if (sbuf) {
      #pragma omp critical (__ensure_sequential_wblog__)
       { flush(stderr,1); }
       delete [] sbuf;
    }};

    SBUF& init(unsigned l0=0) { l=0;
       if (sbuf) { delete [] sbuf; sbuf=NULL; }
       if ((slen=l0)) {
          sbuf = new char[slen];
          if (!sbuf) { sprintf(str,"ERR %s() "
             "failed to allocate sbuf (l=%d)",__FUNCTION__,slen);
             ExitMsg(str);
          }
          memset(sbuf,0,slen);
       }
       return *this;
    };

    int cat(const char *s);
    int catf(const char *fmt, ...);

    void error_bounds(const char *F, int L, const char *fct, const char *fmt);
    void increase_size(const char *F, int L, unsigned n);

    void flush(FILE *fid=stdout, char fflag=0);
    void print(const char *F=0, int L=0, const char *istr=0);

    int skipEscCols(); 

    unsigned l,slen;
    char *sbuf; 

  private:
};

class stdio_buf {  
 public:
    stdio_buf(SBUF &S, FILE *f_, int t_) : s(NULL), fid(0), tid(0), ref(0) {
       if (f_!=stdin && f_!=stdout) {
          fprintf(stderr,"\nERR %s() invalid fid=%p\n'%s'\n\n",
          FCT,f_,S.sbuf?S.sbuf:"");
       }
       if (S.sbuf) {
          s=S.sbuf; fid=f_; tid=t_; ref=1; 
          S.sbuf=0; S.slen=S.l=0;
       }
    };

    stdio_buf(const stdio_buf &b) : s(NULL), fid(0), ref(0) {
        if (!b.ref) fprintf(stderr,"\nERR %s() got b.ref=%d\n\n",FCT,b.ref);
        s=b.s; fid=b.fid; tid=b.tid;
    };

    void print_stdout(FILE *f=NULL) const {
       if (Wb::envVRB &256) {
          unsigned i=0, m=0; for (; s[i]; ++i) {
             if (s[i]!='\n') { if (++m>2) break; } else m=0;
          }
          if (m>2) { i-=2; char *x = (char*)s; 
             if (tid<10) x[i]=tid   +'0'; else
             if (tid<36) x[i]=tid-10+'a'; else
             if (tid<62) x[i]=tid-36+'A'; else { x[i]='*'; }; x[i+1]=' ';
          }
       }
       if (f) { fprintf(f,"%s",s); }
       else {
          if (Wb::envDKT) 
               { mexPrintf("%s",s); }
          else { PRINTF("%s",s); }
       }
    };

   ~stdio_buf() { if (s && !ref) {
       if (Wb::my_caller_tid!=omp_get_thread_num()) {
          fprintf(stderr,"WRN %s() got tid=%d/%d !?\nWRN '%s'\n",
          FCT,omp_get_thread_num(),Wb::my_caller_tid,s);
       }
       else {
          print_stdout();
          delete [] s; 
          fflush(fid);
       }
    }};

  private:
    const char *s;
    FILE *fid;
    int tid;
    char ref;
};

   std::deque<stdio_buf> myIO;

}; 

   void wbSetLogLevel(unsigned l);

   int wblogs(wblog::SBUF &S,
      WBL_COLOR_SCHEME xcol, 
      const char* file, int line,
      const char *fmt, va_list args,
      char Hflag=0 
   );

   int wblogs(
      wblog::SBUF &S, WBL_COLOR_SCHEME xcol,
      const char* file, int line, const char *fmt, ...
   ){
      va_list args; Wb::ARGV wd(&args); 
      va_start(args,fmt);
      return wblogs(S,xcol,file,line,fmt,args);
   };

   int wblogf(FILE *fid, 
      const char* file, int line, const char *fmt, va_list args);

   int wblogf(FILE *fid,
      const char* file, int line, const char *fmt, ...) {

      va_list args; Wb::ARGV wd(&args); 
      va_start(args,fmt);

      int l=wblogf(fid,file,line,fmt,args); 

      return l;
   };

   int wblog1(const char* file, int line, const char *fmt, ...);

#define wblog(...) wblogf(stdout, __VA_ARGS__)

   int wb_printf(const char *fmt, ...) {   

      unsigned l=0, n=256;
      wblog::SBUF S(n);  {
         va_list args; Wb::ARGV wd(&args); 
         va_start(args,fmt);
         l=vsnprintf(S.sbuf,n,fmt,args);
      }
      if (l>=n) fprintf(stderr,
         "\nERR string out of bounds (%d/%d)\nERR %s\n\n",l,n,S.sbuf);

      #pragma omp critical (__ensure_sequential_wblog__)
       { S.flush(); }

      return l;
   };

   unsigned wblog_checktag(const char *fmt, const char *t, char *tag);

   char wblog_findtoken(const char *istr, const char *tok, int maxoffset=-1);

   void banner(
      unsigned n, const char *s,
      const char *istr=NULL, const char *fstr=NULL
   );

#endif

