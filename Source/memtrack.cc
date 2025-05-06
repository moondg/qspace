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

#ifndef __WB_MEM_TRACK_CC__
#define __WB_MEM_TRACK_CC__

// NB! included in wblib.h only if __WBDEBUG__ is set
// see also memlib.hh (always included)

#ifdef __WBDEBUG__
   #pragma message "NB! got WBDEBUG mode"
#endif
#ifdef __WB_MEM_CHECK__
   #pragma message "NB! got WB_MEM_CHECK mode"
#endif
#ifdef EXTENDED_MEM_CHECK
   #pragma message "NB! got EXTENDED_MEM_CHECK mode"
#endif

namespace Wb {

class MREC { 
   public:

#ifdef EXTENDED_MEM_CHECK
      MREC() : id(0), len(0), unit(0), fline(0) {};
     ~MREC() { if (fline) { delete [] fline; fline=0; }};
#else
      MREC() : id(0), len(0), unit(0) {};
     ~MREC() {};
#endif

#ifdef EXTENDED_MEM_CHECK
      template <class T>
      const char* init(
         const char *F, int L,
         size_t i, 
         const T* p QS_UNUSED_VAR 
         size_t l  
      ){
         id=i; unit=sizeof(T); len=l;
         if (F) {
            const std::type_info &tid = typeid(T);
            const char *s=shortFL(F,L), *t=tid.name();
            unsigned n=strlen(s) + strlen(t) + 5;
            if (fline) { delete [] fline; }
            fline = new char[n];
            snprintf(fline,n,"%s (%s)",s,t);
         }
         return fline;
      };
#else
      template <class T>
      const char* init(
         const char *F QS_UNUSED_VAR,
         int L         QS_UNUSED_VAR,
         size_t i, 
         const T* p    QS_UNUSED_VAR, 
         size_t l 
      ){
         id=i; unit=sizeof(T); len=l;
         return 0;
      };
#endif

      void println(size_t i, void *p) const;

      size_t id;
      size_t len;
      unsigned unit; 

#ifdef EXTENDED_MEM_CHECK
      char *fline;
#endif

   protected:
   private:
};

class WbListMTRACK { 
   public:
     WbListMTRACK() : idx(0), totsize(0), maxsize(0) {
#ifndef __WB_MEM_QUIET__
        wblog(FL, "MTR ************ MEM TRACK : ON ********************");
#endif
     };

     template<class T>
     void rm_ptr(const char* F, const int &L, T* &p) { if (p) {
        #pragma omp critical (manage_MTRACK) 
        { map<void*,MREC>::iterator im = M.find((void*)p);
          if (im == M.end()) wblog(F_L,"WRN memory to be freed "
             "not in list (0x%lx; #%d) !!?",p,omp_get_thread_num());
          else {
             totsize -= (im->second.len) * (im->second.unit);
             M.erase(im);
          }
        }
     }};

     template<class T>
     void add_ptr(const char* F, int L, T* &p, size_t l) {

        if (p==NULL) wblog(F_L,"ERR got NULL pointer (out of memory !!?)");

        #pragma omp critical (manage_MTRACK) 
        {  MREC &r=M[(void*)p];

           if (r.id!=0) wblog(F_L,"ERR entry with pointer "
              "%lx (%ld,%d,%d) already exists!!",p,r.id,r.len,r.unit);

           r.init(F_L,++idx,p,l);

           totsize += r.unit*r.len; 
           checkTotSize();
        }
     };

     template<class T>
     bool exists(T* &p) {
        map<void*,MREC>::iterator im;
        #pragma omp critical (manage_MTRACK) 
        { im = M.find(p); }
        return (im!=M.end() ? 1 : 0);
     };

     template<class T>
     void tst_exists(const char* F, int L, T* &p) {
        map<void*,MREC>::iterator im;

        #pragma omp critical (manage_MTRACK) 
        { im = M.find(p); }

        if (im!=M.end()) {
           const MREC &r=im->second;
           wblog(FL,"TST 0x%12lx exists in list (%d @ %d = %d)",
              p, r.len, r.unit, r.len*r.unit);
        }
        else { wblog(FL,"TST 0x%12lx does not exist in list",p); }
     };

     void checkTotSize() { if (totsize>maxsize) {
        size_t stot=(totsize>>30);
        if (stot) {
           if (!(stot>>2)) {
              size_t smax=(maxsize>>30);
              while (smax) { smax>>=1; stot>>=1; }
              if (stot) wblog(FL,
              " *  %s[] reached %s  ",FCT,totStr());
           }
           else {
              size_t smax=(maxsize>>32); stot>>=2;
              if (stot>smax) wblog(FL,
              " *  %s[] reached %s  ",FCT,totStr());
           }
        }
        if (totsize>maxsize) 
        maxsize = totsize;
     }};

     char* totStr(char mflag=0) { size_t n=16, l=0;

        if (mflag) { char s1[n], s2[n];
           Wb::memsize2Str(totsize,s1,n);
           Wb::memsize2Str(maxsize,s2,n); l=snprintf(sbuf,32,"%s (%s)",s1,s2);
        }
        else {
           Wb::memsize2Str(totsize,sbuf,32); l=strlen(sbuf);
        }

        if (l>=32) wblog(FL,
           "ERR %s() string out of bounds (%d/32)",FCT,l);
        return sbuf;
     };

     void printSize(const char *F, int L, char mflag=0) {
        if (mflag) { int n=16; char s1[n], s2[n];
           wblog(F_L,"MEM tracking %d entries @ %s (max. %s) ", M.size(),
           Wb::memsize2Str(totsize,s1,n),
           Wb::memsize2Str(maxsize,s2,n));
        }
        else { int n=16; char s1[n];
           wblog(F_L,"MEM tracking %d entries %s ",
           M.size(), Wb::memsize2Str(totsize,s1,n));
        }
     };

     map<void*, MREC> M;

     unsigned long long idx;
     unsigned long long totsize;
     unsigned long long maxsize;

   protected:
   private:

     char sbuf[32];
};

WbListMTRACK gML;

template<class T> inline
void NEW(const char* F, int L, T* &p, size_t n) {

   if (n) {
      try { p = new T[n]; } catch (...) { p=NULL; }
      if (!p) {
         MemStat(FL); 
         wblog(FL,"ERR out of memory (%ld * %d = %.1fG) !?",
            n, sizeof(T), n*sizeof(T)/double(1<<30)
         );
      }

      gML.add_ptr(F_L,p,n);
   }
   else {
      wblog(FL,"WRN allocate empty space (%s, %d)",shortFL(F,L),n);
      p=NULL;
   }
};

template<class T> inline
void NEW_1(const char* F, int L, T* &p) { 

   try { p = new T; } catch (...) { p=NULL; }
   if (!p) {
      MemStat(FL); 
      wblog(F_L,"ERR failed to allocate instance");
   }
   gML.add_ptr(F_L,p,1);
};

int MemCheck(const char *F="", int L=0,
    const char *task="",
    long dl=0, long db=0
){
    static unsigned long long nbytes=0, iref=0;
    static size_t len=0;

    if (task!=0 ? task[0]==0 : task==NULL) {
       wblog(F,L, "MTR info: %d entries (%s)",
       gML.M.size(), gML.totStr());
       return 0;
    }

    if (!strcmp(task,"start")) {
       len    = gML.M.size();
       nbytes = gML.totsize;
       iref   = gML.idx;
    }
    else if (!strcasecmp(task,"stop")) {
       size_t l=gML.M.size();
       unsigned long long b=gML.totsize;

       int m=0, e = long(l-len) != dl || long(b-nbytes) != db;

       if (e) {
          map<void*,MREC>::iterator im;
          for (im=gML.M.begin(); im!=gML.M.end(); ++im) {
              if (im->second.id>iref) { if (!(m++)) printf("\n"); if (m<20)
                 printf("%6d: %9ld %10lx %8d /%d\n",m,
                 (unsigned long)im->second.id,
                 (unsigned long)im->first, im->second.len, im->second.unit);
              }
          }
          if (m>20) wblog(FL,
             "%NTST %s() %d/%d entries skipped%N",FCT,m-20,l);
          else if (m) printf("\n");
       }

       wblog(F,L,
         "%s  %d entries (%+d) using %s (%+ld)",
          !strcmp(task,"STOP") ? (e ? "ERR" : "SUC") : "MTR",
          l, l-len, gML.totStr(), b-nbytes
       );
    }
    else if (!strcasecmp(task,"info")) {
       size_t l=gML.M.size(); long s=gML.totsize-nbytes;
       int n=16; char ss[n];

       wblog(F,L,
         "MTR %d entries (%+d) using %s (%s%s; %lld)",
          l, l-len, gML.totStr(), s>=0 ? "+":"",
          Wb::memsize2Str(s,ss,n), gML.idx
       );
    }
    else if (!strcasecmp(task,"LIST")) {
       size_t i; map<void*, MREC>::iterator I;
       MemCheck(F,L,"info",dl,db);
          if (!gML.M.size()) { return 0; }

       printf(" entry data_location internal_id  memory_size\n");
       for (i=0, I=gML.M.begin(); I!=gML.M.end() && i<201; ++I) {
           I->second.println(++i,I->first);
       }
       if (I!=gML.M.end())
       printf("\n  ... and %ld others\n\n",long(gML.M.size()-i));
    }
    else wblog(F,L, "ERR MTRACK invalid flag >%s<", task);

    return 0;
};

void MREC::println(size_t i, void *p) const {

    int n=32; char s[n];
    Wb::memsize2Str(len*unit,s,n);

    printf("%6d %12lx %12ld %8d @%3d = %12s %s\n",
    i, (unsigned long)p, id, len, unit, s,
    #ifdef EXTENDED_MEM_CHECK
        r.fline ? r.fline : ""
    #else
        ""
    #endif
    );
};

}; 

#ifndef __WB_MEM_QUIET__

class wbdebug_dummy {
  public:

    wbdebug_dummy() {
       wblog(FL,"MTR ************ DEBUG MODE: ON ********************");
       Wb::MemCheck(FL,"info"); printf("\n");
    };

   ~wbdebug_dummy() {
       if (Wb::gML.M.size()) wblog(FL,
          "WRN %s() got memory leak @ %d entries !?",FCT,Wb::gML.M.size());
       Wb::MemCheck(FL,"LIST");
    };

  protected:
  private:
};

   wbdebug_dummy wbdebug_dummy_var;

#endif

#endif

