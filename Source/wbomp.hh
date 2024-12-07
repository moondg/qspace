/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace OMP routines
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

#ifndef __WB_OMP_HH__
#define __WB_OMP_HH__

#ifdef QS_USING_OMP

   int OMP_NUM_THREADS=0; 

   int QSP_NUM_THREADS=0; 

namespace Wb {

wbstring ompID2Str(char vflag=0);

int ompStatus(const char *F, int L); 

class ompNLock { 

 public:

    ompNLock(const char *s=NULL)
     : level(0), owner(-1), istr(0), nactive(0), ntot(0) {
       omp_init_nest_lock(&_nlk);
       if (s && s[0]) { set_istr(s); }
    };

   ~ompNLock() {
       if (level || nactive) {
          wblog(FL,"WRN %s() got lock @ l=%d, n=%d (id=%d/%d) !?",
             FCT,level,nactive,omp_get_thread_num(),omp_get_num_threads());
          if (level==1 && !nactive) { release(); } 
       }
       if (istr) { delete [] istr; istr=0; }
       owner=-1; 
       omp_destroy_nest_lock(&_nlk);
    };

    void set_istr(const char *s);

    wbstring toStr(const char *istr_=NULL) const; 

    wbstring get_istr() const;

    int acquire();

    int test_acquire();

    int release();

    int level;   
    int owner;   

    char *istr;  

    int nactive; 

    long ntot; 

 private:

    omp_nest_lock_t _nlk;   

    ompNLock(const ompNLock &);
    void operator=(const ompNLock &);

#ifdef LD_CLEBSCH_QS
    wbstring check_BUF_size() const;
#endif
};

class ompGuard { 

 public:

    ompGuard() : lp(0) { };

    ompGuard(const char *F, int L, ompNLock &lk, int lmax=99)
     : lp(0) {
       if (omp_get_num_threads()>1) { 
          acquire(F,L,lk,lmax);
       }
    };

    ompGuard(ompNLock &lk, int lmax=99) : lp(0) {
       if (omp_get_num_threads()>1) { 
          acquire(0,0,lk,lmax);
       }
    };

    void acquire(ompNLock &lk, int lmax=99) { acquire(0,0,lk,lmax); };
    void acquire(const char *F, int L, ompNLock &lk, int lmax=99);

   ~ompGuard() { if (lp) { lp->release(); }};

    int level() { return (lp ? lp->level : -1); };

 private:

    ompNLock *lp;  

    ompGuard(const ompGuard &) {
       wblog(FL,"ERR %s() forbidden copy",FCT); };
    void operator=(const ompGuard &) {
       wblog(FL,"ERR %s() forbidden copy",FCT); };
};

class nthreads_aux { 
 public:
    nthreads_aux() : np_curr(0), np_max(0) { };

    nthreads_aux& increment(const char *fl) {
       int q=omp_get_num_threads();
       if (np_max>q) {
          Wb::termcolor Tc(1); 
          fprintf(stdout,"TST %s%s() having nq=%d->%d%s\n",
             Tc.e1,fl? fl:"(null)",np_max,q,Tc.em);
          np_max=q;
       }
       if (np_max<(++np_curr)) {
          Wb::termcolor Tc(5); 
          fprintf(stdout,"TST %s%s() having nq=%d->%d%s\n",
             Tc.e1,fl? fl:"(null)",np_max,np_curr,Tc.em);
          np_max=np_curr;
       }
       return *this;
    };

    nthreads_aux& decrement(const char *fl) {
       if (--np_curr<0) { Wb::termcolor Tc(5); 
          fprintf(stdout,"%sWRN %s() having nq=%d/%d%s\n",
          Tc.e1,fl? fl:"(null)",np_curr,np_max,Tc.em);
       }
       return *this;
    };

 private:
    int np_curr, np_max;
};

   map< std::string, nthreads_aux >  np_hist;

class checkMaxThreads { 

 public:

    checkMaxThreads(const char *F, int L) { 
       std::string fl=shortFL(F_L);
      #pragma omp critical (checkMaxThreads__)
       { np_hist[fl].increment(fl.c_str()); im=np_hist.find(fl); }
    };

   ~checkMaxThreads() {
      #pragma omp critical (checkMaxThreads__)
       { im->second.decrement(im->first.c_str()); }
    };

 private:

   map< std::string, nthreads_aux >::iterator im;
};

   class thread_xlink;
   map< int, thread_xlink* >  txmap;

class thread_xlink { 
 public:

    thread_xlink(const char *F=NULL, int L=0) : tx(0) {
       t0=omp_get_thread_num();
       nt=omp_get_num_threads();

       if (F) { fl=shortFL(F,L); }
    };

   ~thread_xlink() {
       auto it = txmap.find(t0);
       if (it!=txmap.end()) {
          #pragma omp critical (thread_xlink__)
          txmap.erase(t0);
       }
    };

    void init(const char *F=NULL, int L=0) {
       if (nt<=1) {
          t0=omp_get_thread_num();
          nt=omp_get_num_threads();
       }
       else { 
          unsigned t2=omp_get_thread_num(); if (t2!=t0) {
          wblog(F_L,"WRN %s() thread-id changed %d->%d / %d->%d",
          FCT, t0,t2,nt,omp_get_num_threads());
       }}

       if (!tx) { if (F || fl.empty()) fl=shortFL(F_L);
         #pragma omp critical (thread_xlink__)
          txmap[t0]=this;
       }
       tx=(1<<t0);
    };

    int add(const char *F, int L, size_t h);
    int add(const char *F, int L, const ompNLock *lk);

    static int check_deadlocks();

    unsigned t0, nt; 

    size_t tx;
    std::string fl; 

 private:
};

}; 

   Wb::ompNLock wblog_lock("wblog");     
   Wb::ompNLock rclog_lock("rclog");     
   Wb::ompNLock mxapi_lock("mxapi");     

   Wb::ompNLock XS_buf("XS_buf");

   Wb::ompNLock CS_buf("CS_buf");

   Wb::ompNLock CG_locks("CG_locks");

#else
#endif 

#endif 

