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

#ifndef __WB_OMP_CC__
#define __WB_OMP_CC__

#ifdef QS_USING_OMP

wbstring Wb::ompID2Str(char vflag) {

   unsigned i=omp_get_thread_num(), n=omp_get_num_threads();
   unsigned k=(vflag? 16:0), l=omp_get_level();
   wbvec<char> sx;

   if (!i && n==1 && !l) { sx.init(8+k); sx.cat(FL,"serial"); l=6; }
   else if (l==1)
        { sx.init(12+k); sx.catf(FL,"OMP-%d/%d",   i,n  ); } 
   else { sx.init(16+k); sx.catf(FL,"OMP-%d/%d@%d",i,n,l); } 

   if (k) {
      sx.catf(0,0,", %s",Wb::hostid(vflag).data);
   }

   return sx.data;
};

inline int Wb::omp_parallel() {
   return omp_get_level();
};

unsigned Wb::get_omp_tid_nn(
   unsigned *l_,
   char check 
){
   unsigned tid=0; 
   int q=0, l=omp_get_level(); if (l_) { *l_=l; }
   for (; l>=0; --l) {
      if ((q=omp_get_ancestor_thread_num(l))>0) {
         if (tid) wblog(FL,
            "ERR %s() got nested parallelization: tid=%d/%d (l=%d/%d)",
            FCT,q,tid,l, omp_get_level());
         tid=q; if (!check) { break; }
      }
   }
   return tid;
};

int Wb::ompStatus(const char *F, int L) {
    char fmt[]="  %-30s %2d\n";
    wbvec<char> s(1024);

    s.catf(0,0,fmt,"omp_get_active_level", omp_get_active_level());
    s.catf(0,0,fmt,"omp_get_cancellation", omp_get_cancellation());  
    s.catf(0,0,fmt,"omp_get_default_device", omp_get_default_device()); 
    s.catf(0,0,fmt,"omp_get_dynamic", omp_get_dynamic()); 
    s.catf(0,0,fmt,"omp_get_level",   omp_get_level());   
    s.catf(0,0,fmt,"omp_get_max_task_priority", omp_get_max_task_priority()); 
    s.catf(0,0,fmt,"omp_get_max_threads", omp_get_max_threads());  
    s.catf(0,0,fmt,"omp_get_max_active_levels", omp_get_max_active_levels());
    s.catf(0,0,fmt,"omp_get_num_procs", omp_get_num_procs()); 
    s.catf(0,0,fmt,"omp_get_num_teams", omp_get_num_teams()); 
    s.catf(0,0,fmt,"omp_get_num_threads", omp_get_num_threads());

    s.catf(0,0,fmt,"omp_get_proc_bind", omp_get_proc_bind());
    s.catf(0,0,fmt,"omp_get_team_num", omp_get_team_num());
    s.catf(0,0,fmt,"\n"); 
    s.catf(0,0,fmt,"omp_get_thread_num", omp_get_thread_num()); 
    s.catf(0,0,fmt,"omp_in_parallel", omp_in_parallel());
    s.catf(0,0,fmt,"omp_in_final", omp_in_final());
    s.catf(0,0,fmt,"omp_is_initial_device", omp_is_initial_device());

    fprintf(stdout,"\n%-18s %s -> %s()\n\n",shortFL(F_L),myname,FCT);
    if (s.check_bounds()>=0)
         { fprintf(stdout,"%s\n",s.data); }
    else { fprintf(stdout,"%s  ... (%d/%d)\n",s.data,s.l,s.len); }

    return 0;
};

int Wb::ompNLock::acquire() {

   omp_set_nest_lock(&lk_n);
   owner=omp_get_thread_num(); ++level; ++ntot;

 # if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 12 ) 
   LKF.blogf(FL," +  %s%s",STR(*this), check_BUF_size().data);
 # endif

   return level;
};

int Wb::ompNLock::test_acquire() {
   int l=omp_test_nest_lock(&lk_n);
   if (l) {
      owner=omp_get_thread_num();
      l=++level; ++ntot;

   }
   return l;
};

int Wb::ompNLock::release() {

   if (owner!=omp_get_thread_num()) { 
      wbstring sx=get_istr(); bool q=sx.strlen_()<10;
      wblog(PFL,"WRN %s() "
        "lock%s%s @ l=%d/%d/%d by non-owner %d/%d/%d !?%s%s",
         FCT, q? " ":"", q? sx.data:"", level, nactive, ntot,
         owner, omp_get_thread_num(), omp_get_num_threads(),
         q? "":"\n", q? "":sx.data
      );
   }

   if (level<=0) wblog(FL,"WRN %s() %s",FCT,STR(*this));

   int l=(--level); 
   if (l<=0) owner=-1; 

   #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 12 ) 
    LKF.blogf(FL," -  %s%s",STR(*this),check_BUF_size().data);
   #endif

   omp_unset_nest_lock(&lk_n);

   return (l>0? l : 0);
};

#ifdef LD_CLEBSCH_QS

wbstring Wb::ompNLock::check_BUF_size() const {

   wbstring s; 

   if (istr) { s.init(16);

      if (!strcmp(istr,"CG_locks"))
         snprintf(s.data,s.len," @ %ld", gCS.BUF.size()); else
      if (!strcmp(istr,"XS_buf"))
         snprintf(s.data,s.len," @ %ld",gXS.XBUF.size());
   }
   else { s.init(1); s.data[0]=0; }

   return s;
};

#endif

wbstring Wb::ompNLock::get_istr() const {
   wbstring s(128); 
   if (istr) {
      Wb::termcolor Tc("WRN");
      snprintf(s.data,s.len,"%s%s%s",Tc.e1,istr,Tc.em);
   }
   else {
       snprintf(s.data,s.len,"#%lx",(unsigned long)this);
   }
   return s;
};

wbstring Wb::ompNLock::toStr(const char *istr_) const { 

   wbstring sout(128); 
   unsigned l=0, n=sout.len; char *s=sout.data;

   if ( nactive && !ntot ) l+=snprintf(s,n,"%d:%d", level,nactive); else
   if (!nactive && ntot>9) l+=snprintf(s,n,"%d:%03ld", level,ntot);
   else                    l+=snprintf(s,n,"%d:%d/%ld",level,nactive,ntot);

   l+=Wb::strpad_(s+l,8-l);
   l+=snprintf(s+l,n-l," %02d",owner); 

   if (l<n) {
      if (istr && *istr) { l+=snprintf(s+l,n-l," %s", istr ); }
      else if (istr_)    { l+=snprintf(s+l,n-l,"`%s'",istr_); }
      else {
         l+=snprintf(s+l, n-l," (null)");
         wblog(FL,"WRN %s() got null istr",FCT);
      }
   }
   if (l<n && (!istr || !*istr)) {
      l+=snprintf(s+l, n-l," %p",this); 
   }
   return sout;
};

void Wb::ompNLock::set_istr(const char *s) {

   if ((!istr || !*istr) && s && *s) {
      if (istr) { delete [] istr; }
      istr = new char[strlen(s)+1];
      strcpy(istr,s);
   }
};

void Wb::ompGuard::acquire(const char *F, int L, ompNLock &lk, int lmax) {

   if (!Wb::omp_parallel()) { return; } 

   if (lp) wblog(FL,"ERR %s() already initialized with lock",FCT);

   lk.acquire();
   lp=&lk;

   if (lk.ntot<=3 && F) { Wb::termcolor Tc("WRN");
      fprintf(stdout,"\n  %sTST %s %s() "
        "using ompCriticial `%s'%s\n\n",Tc.e1,shortFL(F_L),FCT,
         lk.istr? lk.istr:"(unnamed ompNLock)",Tc.em);
      fflush(0);
   }

   if (lk.level>lmax) {
      char sx[48]; snprintf(sx,48,"got lk.level=%d/%d",lk.level,lmax);
      if (lk.level<lmax+8)
           { fprintf(stdout,"%s %s WRN %s\n",shortFLT,sx); }
      else { fprintf(stdout,"%s %s ERR %s\n",shortFLT,sx);
         throw Wb::LogException(ERR);
      }
   }
};

wbstring Wb::ompGuard::toStr(const char *istr_) const {
   return (lp ? lp->toStr(istr_) : "(null)");
};

#ifdef LD_CLEBSCH_QS

namespace CG {

   Wb::ompNLock* Guard::find(size_t h) { 
      Wb::ompGuard myLK(CG_locks);
      auto it = CG::lock_map.find(h);
      if (it==CG::lock_map.end())
           { return NULL; }
      else { return it->second; }
   };

   template<class TQ>
   wbstring Guard::Status(const QSet<TQ> &Q, const char *tag) {
      size_t h1=QHash<TQ>()(Q,tag); wbstring istr=Q.toStr();
      if (tag && tag[0]) {
         wbvec<char> sx(strlen(tag)+4);
         sx.catf(0,0," '%s'",tag); istr+=sx.data;
      }
      return Status(h1,istr.data);
   };

   template<class TM>
   wbstring Guard::Status(const cgc_contract_id<TM> &idc) {
      size_t h1=MHash<TM>()(idc); wbstring istr=idc.toStr();
      return Status(h1,istr.data);
   };

   wbstring Guard::Status(size_t h1, const char *istr) {
      wbstring sout; 
      const Wb::ompNLock *lk=Guard::find(h1);
      if (lk) { sout=lk->toStr(istr); }
      else {
         unsigned n=(istr? strlen(istr):0); n=MIN(2U,n)+16; sout.init(n);
         snprintf(sout.data,n,"no entry for %s",istr && istr[0]? istr:"`'");
      }
      return sout;
   };

   int Guard::toStr(wbstring &s, size_t h, const char *istr_) {
      const Wb::ompNLock *lk=Guard::find(h);
      if (lk) { s=lk->toStr(istr_); return 1; } 
      else if (istr_ && istr_[0]) { s.init(128);
         if (h)
              { snprintf(s.data,s.len,"(undefined) %s",istr_); }
         else { snprintf(s.data,s.len,"(undefined; h=0 !?) %s",istr_); }
      }
      else { s.init(48);
         if (h)
              { snprintf(s.data,s.len,"(undefined; h=0x%lX)",h); }
         else { snprintf(s.data,s.len,"(undefined; h=0)"); }
      }
      return 0;
   };

   wbstring Guard::toStr(size_t h, const char *istr_) {
      wbstring s; 
      Guard::toStr(s,h,istr_);
      return s;
   };

   wbstring Guard::toStr() const {
      wbstring sout; if (!hid.len) return sout; 

      unsigned i=0,j, k=0, n=128*hid.len; wbstring s1; char *s;
      sout.init(n); s=sout.data;

      for (; k<hid.len; ++k) {
         Guard::toStr(s1,hid[k]);
         for (j=0; s1.data[j]; ++j, ++i) { if (i<n) s[i]=s1.data[j]; }
         if (i<n && k+1<hid.len) { s[i++]='\n'; }
      }
      s[i<n ? i:n-1]=0;
      if (i>=n) wblog(FL,"ERR string out of bounds (%d/%d)\n%s",i,n,s);
      return sout;
   };
};

int CG::Guard::print_CG_locks(const char *F, int L) {

   Wb::ompGuard myLK(CG_locks);

   if (CG::lock_map.size()) { unsigned i=0;
      fprintf(stdout,"\n%s %s() *tid=%d/%d (systid=%d) got %ld CG::locks total\n"
        "## level/n_active/n_use  owner-thread  source::line "
        "lock-info-string\n\n", shortFL(F_L), FCT,
         omp_get_thread_num(), omp_get_num_threads(),
         Wb::system_tid(), CG::lock_map.size());

      for (auto I=CG::lock_map.begin(); I!=CG::lock_map.end(); ++I) {
         fprintf(stdout,"%2d. %s\n",++i,STR_(I->second));
      }; fprintf(stdout,"\n");
   }

   return 0; 
};

const char* CG::Guard::sprintf_locks( 
   wbstring &sout, const int* got, const char *ss[], const char *indent
){
   if (hid.len) {
      unsigned i=0, l=0, n=256*hid.len; char *s; wbstring sx;
      if (sout.len<n) { sout.init(n); } else { n=sout.len; }
      s=sout.data;

      for (; i<hid.len; ++i) {
         if (!Guard::toStr(sx,hid[i],"")) { sx.init(); }
         if (l<n)
         l+=snprintf(s+l,n-l,"%s%d. %-8s %s%s\n", indent? indent:"",
            i+1, got[i]? "ok.":"rejected", sx.len? "   ":"", ss[i]);
         if (l<n && sx.len)
         l+=snprintf(s+l,n-l,"%s   %s\n", indent? indent:"", sx.data);
      }
      if (l>=n) wblog(FL,
         "WRN %s() string out of bounds (%d/%d)%N%N`%s'%N", FCT,l,n,s);
   }
   else { sout=""; }

   return sout.data;
};

int CG::Guard::deal_with_cond(
   const char *F, int L, size_t h0, const char *cond) {

   int ntry=-1; if (!cond || !cond[0]) return ntry;
   unsigned i=0;

   if (cond[0]=='?') { 
      { Wb::ompGuard myLK(CG_locks);
        const auto lp=CG::lock_map.find(h0);
        if (lp==CG::lock_map.end() || !lp->second->level) { ntry=0; }
      }; i=1;
   }
   if (ntry && cond[i]) {
      if (isdigit(cond[i])) { sscanf(cond+i,"%d",&ntry); }
      else wblog(F_L,"ERR %s() invalid cond='%s'",FCT,cond);
   }

   return ntry;
};

int CG::Guard::acquire(
   const char *F, int L, const CG::FileLock &flk,
   const char *tag, const char *cond 
) {
   if (!Wb::omp_parallel()) { return 0; }

   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   size_t h0;
   wbvec<char> s(64);
   int ntry=-1; const char *sR; const char* ss[1]={s.data};

   flk.ensureFOpen(FL);

   sR=Wb::rstrstr(flk.fname.data,"/RCS"); 
   if (sR) { ++sR; }
   else { const char *s0=flk.fname.data;
      sR=Wb::rstrstr(s0,"Store");
      if (sR>s0) { --sR; }
      else { sR=s0; } 
   }

   h0 = Wb::string_hash(tag ? tag : mytag)
      ^ Wb::string_hash(sR);

   if (cond && cond[0]) {
      ntry=deal_with_cond(F,L,h0,cond);
      if (!ntry) { return 0; }
   }

   if (tag && tag[0]) 
        { s.catf(0,0,"%s", shortFL(F_L,-1,  tag,'*')); }
   else { s.catf(0,0,"%s", shortFL(F_L,-1,mytag,'/')); }
   s.pad(' ',24); 

   s.catf(0,0," %s",sR);
   s.check_bounds(FL,1);

   hid.init(1,&h0);

   return acquire_set(ss,ntry);
};

template<class TQ>
int CG::Guard::acquire(
   const char *F, int L, const QSet<TQ> &Q, const char *tag,
   const char *cond 
) {
   if (!Wb::omp_parallel()) { return 0; } 

   size_t h0=QHash<TQ>()(Q,tag);

   int ntry=-1;
   if (cond && cond[0]) {
      ntry=deal_with_cond(F,L,h0,cond); if (!ntry) return 0; }

   wbvec<char> s(64);
   const char* ss[1]={s.data};

   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   if (tag && tag[0]) 
        { s.catf(0,0,"%s", shortFL(F_L,-1,  tag,'*')); }
   else { s.catf(0,0,"%s", shortFL(F_L,-1,mytag,'/')); }

   s.pad(' ',24); 
   s.catf(F,L," %s",STR(Q));

   hid.init(1,&h0);

   return acquire_set(ss,ntry);
};

template<class TM>
int CG::Guard::acquire(
   const char *F, int L, const cgc_contract_id<TM> &idc,
   const char *cond 
){
   if (!Wb::omp_parallel()) { return 0; } 

   int ntry=-1;
   size_t h0=MHash<TM>()(idc);

   if (cond && *cond) {
      ntry=deal_with_cond(F,L,h0,cond);
      if (!ntry) { return 0; }
   }

   wbvec<char> s(128);
   const char* ss[1]={s.data};
   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   s.catf(0,0,"%-20s %s", 
      shortFL(F_L,-1,mytag,'/'), 
      STR2(idc,0) 
   );

   hid.init(1,&h0);

   return acquire_set(ss,ntry);
};

template<class TQ>
void CG::Guard::acquire(
   const char *F, int L, const QSet<TQ> *Q[], unsigned n
){
   if (!Wb::omp_parallel()) { return; } 

   unsigned i=0, k=0, slen=64;
   wbvec<char*> ss(n);

   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);
   hid.init(n);

   for (i=0; i<n; ++i) { if (Q[i]) {
      if (!Q[i]->isEmpty()) {
         hid[k]=QHash<TQ>()(*Q[i]);
         ss[k] = new char[slen];  
         snprintf(ss[k++],slen,
           "%-20s %s",shortFL(F_L,-1,mytag,'/'), STR_(Q[i])
         );
      }
      else wblog(FL,"WRN %s() got empty QSet (%d/%d)",FCT,i+1,n);
   }}

   if (!k) wblog(FL,
      "ERR Guard() got no relevant QSet data (%d/%d)",k,n);
   hid.len=k; for (; k<n; ++k) { ss[k]=NULL; }

   acquire_set((const char**)ss.data);

   for (i=0; i<n; ++i) { if (ss[i]) { delete [] ss[i]; }}
};

int CG::Guard::acquire_set_iter(
   int *got, Wb::ompNLock *lk[], const char *ss[], unsigned iter
){

   int nlks=0; 
   unsigned i=0, gotnew;

   { Wb::ompGuard myLK(CG_locks);
     for (; i<hid.len; ++i) {
        Wb::ompNLock *&lp=CG::lock_map[hid[i]]; 
        if ((gotnew=(lp ? 0 : 1))) { 
           lp = new Wb::ompNLock(ss[i]);
        }
        lk[i]=lp;

        if (!iter) { ++(lp->nactive); } 
        got[i] = lp->test_acquire();    

        if (got[i]) { ++nlks; }

        #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
        if (got[i] || gotnew) {
           LKF.blogf(FL, gotnew? "add %s":" +  %s", STR2_(lk[i],ss[i]));
        }
        #endif
   } }

   if ((unsigned)nlks<hid.len) {
      for (i=0; i<hid.len; ++i) { if (got[i]) lk[i]->release(); }}

   return ((unsigned)nlks>=hid.len ? nlks: nlks-hid.len);
};

int CG::Guard::acquire_set_wait(
   int *got, Wb::ompNLock *lk[], const char *ss[], int ntry) {

   int nlks=0; unsigned i, missed=0, ndead=0, nlog=0;
   int itry=1;

   double wt=-1,
      wt1=60,  
      wt2=600; 

   if (ntry<=0) {
      if (!ntry) wblog(FL,"WRN %s() got ntry=%d",FCT,ntry);
      ntry=1200;
   }

   Wb::thread_xlink xl; 

   for (; itry<ntry; ++itry) {
      nlks=acquire_set_iter(got,lk,ss,itry); 

      missed=(nlks<=0 ? 1 : 0);
      if (!missed) break; 

      xl.init(FL);
      for (i=0; i<hid.len; ++i) { if (!got[i]) xl.add(FL,lk[i]); }

      if (wt<=0) {
         wt=0.01*std::rand()/double(RAND_MAX); if (wt<1e-4) wt=1e-4; }
      else {
         double x=1+std::rand()/double(RAND_MAX);
         if ((wt*=x)>wt2) { wt=wt2; } 
      }

      i=Wb::thread_xlink::check_deadlocks();
      ndead=(i? ndead+1 : 0);

      if (wt>wt1 || (i && wt>1)) {  
         Wb::termcolor Tc("WRN");
         wbstring sout; sprintf_locks(sout,got,ss," ");

         fprintf(stdout,"\n%-12s %s "
           "itry %d/%d @ dt=%4.3gs %s(th=%d/%d)%s\n%s",
            shortFLT, itry,ntry,wt, Tc.e1,xl.t0,xl.nt,Tc.em, sout.data);

         char *sx1;
         wbvec<char> sx(64);

         sx.catf(0,0,"itry=%d/%d: %d OMP deadlocks? ",itry,ndead,i);
         sx1=sx.current();
         sx.catf(0,0,"(th=%d/%d)",xl.t0,xl.nt);

         if (ndead) { char xflag=(wt>wt1 && ndead>10);
            if (++nlog==1 || xflag) { Wb::termcolor Tc("WRN");
               fprintf(stdout,
                  "%s %s %sWRN %d OMP deadlocks?%s (%ld CG::locks total)\n",
                   shortFLT, Tc.e1,i,Tc.em, CG::lock_map.size());
               print_CG_locks(FL);
            }

            if (xflag)
                 { wblog(FL,"ERR %s %.3g, %d",sx.data,wt,ndead); }
            else { wblog(FL,"WRN %s",sx.data); }
         }
         else wblog(FL,
            "WRN waiting for OMP lock%ss %s", hid.len==1? "":"s",sx1);
      }
      Wb::pause(wt);
   }

   if ((missed && ntry>99)
     #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
      || (itry>1) 
     #endif
   ) {
      char sx[64]; wbstring sout; sprintf_locks(sout,got,ss);

      snprintf(sx,64,"%ld lock%s (%d %s @ %.3g; %d/%d)", hid.len,
         hid.len==1?"":"s",itry,itry!=1?"tries":"try",wt,xl.t0,xl.nt);

     #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
      LKF.blogf(FL, missed? "ERR %s\n%s":"ok. %s\n%s",sx,sout.data);
     #endif

      if (missed) fprintf(stderr,"%-12s %s "
         "ERR failed to acquire %s\n%s",shortFLT,sx,sout.data);
   }

   return nlks; 
};

void CG::Guard::release() {
   if (!hid.len) { active=false; return; }

   unsigned i;
   wbvec<char> err(hid.len);

 { Wb::ompNLock *lp;
   map <size_t, Wb::ompNLock*>::iterator it;
   Wb::ompGuard myLK(CG_locks);

   for (i=0; i<hid.len; ++i) { err[i]=0;
      it=CG::lock_map.find(hid[i]);
      if (it!=CG::lock_map.end()) { lp=it->second;
         if (lp) {
            if (active) {
               lp->release(); 
            }
            else if (lp->owner && lp->owner==omp_get_thread_num()) {
               err[i]|=1;
            }

            if ((--lp->nactive)<=0) {
             # if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
               LKF.blogf(FL,"rm  %s",STR_(lp));
             # endif

               delete lp; CG::lock_map.erase(it);
            }
         } else err[i]|=2;
      } else err[i]|=4;
   }
 }

   for (i=0; i<hid.len; ++i) { if (err[i]) {
      if (i) { err[0]|=err[i]; }
      if (err[i] & 1) { wblog(FL, 
        "ERR %s() still owner of released lock %d/%d #%lx",
         FCT,i+1,hid.len,hid[i]);
      }
      else if (err[i] & 2) { wblog(FL,
        "ERR %s() got null lock %d/%d #%lx",FCT,i+1,hid.len,hid[i]);
      }
      else if (err[i] & 4) {wblog(FL,
        "ERR %s() lock %d/%d #%lx no longer exists",FCT,i+1,hid.len,hid[i]);
      }
   }}
   if (err[0]) { fflush(0); }

   hid.init();
   active=false;
};

void CG::Guard::print(const char *F, int L) {

   if (!hid.len) {
      wblog(F_L," *  %s() no locks in place (#hid=%d)",FCT,hid.len);
      return;
   }

   if (hid.len==1)
        wblog(F_L," *  CG::Guard() got single lock entry");
   else wblog(F_L," *  CG::Guard() contains %d locks",hid.len,hid.len);

   Wb::ompGuard myLK(CG_locks);
   for (unsigned i=0; i<hid.len; ++i) {
      auto it=CG::lock_map.find(hid[i]);
      wblog(FL," *   %d. %-24s #%lx",i+1,
         it==CG::lock_map.end() ? "(not in locks)" :
         (it->second ?  it->second->istr : "(null)"), hid[i]
      );
   }
};

int Wb::thread_xlink::add(const char *F, int L, size_t h) {
   Wb::ompGuard myLK(CG_locks);
   const auto im = CG::lock_map.find(h);
   if (im==CG::lock_map.end()) wblog(FL,
      "ERR %s() non existing lock CG::lock_map[%x]",FCT,h);
   return add(F_L,im->second);
};

#endif 
namespace Wb {

int thread_xlink::add(const char *F, int L, const ompNLock *lk) {

   unsigned t2=0, n=8*sizeof(tx); size_t b=1;

   assert(lk!=NULL);

   t2=lk->owner;
   if (int(t2)<0) { 
      return (lk->level<=0 ? -1 : 0);
   }

   if (t0>=n || t2>=n) wblog(F_L,
      "ERR %s() thread ID out of bounds %d/%d/%d/%d\n%s",
      FCT,t0,t2,nt,n,STR_(lk));

   if (tx & (b<<=t2)) { return 0; } 
   else { tx|=b; return 1; }
};

int check_deadlock_rec(unsigned i0, size_t tt, const thread_xlink *x) {

   unsigned r=0, i=0; size_t b=1, q=x->tx;

   for (; q && !r; ++i, b<<=1, q>>=1) { if (q & 1) { 
      if (tt & b) {
         if (i==i0) ++r;
      }
      else {
         const auto it=txmap.find(i);
         if (it!=txmap.end()) { 
            r+=check_deadlock_rec(i0, tt | b, it->second);
         }
      }
   }}

   return r;
};

int thread_xlink::check_deadlocks() { 

   unsigned i, r=0, n=8*sizeof(size_t);

  #pragma omp critical (thread_xlink__)
   for (auto it=txmap.begin(); it!=txmap.end(); ++it) {
      if ((i=it->first)<n)
         r+=check_deadlock_rec(i, (size_t(1))<<i, it->second);
      else wblog(FL,
        "WRN %s() thread ID out of bounds (%d/%d)",FCT,i,n);
   }

   return r;
};

}; 

#endif 

#endif 

