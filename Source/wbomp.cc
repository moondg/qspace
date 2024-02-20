/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace OMP routines
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

#ifndef __WB_OMP_CC__
#define __WB_OMP_CC__

#ifdef QS_USING_OMP

wbstring Wb::ompID2Str(char vflag) {
   wbstring s(8); 
   unsigned k=(vflag ? 16 : 0),
   i=omp_get_thread_num(), n=omp_get_num_threads(), l=omp_get_level();

   if (!i && n==1 && !l) { s.init(8+k); strcpy(s.data,"serial"); l=6; }
   else if (l!=1) 
        { s.init(16+k); l=snprintf(s.data,s.len,"OMP-%d/%d@%d",i,n,l); }
   else { s.init(12+k); l=snprintf(s.data,s.len,"OMP-%d/%d",   i,n  ); }

   if (k && l<s.len) {
      l+=snprintf(s.data+l, s.len-l,", %s",Wb::hostid(vflag).data);
   }

   return s;
};

int Wb::ompStatus(const char *F, int L) {
    unsigned l=0, n=1024;
    char fmt[]="  %-30s %2d\n", s[n];

    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_active_level", omp_get_active_level()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_cancellation", omp_get_cancellation());  
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_default_device", omp_get_default_device()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_dynamic", omp_get_dynamic());   
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_level", omp_get_level());       
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_max_task_priority", omp_get_max_task_priority()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_max_threads", omp_get_max_threads());  
    if (l<n) l+=snprintf(s+l,n-l,fmt,   
      "omp_get_max_active_levels", omp_get_max_active_levels());
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_num_procs", omp_get_num_procs()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_num_teams", omp_get_num_teams()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_num_threads", omp_get_num_threads());  

    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_proc_bind", omp_get_proc_bind());
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_team_num", omp_get_team_num());
    if (l<n) l+=snprintf(s+l,n-l,"\n"); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_get_thread_num", omp_get_thread_num()); 
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_in_parallel", omp_in_parallel());
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_in_final", omp_in_final());
    if (l<n) l+=snprintf(s+l,n-l,fmt,
      "omp_is_initial_device", omp_is_initial_device());

    fprintf(stdout,"\n%-18s %s -> %s()\n\n",shortFL(F_L),myname,FCT);
    if (l<n)
         fprintf(stdout,"%s\n",s);
    else fprintf(stdout,"%s  ... (%d/%d)\n",s,l,n);

    return 0;
};

int Wb::ompNLock::acquire() {

   omp_set_nest_lock(&_nlk);
   owner=omp_get_thread_num(); ++ntot; ++level;

#if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 12 ) 
   LKF.blogf(FL," +  %s%s",STR_(this), check_BUF_size().data);
#endif

   return level;
};

int Wb::ompNLock::test_acquire() {
   int l=omp_test_nest_lock(&_nlk);
   if (!l) return l;

   owner=omp_get_thread_num(); ++ntot; ++level;

   return level;
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

   if (level<=0) wblog(FL,"WRN %s() %s",FCT,STR_(this));

   int l=(--level); 
   if (l<=0) owner=-1; 

   #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 12 ) 
    LKF.blogf(FL," -  %s%s",STR_(this),check_BUF_size().data);
   #endif

   omp_unset_nest_lock(&_nlk);

   return (l>0? l : 0);
};

#ifdef LOAD_CGC_QSPACE

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

   if ( nactive && !ntot ) l+=snprintf(s,n,"%d/%d", level,nactive); else
   if (!nactive && ntot>9) l+=snprintf(s,n,"%d/%03ld", level,ntot );
   else                    l+=snprintf(s,n,"%d/%d/%ld",level,nactive,ntot);

   l+=Wb::strpad_(s+l,8-l);
   l+=snprintf(s+l,n-l," %2d",owner); 

   if (l<n) {
      if (istr && istr[0]) l+=snprintf(s+l, n-l," %s", istr );
      else if (istr_)      l+=snprintf(s+l, n-l,"`%s'",istr_);
      else {
         l+=snprintf(s+l, n-l," (null)");
         wblog(FL,"WRN %s() got null istr",FCT);
      }
   }
   if (l<n && (!istr || !istr[0])) {
      l+=snprintf(s+l, n-l," %p",this); 
   }
   return sout;
};

void Wb::ompNLock::set_istr(const char *s) {

   if ((!istr || !istr[0]) && s && s[0]) {
      if (istr) { delete [] istr; }
      istr = new char[strlen(s)+1];
      strcpy(istr,s);
   }
};

void Wb::ompGuard::acquire(const char *F, int L, ompNLock &lk, int lmax) {

   if (lp) wblog(FL,"ERR %s() already initialized with lock",FCT);

   if (lk.level<lmax) { lp=&lk;
      lk.acquire();

      if (lk.ntot<=1 && F) {
         Wb::termcolor Tc("WRN");
         fprintf(stdout,"\n  %sTST %s %s() "
           "using ompCriticial `%s'%s\n\n",Tc.e1,shortFL(F_L),FCT,
            lk.istr? lk.istr:"(unnamed ompNLock)",Tc.em);
         fflush(0);
      }
   }
   else if (lk.level>lmax) fprintf(stdout,"%s %s "
   "WRN got lk.level=%d/%d\n",shortFLT,lk.level,lmax);
};

#ifdef LOAD_CGC_QSPACE

namespace CG {

   Wb::ompNLock* Guard::find(size_t h) {
      auto it = CG::lock_map.find(h);
      if (it==CG::lock_map.end())
           { return NULL; }
      else { return it->second; }
   };

   template<class TQ>
   wbstring Guard::Status(const QSet<TQ> &Q, const char *tag) {
      size_t h1=QHash<TQ>()(Q,tag); wbstring istr=Q.toStr();
      if (tag && tag[0]) {
         unsigned l=strlen(tag)+4; char sx[l];
         snprintf(sx,l," '%s'",tag); istr+=sx;
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
      else { s.init(32);
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
         if (i<n && k<hid.len) { s[i++]='\n'; }
      }
      s[i<n ? i:n-1]=0;
      if (i>=n) wblog(FL,"ERR string out of bounds (%d/%d)\n%s",i,n,s);
      return sout;
   };
};

int CG::Guard::print_CG_locks() { 
   unsigned i=0;

   fprintf(stdout,"level/nactive/nuse, owner_thread    info_string\n");
   for (auto I=CG::lock_map.begin(); I!=CG::lock_map.end(); ++I) {
      fprintf(stdout,"%2d. %s\n",++i,STR_(I->second));
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
      { Wb::ompGuard gLK(CG_locks);
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
   if (omp_get_num_threads()<=1) { return 0; }

   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   size_t h0; unsigned l,l2,n=64;
   int ntry=-1; const char *sR; char s[n]; const char* ss[1]={s};

   flk.ensureFOpen(FL);

   sR=Wb::rstrstr(flk.fname.data,"/RCS"); 
   if (sR) { ++sR; }
   else { const char *s0=flk.fname.data;
      sR=Wb::rstrstr(s0,"Store");
      if (sR>s0) { --sR; }
      else { sR=s0; } 
   }
   l2=strlen(sR);

   h0 = Wb::string_hash(tag ? tag : mytag)
      ^ Wb::string_hash(sR);

   if (cond && cond[0]) {
      ntry=deal_with_cond(F,L,h0,cond);
      if (!ntry) { return 0; }
   }

   if (tag && tag[0]) 
        { l=snprintf(s,n,"%s", shortFL(F_L,-1,  tag,'*')); }
   else { l=snprintf(s,n,"%s", shortFL(F_L,-1,mytag,'/')); }
   for (; l<24; ++l) { s[l]=' '; } 

   if (l+l2+1>=n) wblog(FL,
      "ERR %s() string out of bounds (%d+%d/%d)\n%s %s",FCT,l,l2,n,s,sR);
   l+=snprintf(s+l,n-l," %s",sR);

   hid.init(1,&h0);

   return acquire_set(ss,ntry);
};

template<class TQ>
int CG::Guard::acquire(
   const char *F, int L, const QSet<TQ> &Q, const char *tag,
   const char *cond 
) {
   if (omp_get_num_threads()<=1) { return 0; }

   size_t h0=QHash<TQ>()(Q,tag);

   int ntry=-1;
   if (cond && cond[0]) {
      ntry=deal_with_cond(F,L,h0,cond); if (!ntry) return 0; }

   unsigned l, n=64; char s[n]; const char* ss[1]={s};
   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   if (tag && tag[0]) 
        { l=snprintf(s,n,"%s", shortFL(F_L,-1,  tag,'*')); }
   else { l=snprintf(s,n,"%s", shortFL(F_L,-1,mytag,'/')); }

   if (l<n) {
      for (; l<24; ++l) { s[l]=' '; } 
      l+=snprintf(s+l,n-l," %s",STR(Q));
   }
   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,s);

   hid.init(1,&h0);

   return acquire_set(ss,ntry);
};

template<class TM>
int CG::Guard::acquire(
   const char *F, int L, const cgc_contract_id<TM> &idc,
   const char *cond 
){
   if (omp_get_num_threads()<=1) return 0; 

   int ntry=-1; size_t h0=MHash<TM>()(idc);
   if (cond && cond[0]) {
      ntry=deal_with_cond(F,L,h0,cond); if (!ntry) return 0;
   }

   unsigned n=128; char s[n]; const char* ss[1]={s};
   if (hid.len) wblog(FL,
      "ERR Guard() already initialized (len=%d)",hid.len);

   snprintf(s,n,"%-20s %s", 
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
   if (omp_get_num_threads()<=1) return; 

   unsigned i=0, k=0, slen=64; char* ss[n];
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

   acquire_set((const char**)ss);

   for (i=0; i<n; ++i) { if (ss[i]) { delete [] ss[i]; }}
};

int CG::Guard::acquire_set_iter(
   int *got, Wb::ompNLock *lk[], const char *ss[], unsigned iter
){

   int nlks=0; 
   unsigned i=0, gotnew;

   { Wb::ompGuard gLK(CG_locks);
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
   int iter=1;

   double wt=-1,
      wt1=60,  
      wt2=600; 

   if (ntry<=0) {
      if (!ntry) wblog(FL,"WRN %s() got ntry=%d !?",FCT,ntry);
      ntry=1200;
   }

   Wb::thread_xlink xl; 

   for (; iter<ntry; ++iter) {
      nlks=acquire_set_iter(got,lk,ss,iter); 

      missed=(nlks<=0 ? 1 : 0);
      if (!missed) break; 

      xl.init(FL);
      for (i=0; i<hid.len; ++i) { if (!got[i]) xl.add(FL,lk[i]); }

      if (wt<=0) {
         wt=0.01*std::rand()/double(RAND_MAX); if (wt<1E-4) wt=1E-4; }
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
           "itry %d/%d @ dt=%4.3gs %s(ith=%d/%d)%s\n%s",
            shortFLT, iter,ntry,wt, Tc.e1,xl.t0,xl.nt,Tc.em, sout.data);

         unsigned l, n=64; char sx[n];
         l=snprintf(sx,n,"%d/%d OMP deadlocks? ",i,ndead);
         snprintf(sx+l,n-l,"(ith=%d/%d)",xl.t0,xl.nt); 

         if (ndead) { char xflag=(wt>wt1 && ndead>10);
            if (++nlog==1 || xflag) {
               Wb::termcolor Tc("WRN");
               fprintf(stdout,
                  "%s %s %sWRN %d deadlocks?%s (%ld locks total)\n\n",
                   shortFLT, Tc.e1,i,Tc.em, CG::lock_map.size());
               print_CG_locks();
            }

            if (xflag)
                 { wblog(FL,"ERR %s %.3g, %d",sx,wt,ndead); }
            else { wblog(FL,"WRN %s",sx); }
         }
         else wblog(FL,
            "WRN waiting for OMP lock%ss %s", hid.len==1? "":"s",sx+l);
      }
      Wb::pause(wt);
   }

   if ((missed && ntry>99)
     #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
      || (iter>1) 
     #endif
   ) {
      char sx[64]; wbstring sout; sprintf_locks(sout,got,ss);

      snprintf(sx,64,"%ld lock%s (%d %s @ %.3g; %d/%d)", hid.len,
         hid.len==1?"":"s",iter,iter!=1?"tries":"try",wt,xl.t0,xl.nt);

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

   unsigned i=0; char err[hid.len];

   map <size_t, Wb::ompNLock*>::iterator it;
   Wb::ompNLock *lp;
 { Wb::ompGuard gLK(CG_locks);

   for (i=0; i<hid.len; ++i) { err[i]=0;
      it=CG::lock_map.find(hid[i]);
      if (it!=CG::lock_map.end()) { lp=it->second;
         if (lp) {
            if (active) { lp->release(); } 
            else if (lp->owner==omp_get_thread_num()) { err[i]|=1; }

            if ((--lp->nactive)<=0) {
               #if defined(DBG_GCX_LOCKS) && ( DBG_GCX_LOCKS & 8 )
                LKF.blogf(FL,"rm  %s",STR_(lp));
               #endif

               delete lp; CG::lock_map.erase(it);
            }
         } else err[i]|=2;
      } else err[i]|=4;
  }}

   for (i=0; i<hid.len; ++i) { if (err[i]) { if (i) { err[0]|=err[i]; }
      if (err[i] & 1) wblog(FL,"ERR %s() " 
         "still owner of released lock %d/%d #%lx",FCT,i+1,hid.len,hid[i]);
      else if (err[i] & 2) wblog(FL,"ERR %s() "
         "got null lock %d/%d #%lx",FCT,i+1,hid.len,hid[i]);
      else if (err[i] & 4) wblog(FL,"ERR %s() "
         "lock %d/%d #%lx no longer exists !?",FCT,i+1,hid.len,hid[i]);
   }}; if (err[0]) fflush(0);

   hid.init(); active=false;
};

void CG::Guard::print(const char *F, int L) {

   if (!hid.len) {
      wblog(F_L," *  %s() no locks in place (#hid=%d)",FCT,hid.len);
      return;
   }

   if (hid.len==1)
        wblog(F_L," *  CG::Guard() got single lock entry");
   else wblog(F_L," *  CG::Guard() contains %d locks",hid.len,hid.len);

   Wb::ompGuard gLK(CG_locks);
   for (unsigned i=0; i<hid.len; ++i) {
      auto it=CG::lock_map.find(hid[i]);
      wblog(FL," *   %d. %-24s #%lx",i+1,
         it==CG::lock_map.end() ? "(not in locks)" :
         (it->second ?  it->second->istr : "(null)"), hid[i]
      );
   }
};

int Wb::thread_xlink::add(const char *F, int L, size_t h) {
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

