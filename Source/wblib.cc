#ifndef __WB_WBLIB_CC__
#define __WB_WBLIB_CC__

/* -------------------------------------------------------------------- */

Wb::gpara__::gpara__() {
#ifdef __WBDEBUG__
// #WBLIB_SRAND_INIT
   wblog(FL,"TST %s() setting up global %s",myname,FCT);
#endif
// ensure first-time initialization with init()! // Wb,Aug28,25
   tlast=0; // HEADER_INIT_DT

   wb_srand(); 

#ifdef __WB_MPFR_HH__
   gmp_randinit_default(Wb::wb_rstate);
#endif

   envDKT=got_DESKTOP();
   my_caller_tid=omp_get_thread_num(); 

   init();
};

void Wb::gpara__::init(char force) {
   time_t tnow=time(NULL), dt=tnow-tlast; 

   #ifdef __WBDEBUG__
      wblog(FL,"TST gpara__::%s() %s ENV (dt=%ld)", FCT,
      dt>=1 ? "checking":"skipping", dt);
   #endif

   wbl::status=0; 

   envVRB=get_WB_VERBOSE();  

   envDBG=got_DBSTOP();

   GetEnvInt(0,0,"QS_FERM",&envFERM,"envFERM"); 

   GetEnvInt(0,1,"QS_FULL_OM",&envFullOM,"envFullOM");

   if (!force && dt<1) { return; }

   memset(str,0,STRLEN+1);   

   #ifdef QS_USING_OMP
      Wb::GetNumThreads(FL_q(tlast),OMP_NUM_THREADS,"OMP_NUM_THREADS");
      Wb::GetNumThreads(FL_q(tlast),QSP_NUM_THREADS,"QSP_NUM_THREADS");

      sp_num_threads=MAX(OMP_NUM_THREADS,QSP_NUM_THREADS);

    # ifdef __APPLE__
      static unsigned ncall=0;
      if (dt<(1<<20)) { 
      if (QSP_NUM_THREADS>1 && OMP_NUM_THREADS!=1 && ++ncall<=1) {
         PRINTF("\n%s \e[31m ERR nested parallelization\n\n"
            "   NB! macOS may not support nested parallelization\n"
            "   got QSP_NUM_THREADS=%d with OMP_NUM_THREADS=%d%s\n"
            "   Hint: use setNumThreads(1) or (un)set QSP_NUM_THREADS<=1\n"
            "\e[0m\n", PROG, QSP_NUM_THREADS,
            OMP_NUM_THREADS, OMP_NUM_THREADS ? "":"(=auto)");
         ExitMsg("no nested parallelization on macOS (system may hang)");
      }}
    # endif
   #endif

   tlast=tnow;

   my_caller_tid=omp_get_thread_num(); 
};

void Wb::gpara__::info(const char *F, int L) const {

   struct tm *tblock=localtime(&tlast);
   char sx[32]; strftime(sx,32,"%D %T",tblock);

   wblog(F_L,"<i> gpara::%s()%N",FCT);
   PRINTF("  %-16s %s\n","tlast",sx);
   PRINTF("  %-16s %d\n","wbl::status",wbl::status);
   PRINTF("  %-16s %d\n","envFERM",envFERM);
   PRINTF("  %-16s %d\n","envVRB",envVRB);
   PRINTF("  %-16s %d\n","envDBG",envDBG);

 # ifdef QS_USING_OMP
   PRINTF("  %-16s %d\n","OMP_NUM_THREADS",OMP_NUM_THREADS);
   PRINTF("  %-16s %d\n","QSP_NUM_THREADS",QSP_NUM_THREADS);
   PRINTF("  %-16s %d\n","sp_num_threads",sp_num_threads);
 # endif
   PRINTF("\n");
};

#ifdef __APPLE__
int Wb::system_tid() {
   uint64_t tid=0;
   pthread_threadid_np(nullptr, &tid);
   return tid;
}
#else
int Wb::system_tid() { 
   return syscall(SYS_gettid);
}
#endif

void Wb::CleanUp::Check() {
   wbl::check_ERR_pending();
};

Wb::CleanUp::~CleanUp() {
   wbl::myIO.clear();
   if (gwb_Profs.size()) { Wb::save_and_clear_Profiling(); }
};

template<class T>
Wb::tmpSet__<T>::~tmpSet__() { if (v) {
   *v=x_done; v=NULL;
}};

template<class T>
Wb::tmpSet__<T>& Wb::tmpSet__<T>::set(T x, char op) {
   if (!v) wblog(FL,"ERR %s() not yet initialized",FCT);
   switch (op) {
      case '=': *v  = x; break;
      case '|': *v |= x; break;
      case '&': *v &= x; break;
      case '+': *v += x; break;
      default: wblog(FL,"ERR %s() invalid switch %s",FCT,cSTR(op));
   }
   return *this;
};

template<class T>
Wb::iterLevel<T>::iterLevel(T *v_, T dx_)
 : v(v_), dx(dx_), file(NULL), line(0) { (*v)+=dx; };

template<class T>
Wb::iterLevel<T>::iterLevel(
   const char *F, int L, const char *fct, T *v_, const char *istr, T dx_)
 : v(v_), dx(dx_), file(NULL), line(0) {

   if (F && *F) { unsigned n=strlen(F), i=n-1;
      for (; i<n; --i) { if (F[i]=='/') { break; }}
      if (++i<n) {
         file = new char[n-i+1]; strcpy(file,F+i);
         line=L;
         ::wblogf(stdout,F,L,"TST ┌ %s() %s", fct?fct:FCT, istr?istr:"");
      }
   }; (*v)+=dx;
};

template<class T>
Wb::iterLevel<T>::~iterLevel() {
   if (file) {
      if ((void*)v==(void*)&wbl::level) {
         (*v)-=dx; v=NULL; 
         wblog(file,line,"TST └");
      }
      else if (v) {
         wbvec<char> sx(32);
         if (*v>=256 && T(int(*v))==*v) 
              { int b=(*v)/256, a=(*v)-b*256; sx.catf(0,0,"(%d,%d)",a,b); }
         else { sx.catf(0,0,"%g",double(*v)); }
         (*v)-=dx;
         if (*v>=256 && T(int(*v))==*v)
              { int b=(*v)/256, a=(*v)-b*256; sx.catf(0,0," -> (%d,%d)",a,b); }
         else { sx.catf(0,0," -> %g",double(*v)); }
         wblog(file,line,"WRN %s() %s",FCT,sx.data);
         v=NULL;
      }
      delete [] file; file=NULL;
   }

   if (v) { (*v)-=dx; v=NULL; }
   line=0; 
};

#endif

