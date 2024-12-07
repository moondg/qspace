
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to wbhash.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "wbhash"
#endif

#define LD_CLEBSCH_QS
#define QS_SKIP_MPFR

#include "wblib.h"

#include <limits>
#include <unordered_map>

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "wbhash"
#endif

class QHash_ {
   public:
      size_t operator()(const qset<double> &x) const {
         unsigned i=0; long q; unsigned long l, h=5381U;
         for (; i<x.len; ++i) { q=long(x[i]);
            l=(unsigned long)(q<0 ? -1-q : q);
            h ^= ((h<<6) + (h>>2)) + l;
            h ^= ((h<<6) + (h>>2)) + l;
         }
         return h;
      };
};

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin,  const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned i,j,n,m;
   unsigned long h, N, offset=5381, mult=33;
   char type=0;

   MX_CHECK_HELPER_NARGS(2,-1,2); 
   if (!mxIsDouble(argin[0])) usage(FL,
      "invalid usage (numeric array required as first argument)");

   mxGetNumber(argin[1],N); {
      if (N<2 || N>(1<<30))
      wblog(FL,"ERR %s() invalid N=%d",FCT,N);
   }

   if (nargin>2) {
      OPTS opts(argin+2,nargin-2);
      if (opts.getOpt("-stl")) { type='s'; } else
      if (opts.getOpt("-wbb")) { type='b'; } else
      if (opts.getOpt("-ump")) { type='u'; }
      else {
         opts.getOpt("mult",mult);
         opts.getOpt("offset",offset);
      }
      opts.checkAnyLeft();
   }

   wbarray<double> dd(FL,argin[0],'r');
   const double *x=dd.data;

   if (dd.SIZE.len!=2) wblog(FL,
      "ERR %s() expecting (integer) matrix as first argument",FCT);
   for (n=dd.numel(), i=0; i<n; ++i) {
      if (double(char(x[i]))!=x[i]) wblog(FL,
     "ERR %s() expecting char array with STL flag (%g)",FCT,x[i]);
   }

   n=dd.SIZE[0]; m=dd.SIZE[1];

   wbvector<double> hk(n);
   wbstring s_(m+1); char *s=s_.data;

   if (type=='s') { 
      hash<std::string> string_hash;

      for (i=0; i<n; ++i) {
      for (j=0; j<m; ++j) {
         char c = char(x[i+j*n]); 
         s[j]=(c>0 ? c : -1-c);   
      }; hk[i]=(string_hash(s)%N); }

      argout[0]=hk.toMx('t');     
   }
   else if (type=='b') { 
      unsigned long l; long q;
      for (i=0; i<n; ++i) { h=offset; 
      for (j=0; j<m; ++j) { q=long(x[i+j*n]); 
         l=(unsigned long)(q<0 ? -1-q : q);
         h ^= ((h<<6) + (h>>2)) + l;
         h ^= ((h<<6) + (h>>2)) + l; 
      }; hk[i]=(h%N); }

      argout[0]=hk.toMx('t');     
   }
   else if (type=='u') { 

      unordered_map< qset<double>, wbvector<double>, QHash_ > X;
      qset<double> q(m);

      if (!nargout)
      wblog(FL,"<i> max_load_factor=%g",X.max_load_factor());

      for (i=0; i<n; ++i) {
         for (j=0; j<m; ++j) { q[j]=x[i+j*n]; } 
         if (!nargout) {
            if ((j=X.bucket_size(X.bucket(q)))) wblog(FL,
            " -> got collision in bucket #%4d / %4d @ %d",
            X.bucket(q),X.bucket_count(),j+1);
         }
         X[q].init(q);
      }

      unsigned nb=X.bucket_count();

      if (!nargout) {
      wblog(FL,"<i> got %d buckets (%d;%d)",nb,n,N); }
      if (nb<N) {
         X.rehash(N); nb=X.bucket_count(); if (!nargout) {
         wblog(FL," => got %d buckets",nb);  }
      }

      if (!nargout) {
         wblog(FL," *  max_bucket_count=%d",X.max_bucket_count());
         wblog(FL," *  max_size=%d",X.max_size());
      }

      wblog(FL," *  load factor %d/%d = %.3f = %.3f",
      X.size(), nb, X.size()/double(nb), X.load_factor());

      n=X.bucket_count(); m=0;
      wbvector<unsigned> nc(8);
      for (i=0; i<n; ++i) { j=X.bucket_size(i);
         if (m<j) m=j;
         if (j<nc.len) nc[j]++;
      }
      if (m>2)
      wblog(FL,"<i> collisions [%s] <= %d",nc.toStr().data,m);

      i=0;
      for (auto it=X.begin(); it!=X.end(); ++it, ++i) {
         hk[i]=(QHash_()(it->first)%nb);

         if (!nargout) printf(
          " %3d. %8d/%4d [%20s] => [%20s]\n", i+1,
            long(hk[i]), X.bucket(it->first),
            it->first.toStr().data, it->second.toStr().data
         );
      }

      argout[0]=hk.toMx('t');     
   }
   else {

      if (mult==33U) {
         for (i=0; i<n; ++i) { h=offset;
         for (j=0; j<m; ++j) {
            h = ((h<<5) + h) + (unsigned long)x[i+j*n]; 
         }; hk[i]=(h%N); }
      }
      else {
         for (i=0; i<n; ++i) { h=offset;
         for (j=0; j<m; ++j) {
            h = mult*h + (unsigned long)x[i+j*n]; 
         }; hk[i]=(h%N); }
      }

      argout[0]=hk.toMx('t');     
      if (nargout>1) { argout[1]=
         MXPut().add(offset,"offset").add(mult,"mult").add(N,"N").toMx();
      }
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbhash"); }
   aclu.Check();
};

