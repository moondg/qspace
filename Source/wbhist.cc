
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to wbhist.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "wbhist"
#endif

#include "wblib.h"

class wbhist { 
 public:

   wbhist() : ndat(0), ncalls(0) {};
   wbhist(const mxArray *a) : ndat(0), ncalls(0) { init(a); }

   wbhist& init(const mxArray *a) {
      unsigned i,n;

      if (!a || !mxIsWbvector(FL,a)) wblog(FL,
         "ERR %s() invalid setup up for binning",PROG);

      xbin.init(FL,a); xbin.Sort(); n=xbin.len; ncalls=ndat=0;
      ybin.init(n); 

      if (n<2) { xbin.init(); ybin.init();
         wblog(FL,"ERR %s() invalid xbin (%d)",PROG,xbin.len);
      }

      xbdr.init(--n);
      double *x=xbdr.data, *x0=xbin.data;
      for (i=0; i<n; ++i) { x[i]=0.5*(x0[i] + x0[i+1]); }

      return *this;
   };

   wbhist& add(const char *F, int L, const mxArray *ax, const mxArray *ay){

      wbvector<double> X,Y; wbperm P;
      X.init(F_L,ax); Y.init(F_L,ay);

      if (xbin.len!=ybin.len) wblog(F_L,"ERR %s() "
         "invalid binning data (%d/%d)",PROG,xbin.len,ybin.len);
      if (Y.len && X.len!=Y.len) wblog(F_L,"ERR %s() "
         "invalid input data (length=%d/%d)",PROG,X.len,Y.len);

      X.Sort(P); if (Y.len) Y.Permute(P);

      unsigned i,n=xbdr.len, ib=0; bool goty=(Y.len!=0);
      const double *xl=xbdr.data, *x=X.data, *y=Y.data;
      double *yb=ybin.data;

      for (i=0; i<X.len; ++i) {
         while (xl[ib]<x[i] && ib<n) { ++ib; }
         yb[ib]+=(goty ? y[i] : 1);
      }

      ++ncalls; ndat+=n;

      return *this;
   };

   long unsigned ndat, ncalls;
   wbvector<double> xbdr, xbin, ybin;

 protected:
 private:

};

   map<string, wbhist> gHist;

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin,  const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   unsigned k=0,i=0; string id("default");

   MX_CHECK_HELPER_NARGS(0,-1,-1);

   while (nargin && mxIsChar(argin[0]) && k<2) {
      int i=mxGetString(argin[0],str,32);
      if (i) wblog(FL,"ERR %s() failed to read string (%d)",PROG,k+1);

      if (!strcmp(str,"remove")) {
         if (nargin>1 || nargout) wblog(FL,
            "ERR %s('%s') invalid usage",PROG,id.data());
         gHist.erase(id.data()); return;
      }
      else if (!strcmp(str,"clear")) {
         if (nargin>1 || nargout) wblog(FL,
            "ERR %s('%s') invalid usage",PROG,id.data());
         gHist.clear(); return;
      }
      else if (!strcmp(str,"stat")) {
         if (nargin>1 || nargout) wblog(FL,
            "ERR %s('%s') invalid usage",PROG,id.data());

         if (gHist.size()<1) {
            printf("\n   wbhist status: (empty)\n\n");
            return;
         }

         printf("\n   wbhist status:\n\n");

         for (map<string,wbhist>::iterator
             im=gHist.begin(); im!=gHist.end(); ++im
         ){
             const wbhist &H=im->second;
             if (H.xbin.len) {
                printf("     %-12s: %.3g .. %.3g (%d; {%ld@%ld})\n",
                im->first.data(), H.xbin.data[0],H.xbin.data[H.xbin.len-1],
                H.xbin.len,H.ncalls,H.ndat);
             }
             else { printf("     %-12s: (empty)\n",im->first.data()); }
         }
         printf("\n"); return;
      }
      else if (k==0) { id=str; }
      else wblog(FL,"ERR %s() invalid task '%s'",PROG,str);

      ++k; ++argin; --nargin;
   }

   wbhist &H=gHist[id];

   if (nargin) {

      if (nargin<2 || nargin>3 || nargout>1)
      usage(FL,"invalid number of I/O arguments");

      for (i=0; int(i)<nargin; ++i) {
         if (argin[i] && mxGetNumberOfElements(argin[i])) {
         if (!mxIsWbvector(FL,argin[i])) wblog(FL,
            "ERR %s('%s') invalid usage (expecting vector as arg #%d)",
             PROG,id.data(),k+i+1);
         }
      }

      if (nargin==3) { H.init(argin[2]); }
      if (!H.xbin.len) wblog(FL,"ERR %s('%s') "
         "got empty xbin (not initialized yet!?)",PROG,id.data());

      H.add(FL,argin[0],argin[1]);
   }

   if (nargout) {
      if (nargout==1) { argout[0]=H.ybin.toMx(); return; } else
      if (nargout==2) {
         argout[0]=H.xbin.toMx();
         argout[1]=H.ybin.toMx(); return;
      }
      else if (nargout==3) {
         argout[0]=H.xbin.toMx();
         argout[1]=H.ybin.toMx();
         argout[2]=H.xbdr.toMx(); return;
      }
   }
   return;

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbhist"); }
   aclu.Check();
};

