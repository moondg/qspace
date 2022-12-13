
/* CHANGE LOG / COMMENTS ============================================= *
   see also $MLIB/filestat_old.m
   using [i,s]=system('stat ...') which however, turned out unreliable!
   since the returned string s sometimes can occur with the next(!?!)
   call to system() )#))@@ // Wb,Jan28,16
   NB! this is also BY FAR faster than the system calls above!
 * =================================================================== */

char USAGE[]=""; // outsourced to filestat.m // Wb,Feb14,19

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "filestat"
#endif

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(2,3,1); 

   wbstring fname(argin[0]), w(argin[1]); char sflag=0;
   if (nargin==3) {
      if (!mxIsChar(argin[2]) || !Mx::IsEqual(argin[2],"-s")) wblog(FL,
         "ERR %s() invalid usage (arg #3)",PROG);
      sflag=1;
   }

   if (strlen(w.data)==1) { char q=tolower(w[0]); if (q=='t') { q='m'; }
   if (q=='m' || q=='c' || q=='a') {

      double t=Wb::getFileTime(fname.data,q);

      if (!t) wblog(FL,"ERR %s() "
         "invalid file name !?\n'%s' (%s)",PROG,fname.data,w.data);

      if (!sflag) {
         t=(t+62167309200.00)/double(24*3600);

         argout[0]=numtoMx(t); 
      }
      else {
         time_t it=(time_t)t;
         struct tm *S=localtime(&it);
         unsigned l=0, n=32; wbstring s_(n); char *s=s_.data;
         double dt=t-floor(t);

         l=strftime(s,n,"%d-%b-%Y %H:%M:%S",S); 

         if (l>2 && s[l-1]==10) { s[--l]=0; } 
         if (l<n && dt) {
            char c=0; if (l) { c=(s[--l]-'0'); }
            l+=snprintf(s+l,n-l,"%.5f",c+dt);
         }
         if (l>=n) wblog(FL,
            "ERR %s() string out of bounts (%d/%d)\n%s",PROG,l,n,s);
         argout[0]=s_.toMx(); 
      }

      return;
   }}

   wblog(FL,"ERR %s() got invalid switch '%s' (arg #2)",PROG,w.data);

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in filestat"); }
};

