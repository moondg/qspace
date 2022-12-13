
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to uniquerows.m // Wb,Feb14,19

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   WBINDEX D;
   wbperm P; mxArray *a=0;
   char uflag=0; 

   MX_CHECK_HELPER_NARGS(1,-1,3); 
   if (nargin>1) {
      wbstring s(argin[1]);
      if (s=="-1") { uflag=1; --nargin; }
   }
   if (nargin!=1) usage(FL,"Invalid number of input arguments");

   if (mxIsDouble(argin[0])) {
      if (mxIsComplex(argin[0])) { 
         wbMatrix<wbcomplex> Q; Q.init(argin[0]);
         Q.groupRecs(P,D);
         a=Q.toMx();
      }
      else {
         wbMatrix<double> Q; Q.init(argin[0]);
         Q.groupRecs(P,D);
         a=Q.toMx();
      }
   }
   else if (mxIsChar(argin[0])) {
      wbMatrix<char> Q; Q.init(argin[0]);
      Q.groupRecs(P,D);
      a=Q.toMx();
   }
   else {
      wblog(FL,"ERR invalid type of input argument #1 (%s)",
      mxGetClassName(argin[0]));
   }

   argout[0]=a; 

   if (nargout>1) {
      if (uflag) {
         unsigned i,l,d; wbindex I(D.len);
         for (l=i=0; i<D.len; i++, l+=d) { d=D[i]; I[i]=P[l]; }
         argout[1]=I.toMx();
      }
      else {
         unsigned i,l,d; WBINDEX I;
         mxArray *S = mxCreateCellMatrix(1,D.len);

         for (l=i=0; i<D.len; i++, l+=d) {
            d=D[i]; I.init2ref(d, P.data+l); I+=1; 
            mxSetCell(S,i,I.toMx());
         }
         argout[1]=S;
      }
   }

   if (nargout>2) argout[2]=D.toMx();

   return;

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in uniquerows"); }
};

