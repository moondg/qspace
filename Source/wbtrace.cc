
/* COMMENTS / CHANGE LOG ============================================= *
 * =================================================================== */

char USAGE[]=""; // outsourced to wbtrace.m // Wb,Feb14,19

#include "wblib.h"

template<class T>
mxArray* MEX_FUNCTION(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[],
    wbarray<T> &A,
    wbarray<T> &Aout
);

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; try { 

   MX_CHECK_HELPER_NARGS(2,-1,1); 
   if (nargin>3) usage(FL,"ERR invalid number of I/O arguments");

   if (Mx::IsNumArray(0,1,argin[0],-1,'c')<=0) wblog(FL, 
      "ERR invalid usage (double input array expected)\n%s", str);

   if (!mxIsComplex(argin[0])) {
      wbarray<double> A(argin[0]), Aout; 
      argout[0]=MEX_FUNCTION(nargout,argout,nargin,argin,A,Aout);
   }
   else {
      wbarray<wbcomplex> A(argin[0],0), Aout; 
      argout[0]=MEX_FUNCTION(nargout,argout,nargin,argin,A,Aout);
   }

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in wbtrace"); }
   aclu.Check();
};

template<class T>
mxArray* MEX_FUNCTION(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[],
   wbarray<T> &A,
   wbarray<T> &Aout
){
   unsigned i,j,m,n; double k;

   wbMatrix<double> D12;
   wbMatrix<unsigned> I12;
   wbMatrix<unsigned> mark;

   if (!A.isEmpty()) {
      D12.init(argin[1]); I12.initDim(D12);

      m=I12.dim1; n=A.SIZE.len;
      mark.init(m,2); 

      if (nargin>2) {
         double r;

         if (!Mx::IsDblScalar(FL,argin[2])) wblog(FL,
            "ERR invalid arg #3 (rank)");
         if (mxGetNumber(argin[2],r)) wblog(FL,"ERR %s",str);

         if (r<n) wblog(FL,"ERR invalid rank (%g/%d)",r,n);
         else if (r>n) {
            A.appendSingletons((unsigned)r);
            n=A.SIZE.len;
         }
      }

      for (i=0; i<m; i++)
      for (j=0; j<2; j++) { k=D12(i,j);
         if (k<1 || k>n) wblog(FL,
            "ERR trace index out of bounds: %g/%d (%d,%d; %d)\n"
            "Hint: use 3rd argument to specify rank.",
             k,n, i,j, mark(i,j)
         );
         else if (k!=double(unsigned(k)) || mark(i,j)++) wblog(FL,
            "ERR Invalid trace index set (%d,%d: %g/%d; %d).",
             i,j,D12(i,j),n,mark(i,j));
         else

         I12(i,j)=unsigned(D12(i,j))-1; 
      }

      A.contract(I12,Aout);
   }
   else {
      Aout.init(1,1); Aout.data[0]=0;
   }

   return Aout.toMx();
};

