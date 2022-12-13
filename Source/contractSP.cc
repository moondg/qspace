char USAGE[]=
/* ================================================================== */
"   Usage: C=contractSP(A,ica,B,icb [,P]);                        \n\
                                                                  \n\
      wrapper routine to contract arbitrary-dimensional           \n\
      [wb]sparse objects                                          \n\
                                                                  \n\
   Wb,Jan24,12                                                    \n\
";
/* ================================================================== */

#ifdef MATLAB_MEX_FILE
   #define PROG mexFunctionName()
#else
   #define PROG "contractSP"
#endif

// #define CG_CHECK_MW_PERM
   #define WB_SPARSE_CLOCK

#include "wblib.h"

void mexFunction(
   int nargout, mxArray *argout[],
   int nargin, const mxArray *argin[]
){
   MX_CHECK_HELPER_NARGS(4,-1,1); 
   if (nargin>5) wblog(FL,
      "ERR invalid number of input arguments (%d)",nargin);

   wbsparray<unsigned,double> A,B,C;
   wbindex ica, icb;
   wbperm P;

   if ((!mxIsWbsparray(argin[0]) && !mxIsDouble(argin[0])) ||
       (!mxIsWbsparray(argin[2]) && !mxIsDouble(argin[2]))
     ) wblog(FL,"ERR %s() invalid numeric/sparse objects A or B",FCT);

   A.init(FL,argin[0]); ica.init(FL,argin[1],1);
   B.init(FL,argin[2]); icb.init(FL,argin[3],1);

   A.info("A");
   B.info("B");

   if (nargin>4) {
      if (!mxIsDouble(argin[4])) wblog(FL,
         "ERR %s() invalid permutation (arg #5)",FCT);
      P.wbvector<unsigned>::init(FL,argin[4]);
      P-=1; 
   }

   A.contract(FL,ica,B,icb,C,P);
   argout[0]=C.toMx(); 

};

