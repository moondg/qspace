
char USAGE[] =
/* =============================================================
 * clebsch.cc */

"   Usage:                                                    \n\
                                                              \n\
   S=getCG([opts], 'qtype',q1,q2 [,q]);                       \n\
                                                              \n\
      get CData for given symmetries. If q is not specified,  \n\
      all possible q's resulting from q1 and q2 are returned. \n\
                                                              \n\
   Options                                                    \n\
                                                              \n\
     '-h'  display this usage and exit                        \n\
     '-v'  verbose                                            \n\
                                                              \n\
   (C) Wb,Oct12,09 ; Wb,Apr26,10                              \n";

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){
   unsigned k=0; mxArray *a=NULL;
   char vflag=0;

   MX_CHECK_HELPER_NARGS(3,-1,-1); 
   if (nargin && Mx::IsEqual(argin[k],"-v")) { vflag=1; nargin--; k++; }

   if (nargin>4 || !mxIsChar(argin[k]) || 
       Mx::IsVector(argin[k+1])<=0 ||
       Mx::IsVector(argin[k+2])<=0 || (nargin>3 &&
       Mx::IsVector(argin[k+3])<=0    )
   ){
       if (nargin || k || nargout) wblog(FL,"ERR invalid usage");
       else { usage(); return; }
   }

   QVec qvec(FL,argin[k]);
   qset<gTQ> J1(FL,argin[k+1],&qvec), J2(FL,argin[k+2],&qvec), J;

   if (nargin>3) J.init(FL,argin[k+3],&qvec);

   if (J1.len<=1) {
      if (J1.len==0 || J2.len!=J1.len || qvec.len!=J1.len)
      wblog(FL,"ERR invalid usage (%d/%d/%d)",J1.len,J2.len,qvec.len);

      if (J.isEmpty()) {
         wbvector < CData<gTQ,unsigned,double>* > S;
         gCS.getCData(FL,qvec[0],J1,J2,S);
         a=S.toMxP(); 
      }
      else {
         const CData <gTQ,unsigned,double> &S=gCS.getCData(FL,qvec[0],J1,J2,J);
         a=S.toMx(); 
      }
   }
   else wblog(FL,
   "ERR multiple Q-labels not implemented yet (%d)",J1.len);

   if (vflag) gCS.printStatus(FL);

   argout[0]=a; 
}

