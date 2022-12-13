char USAGE[] =
/* =============================================================
 * printfc.cc */

"   Usage:                                                    \n\
                                                              \n\
      MEX program that emulates matlabs input()               \n\
      with the extension that it allows for a timeout.        \n\
                                                              \n\
   AW (C) Jan 2006                                             \n";

/*    e.g. that also allows to use backspace on command line output
 *    for i=1:999, printfc([zeros(1,10)+8, sprintf('i=%03d',i)]); end
 *
 * This is a MEX-file for MATLAB.
 * ============================================================= */

#include <mex.h>

#include "wblib.h"

void mexFunction(
    int nlhs, mxArray *plhs[],
    int nrhs, const mxArray *prhs[]
){
    unsigned i,l,m,n; char *s, *s2;

    mxArray *mstr, *aa[nrhs]; memcpy(aa,prhs,nrhs*sizeof(mxArray*));

    if (nlhs>0) usage(FLINE, "No output arguments available!");

    if (nrhs) if (isHelpIndicator(prhs[0])) { usage(); return; }

    if (nrhs==1)
    if (mxIsChar(prhs[0])) {

       m=1+mxGetNumberOfElements(prhs[0])*sizeof(mxChar);

       s = new char[m];
       mxGetString(prhs[0],s,m);

       s2 = new char[m]; n=strlen(s);
       strcpy(s2,s);

       for (l=i=0; i<n; i++) {
          if (s[i]!='\\') { s2[l++]=s[i]; continue; }
          switch (s[++i]) {
             case 'n': s2[l++]='\n'; break;
             case 't': s2[l++]='\t'; break;
             case 'r': s2[l++]='\r'; break;
             case 'b': s2[l++]='\b'; break;
             case 'e': s2[l++]='\e'; break;

             default : wblog(FL,
             "unknown escape sequence %c%c<%d>", s[i-1], s[i], s[i]);
          }
       }

       s2[l]=0;
       printf(s2);

       delete [] s; delete [] s2; fflush(0);
       return;
    }

    i=mexCallMATLAB(1, &mstr, nrhs, aa, "sprintf");
    if (i!=0) return;

    m=1+mxGetNumberOfElements(mstr)*sizeof(mxChar);

    s = new char[m];

    mxGetString(mstr,s,m); 
    mexPrintf(s);

    mxDestroyArray(mstr); delete [] s; fflush(0);

    return;
}

