char USAGE[]=
/* ==================================================================
 * var2mat.cc
 *
 *  C++ routine to set/replace variable in mat file. */

"  Usage: var2mat [-q] <matfile> <varname> <value> \n\
                                                   \n\
     Flags: -q  quite, i.e. no output unless error.\n\
                                                   \n\
  Wb,Aug21,06                                      \n\
";

/* This is a MEX-file for MATLAB.
 * ================================================================== */

#define MAIN

#include "wblib.h"

int main (int nargin, char **argin) {

  int i, q=0, ex=0, m=argin[1][0]=='-' ? 1 : 0;

  { if (nargin>1) if (checkHelpVersion(argin[1])) { return 0; }

    if (nargin!=4+m) {
       fprintf(stderr,
       "ERR Invalid number of arguments (%d/%d)\n", nargin-1-m, 3);
       exit(1);
    }

    if (m) { q=1;
       if (strcmp(argin[1],"-q")) {
          fprintf(stderr,"ERR Invalid flag `%s'\n",argin[1]);
          exit(1);
       }
  } }

  double x;
  mxArray *a;
  const char *fname=argin[1+m], *vname=argin[2+m];

  Wb::matFile F;

  if (charGetNumber(FL, argin[3+m], x)) wblog(FL,"ERR %s",str);

  F.Open(FL,fname,"u"); 

  a=matGetVariableInfo(F.mfp,vname);
  if (a) { ex=1; mxDestroyArray(a); }

  i=matPutVariable(F.mfp, vname, numtoMx(x));

  if (i) {
     fprintf(stderr, "\n  ERR Failed to %s %s=%g %s file `%s' (%d)\n\n",
     ex ? "replace" : "write", vname, x, ex ? "in" : "to", fname, i);
  }
  else if (!q) {
     fprintf(stdout, "\n Successfully %s value %s=%g %s file `%s'\n",
     ex ? "replaced" : "wrote", vname, x, ex ? "in" : "to", fname);
  }

  return 0;
}

