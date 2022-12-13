/* =================================================================== */
/* dummy version of bracketing out MEX library if
   (i) it is not available and
   (ii) also not required in the code per se!
   in case any of the mex-routines is required
   an error message issued and the program terminates.

   NB! return types taken from Matlab R2013a documentation
   Wb,Apr02,13 */
/* =================================================================== */

typedef char mxArray;

// taken from Matlab R2013a documentation
enum mxComplexity {mxREAL=0, mxCOMPLEX};
// #define mxREAL    0
// #define mxCOMPLEX 1

#define mxDOUBLE_CLASS 10

#define MATFile FILE
#define mxChar  short
#define mwSize  size_t
#define mwIndex int

#define ERRMSG_NOMEX wblog(__FILE__,__LINE__, "ERR NOMEX (%s)", __FUNCTION__)

int wblog(const char* file, int line, const char *fmt, ...);

int mexCallMATLAB(
    unsigned nout, mxArray** aout,
    unsigned nin, mxArray** ain, const char* fstr)
{ ERRMSG_NOMEX; return 0; };

int mxIsChar   (const mxArray* a) { ERRMSG_NOMEX; return 0; };
int mxIsDouble (const mxArray* a) { ERRMSG_NOMEX; return 0; };
int mxIsComplex(const mxArray* a) { ERRMSG_NOMEX; return 0; };
int mxIsCell   (const mxArray* a) { ERRMSG_NOMEX; return 0; };
int mxIsStruct (const mxArray* a) { ERRMSG_NOMEX; return 0; };
int mxIsSparse (const mxArray* a) { ERRMSG_NOMEX; return 0; };

int mxIsEmpty  (const mxArray* a) { ERRMSG_NOMEX; return 0; };
size_t mxGetM  (const mxArray* a) { ERRMSG_NOMEX; return 0; };
size_t mxGetN  (const mxArray* a) { ERRMSG_NOMEX; return 0; };

int mxGetNzmax (const mxArray* a) { ERRMSG_NOMEX; return 0; };
int* mxGetIr   (const mxArray* a) { ERRMSG_NOMEX; return (int*)0; };
int* mxGetJc   (const mxArray* a) { ERRMSG_NOMEX; return (int*)0; };

int mxDestroyArray   (mxArray* a) { ERRMSG_NOMEX; return 0; };

void mexErrMsgTxt(const char* estr) {
   printf("\nmexErrMsgTxt:%s\nExit.\n\n", estr);
   exit(66);
};

char* mxGetClassName(const mxArray* a) {
   static char ostr[]="(none)";
   ERRMSG_NOMEX; return ostr;
};

size_t mxGetNumberOfElements(const mxArray* a) {
   ERRMSG_NOMEX; return 0;
};

int mxGetNumberOfFields(const mxArray* a) {
   ERRMSG_NOMEX; return 0;
};

unsigned mxGetNumberOfDimensions(const mxArray* a) {
   ERRMSG_NOMEX; return 0;
};

const size_t* mxGetDimensions(const mxArray* a) {
   ERRMSG_NOMEX; return (const size_t*)NULL;
};

mxArray* mxGetFieldByNumber(const mxArray *a, unsigned i, unsigned j) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

const char* mxGetFieldNameByNumber(const mxArray *a, int fid) {
   ERRMSG_NOMEX; return (const char*)NULL;
};

int mxGetFieldNumber(const mxArray *a, const char *fld) {
   ERRMSG_NOMEX; return 0;
};

int mxSetFieldByNumber(mxArray* S, unsigned i, unsigned j, mxArray* a) {
   ERRMSG_NOMEX; return 0;
};

mxArray* mxGetField(const mxArray *a, mwIndex idx, const char *fld) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

int mxAddField(const mxArray *a, const char *fld) {
   ERRMSG_NOMEX; return 0;
};

int mxGetString(const mxArray* a, char *ostr, int len) {
   ERRMSG_NOMEX; return 0;
};

mxArray* mxCreateString(const char *s) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateDoubleMatrix(size_t d1, size_t d2, char mxreal) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateNumericArray(int d, size_t *dims, int cflag, int rflag) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateSparse(size_t m, size_t n, size_t nzmax, mxComplexity cflag) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateCharArray(int d, size_t *dims) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateStructMatrix(size_t d1, size_t d2, int nf, const char** fnames) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxCreateCellMatrix(size_t d1, size_t d2) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* mxGetCell(const mxArray* S, unsigned idx) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

int mxSetCell(mxArray* S, unsigned idx, mxArray* b) {
   ERRMSG_NOMEX; return 0;
};

int mexPutVariable(const char *ws, const char *vname, mxArray *a) {
   ERRMSG_NOMEX; return 0;
};

int matPutVariable(MATFile *mfp, const char *vname, const mxArray *a) {
   ERRMSG_NOMEX; return 0;
};

mxArray* matGetVariable(MATFile *mfp, const char *vname) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

mxArray* matGetVariableInfo(MATFile *mfp, const char *vname) {
   ERRMSG_NOMEX; return (mxArray*)NULL;
};

int matDeleteVariable(MATFile *mfp, const char *name) {
   ERRMSG_NOMEX; return 0;
};

MATFile* matOpen(const char *file, const char *mode) {
   ERRMSG_NOMEX; return 0;
};

int matClose(MATFile *mfp) {
   ERRMSG_NOMEX; return 0;
};

