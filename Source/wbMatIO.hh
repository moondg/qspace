#ifndef __WB_MatIO_HCC__
#define __WB_MatIO_HCC__

// IO handling of mat files (MatLab binaries)

namespace Wb {

class matFile { 
   public:

      matFile() : mfp(nullptr) {};

      matFile(const char *F, int L, const char *fname, const char *rw
      ) : mfp(nullptr) { open(F,L,fname,rw,1); }; 

      matFile(const matFile &X) : mfp(nullptr) {
         wblog(FL,"ERR do not copy matFile object! (0x%lX)",&X);
      };

     ~matFile() {
         init();
      };

      void close() { init(); };
      void init() {
         if (mfp) if (matClose(mfp)) {
            wblog(F.data,L,"WRN Could not close file %s%s /%s !?",
            strlen(fname.data)>20 ? "\n" : "", fname.data, rw.data);
         }
         mfp=nullptr; F.init(); L=0; fname.init(); rw.init();
      };

      int Open(const char *F, int L, const char *FName, const char *RW
      ){ return open(F,L,FName,RW,1); }; 

      int open(
         const char *file, int line, const char *FName, const char *RW,
         char xflag=0 
      ){
         if (mfp) init();
         F=file; L=line; fname=FName;
         if (!RW || !RW[0]) wblog(FL,
            "ERR %s() got empty read/write flag",FCT);

         if (!strcmp(RW,"w")) {
            rw="w7.3";
         }
         else if (!strcmp(RW,"w:std")) {
            rw="wz";
         }
         else { rw=RW; }

         mfp=matOpen(fname.data, rw.data);
         if (!mfp) {
            if (xflag) { unsigned n=strlen(fname.data);
               if (n<50)
                  wblog(F.data,L,"ERR failed to %s access%s%s"
                  "\n(hint: non-existing directory?)",
                  rw.data, n<30?" ":" file\n",fname.data);
               else {
                  char *x,cwd[256]; cwd[0]=0; x=getcwd(cwd,256);
                  wblog(F.data,L,"ERR failed to access (%s) ... %N%N"
                  "   mat: %s !?%N"
                  "   cwd: %s%N"
                  " (hint: non-existing directory or no rw permission)%N",
                  rw.data, fname.data, x ? cwd : "(pwd)");
               }
            }
            init(); return 0;
         }
         return 1;
      };

      matFile& put(
         const char *file, int line,
         const char *vname, const mxArray* a
       ){
         int i=matPutVariable(mfp,vname,a);
         if (i) {
            wblog(file ? file : F.data, line ? line : L,
              "ERR matPutVariable() returned %d \n(%s,%s; 0x%lX) !?\n"
              "hint: possibly out of disc space", i,
               fname.data ? fname.data : "(undefined)",
               vname ? vname : "(empty)", mfp
            );
         }
         return *this;
      };
      matFile& put(const char *vname, const mxArray* a){
         return put(0,0,vname,a);
      };

      int existsVar(const char* vn) const {
         mxArray *a=matGetVariableInfo(mfp,vn);
         int r=(a!=nullptr);
         mxDestroyArray(a);
         return r;
      };

      int deleteVar(const char* vn) const {
         int i=matDeleteVariable(mfp,vn);
         return i;
      };

      int deleteVar(const char *f, int L, const char* vn) const {
         int i=matDeleteVariable(mfp,vn);
         if (i) wblog(f,L,
            "ERR Failed to delete variable `%s'\n"
            "from file `%s' (0x%lX; %d/%d ?)",
         vn, fname.data, mfp, existsVar(vn), i);
         return i;
      };

      MATFile *mfp;

      wbstring F, fname, rw; 
      int L;

   protected:
   private:
};

} 

#endif

