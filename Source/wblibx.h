#ifndef __WB_WBLIBX_HH__
#define __WB_WBLIBX_HH__

/* -------------------------------------------------------------------- */
/* -------------------------------------------------------------------- */

   void usage(const char *file=NULL, int line=0, const char* estr=NULL);

   int isHelpIndicator(const mxArray *a0);
   int checkHelpVersion( const mxArray *a0, mxArray **argout=NULL);

   void wbdie(const char *file, int line, const char* istr);

   void dbstop(const char* file, int line);

   char* memsize2Str(double x, char *s, unsigned l=-1);

   void version_info_toMx(mxArray *&S); 

#ifndef MEX_EXT
#if __linux__
   #define MEX_EXT mexa64
#elif __unix__
   #define MEX_EXT mexa64
#elif __APPLE__
   #define MEX_EXT mexmac*
#else
   #define MEX_EXT unknown
#endif
#endif

namespace Wb {
class VersionInfo { 
  public:
   VersionInfo() {
      fctn[0]=tag[0]=matlab[0]=mpfr[0]=compiler[0]=os[0] = 0;
      qspace[0]=git[0]=flags[0]=compiled[0] = 0;
      init();
   };

   VersionInfo& init();
   wbstring toStr() const;

   void print() const;
   mxArray* toMx() const;

   char fctn[32], tag[8], matlab[16], mpfr[16], compiler[16], os[24];
   char qspace[16], git[16], flags[128], compiled[64];

  private:

   inline unsigned str_cpy(char *s, unsigned n, const char *a);
};
}

template <class T>
class is_pointer_ { public: char q() { return 0; }; };

template <class T>
class is_pointer_<T*> { public: char q() { return 1; }; };

template <class T>
class is_pointer_<T**> { public: char q() { return 2; }; };

template <class T>
class is_pointer_<T***> { public: char q() { return 3; }; };

template <class T>
inline char ispointer(const T& x) { return is_pointer_<T>().q(); };

#endif

