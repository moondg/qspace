#ifndef __WB_WBLIB_HH__
#define __WB_WBLIB_HH__

// forward declarations

template <class T>  class wbvector;
template <class T>  class wbvec;
template <class T>  class wbMatrix;
template <class T>  class wbarray;
template <class TD> class wbsparray;
template <class T>  class sparseIndex2D;

template <class TQ, class TD> class QSpace;

class bitset;
class wbcomplex;
class wbindex;
class wbperm;
class wbstring;

class iTags;
class ctrIdx;
class QType;
class QVec;
class QDir;

namespace CG {
   class FileLock;
}

namespace Wb {
   class Clock;
   class ClockSet;
}

using namespace std;

#ifdef MATLAB_MEX_FILE
class MXPut;
#endif

#define ENABLE_IF_isPOD(T__) typename \
   std::enable_if<  WbUtil<T__>::isPOD(), T__ >::type* = nullptr

#define ENABLE_IF_noPOD(T__) typename \
   std::enable_if< !WbUtil<T__>::isPOD(), T__ >::type* = nullptr

#define ENABLE_IF_isINT(T__) typename \
   std::enable_if<std::is_integral<T__>::value>* = nullptr

#define ENABLE_IF_isInt(T__) typename \
   std::enable_if< WbUtil<T__>::isInt(), T__ >::type* = nullptr

#define ENABLE_IF_isFloat(T__) typename \
   std::enable_if< (WbUtil<T__>::isFloat()&3)!=0, T__ >::type* = nullptr

#define ENABLE_IF_isFloat0(T__) typename\
   std::enable_if< (WbUtil<T__>::isFloat()&3)==0, T__ >::type* = nullptr

#define ENABLE_IF_not_isFloat(T__) typename \
   std::enable_if< (WbUtil<T__>::isFloat()&3)==0, T__ >::type* = nullptr

#define ENABLE_IF_isFloat1(T__) typename\
   std::enable_if< (WbUtil<T__>::isFloat()&3)==1, T__ >::type* = nullptr

#define ENABLE_IF_isFloat2(T__) typename\
   std::enable_if< (WbUtil<T__>::isFloat()&3)==3, T__ >::type* = nullptr

#define ENABLE_IF_isComplex(T__) typename \
   std::enable_if<  WbUtil<T__>::isComplex(), T__ >::type* = nullptr

namespace Wb {
   class myname__ {
     public:
       myname__() {
          #ifdef PROG
             sprintf_str("%s",PROG);
          #elif defined MATLAB_MEX_FILE
             sprintf_str("%s",mexFunctionName());
          #else
             sprintf_str("(?myname?)");
          #endif

           unsigned i=0, k=0, n=24;
           for (; str[i]; ++i) { if (str[i]=='/') k=i+1; }

           memset(data,0,32);
           strcpy(data+20,"tag=`"); 
           strncpy(data,str+k,n);                 

          #ifdef PROG_TAG
             snprintf(data+26,4,"%.3s",PROG_TAG); 
          #else
             memcpy(data+26,data,3); data[29]=0;
          #endif
       };

       char data[32];

     private:
   };

   myname__ myName; 

#define myname Wb::myName.data
#define mytag  Wb::myName.data+26 

   WBL_COLOR_SCHEME useCol = WLC_DARK; 

   int get_WB_VERBOSE(const char *F=0, int L=0);
   int got_DBSTOP(const char *F=0, int L=0);
   int got_DESKTOP();
   int is_DEPLOYED(const char *F=0, int L=0);

   int envVRB= 0; 
   int envDKT=-1; 

   int envDBG= 0; 

   int envFERM=0;   

   int envFullOM=0; 

   int my_caller_tid=0; 

}; 

void dbstop(const char* F, int L);
void ExitMsg(const char* s="", char xflag=0);
void doflush();

namespace Wb {
   const char* class_function(const char *s, char vflag=0);

   int system_tid();
};

   void usage(const char *F, int L, const char* fmt=NULL, ...);
   void usage(const char *F=NULL, int L=0) { usage(F,L,NULL); };

   int isHelpIndicator(const char *s);

#ifdef MATLAB_MEX_FILE
   int isHelpIndicator(const mxArray *a0);
   int checkHelpVersion(const mxArray *a0, mxArray **argout=NULL);
#else
   int checkHelpVersion(const char *s);
#endif

   void wbdie(const char *F, int L, const char* istr);

   void dbstop(const char* F, int L);

namespace Wb {
   char* memsize2Str(double x, char *s, unsigned l=-1);
}

namespace wbl { 
   const int SLEN=256;  

   int status=0;   
};

namespace Wb {

   class gpara__ {
     public:
       gpara__();
       void init(char force=0);
       void info(const char *F=NULL, int L=0) const; 

    private:
       time_t tlast;
   };

   gpara__ gpara;

class CleanUp { 
  public:
    CleanUp() { 
       gpara.init(); 
    };

   void Check();
   ~CleanUp();
 private:
};
}; 

   void version_info_toMx(mxArray *&S); 

#ifndef MEX_EXT
#  if __linux__
#     define MEX_EXT mexa64
#  elif __unix__
#     define MEX_EXT mexa64
#  elif __APPLE__
#     define MEX_EXT mexmac*
#  else
#     define MEX_EXT unknown
#  endif
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
inline char ispointer(const T& x QS_UNUSED_VAR) { return is_pointer_<T>().q(); };

namespace Wb {
   class ARGV { 
     public:
       ARGV(va_list *vl=NULL) : args(vl) {};
      ~ARGV() { if (args) { va_end(*args); args=NULL; }};

     private:
       va_list *args;
   };

template<class T>
class tmpSet__ { 
  public:
    tmpSet__(T &v_) : v(&v_), x_done(v_) {};

    tmpSet__(T &v_, T x) : v(&v_), x_done(v_) { *v=x; };
    tmpSet__(T &v_, T x, T xd) : v(&v_), x_done(xd) { *v=x; };

   ~tmpSet__();

    tmpSet__& set(T x, char op='=');
    tmpSet__& set(T x, char op, T xd) {
       set(x,op); x_done=xd; return *this; };

    T *v, x_done; 
};

template<class T>
class iterLevel { 
  public:
    iterLevel(T *v_, T dx_=1);

    iterLevel(const char *F, int L, const char *fct,
       T *v_, const char *istr=NULL, T dx_=1);

   ~iterLevel();  

  private:
    T *v, dx;   

    char *file; 
    int line;

}; 

}; 

#endif

