#ifndef __WB_BLAS_HH__
#define __WB_BLAS_HH__

// ----------------------------------------------------------------- //
// NB! matlab's mwblas/mwlapack have their own header files
// specifically, starting with matlab-R2009, the LDA type
// (e.g. in DGEMM) is no longer int but ptrdiff_t
// for compatibility reasons with 64bit linux (*@&$)!
//
// trying to avoid matlabs blas.h and lapack.h as it defines
// other troublesome stuff, eg. second as wrapper for fortran _second,
// but the latter is defined/used with maps!)*&#

// #define pINT int       // for matlab <= R2008
   #define pINT ptrdiff_t // for matlab >= R2009
// NB! all int, also info variable, replaced by pINT (!)

extern "C" {

#define dgemm dgemm_         
void dgemm (   

   const char   &,    
   const char   &,    
   const pINT   &,    
   const pINT   &,    
   const pINT   &,    
   const double &,    
         double [],   
   const pINT   &,    
         double [],   
   const pINT   &,    
   const double &,    
         double [],   
   const pINT   &     
);

#define zgemm zgemm_         
void zgemm (   

   const char   &,    
   const char   &,    
   const pINT   &,    
   const pINT   &,    
   const pINT   &,    
   const wbcomplex &, 
         wbcomplex[], 
   const pINT   &,    
         wbcomplex[], 
   const pINT   &,    
   const wbcomplex &, 
         wbcomplex[], 
   const pINT   &     
);

#define ilaenv ilaenv_

   int ilaenv(
       const pINT &,    
       const char *,    
       const char *,    
       const pINT &,    
       const pINT &,    
       const pINT &,    
       const pINT &     
   );

#define dgetri dgetri_
void dgetri(   
   const pINT   &,    
         double [],   
   const pINT   &,    
         pINT   [],   
         double [],   
   const pINT   &,    
         pINT   &     
);                    

#define dgetrf dgetrf_
void dgetrf(   
   const pINT   &,    
   const pINT   &,    
         double [],   
   const pINT   &,    
         pINT   [],   
         pINT   &     
);                    

#define dsyevd dsyevd_ 

void dsyevd(  
 const char   &,    
 const char   &,    
 const pINT   &,    
       double [],   
 const pINT   &,    
       double [],   
       double [],   
 const pINT   &,    
       pINT   [],   
 const pINT   &,    
       pINT   &     
);

#define zheev zheev_

void zheev (  
 const char   &,    
 const char   &,    
 const pINT   &,    
       wbcomplex[], 
 const pINT   &,    
       double [],   
 wbcomplex[], const pINT &,  
       double [],            
       pINT   &     
);

#define zheevd zheevd_

void zheevd(  
 const char   &,    
 const char   &,    
 const pINT   &,    
       wbcomplex[], 
 const pINT   &,    
       double [],   
 wbcomplex[], const pINT &, 
 double [],   const pINT &, 
 pINT   [],   const pINT &, 
       pINT   &     
);

#define dgeev dgeev_
void dgeev (  

 const char   &,    
 const char   &,    
 const pINT   &,    
       double [],   
 const pINT   &,    
       double [],   
       double [],   
       double [],   
 const pINT   &,    
       double [],   
 const pINT   &,    
       double [],   
 const pINT   &,    
       pINT   &     
);

#define zgeev zgeev_
void zgeev (  

 const char   &,    
 const char   &,    
 const pINT   &,    
    wbcomplex [],   
 const pINT   &,    
    wbcomplex [],   
    wbcomplex [],   
 const pINT   &,    
    wbcomplex [],   
 const pINT   &,    
    wbcomplex [],   
 const pINT   &,    
    double    [],   
       pINT   &     
);

#define dgesvd dgesvd_ 

void dgesvd (
 const char   &,    
 const char   &,    
 const pINT   &,    
 const pINT   &,    
       double [],   
 const pINT   &,    
       double [],   
       double [],   
 const pINT   &,    
       double [],   
 const pINT   &,    

       double [],   
 const pINT   &,    
       pINT   &     
);

#define dgesdd dgesdd_ 

void dgesdd (
 const char   &,    
 const pINT   &,    
 const pINT   &,    
       double [],   
 const pINT   &,    
       double [],   
       double [],   
 const pINT   &,    
       double [],   
 const pINT   &,    

       double [],   
 const pINT   &,    
       pINT   [],   
       pINT   &     
);

#define zgesvd zgesvd_ 

void zgesvd (
 const char   &,    
 const char   &,    
 const pINT   &,    
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
       double [],   
     wbcomplex[],   
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
     wbcomplex[],   
       pINT   &     
);

#define zgesdd zgesdd_ 

void zgesdd (
 const char   &,    
 const pINT   &,    
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
       double [],   
     wbcomplex[],   
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
     wbcomplex[],   
 const pINT   &,    
       double[],    
       pINT  [],    
       pINT   &     
);

}  

class STAT_DGEMM {
  public:
    STAT_DGEMM() { nn.init(32); ss.init(32); }
   ~STAT_DGEMM() {
       unsigned i; double fac=nn.max(); fac=1./(fac!=0 ? fac:1);
       if (nn.max()<1) { return; }

       wblog(FL,"%NDGEMM statistics");
       printf("\n idx_n  avg.dim        count    cost (^2.5)\n\n");
       for (i=0; i<nn.len; i++) if (nn[i])
          printf("  %3d  %8.2f %12d   %9g\n", i,
          pow(double(ss[i])/nn[i], 1/3.), nn[i],
          pow(double(ss[i])/nn[i], 2.5/3) * nn[i]*fac);
       printf("\n");
    };

    void account(unsigned n0) {
       unsigned i,n=n0;
       for (i=0; i<nn.len; i++) { n>>=1; if (n==0) break; }
       if (i>=nn.len) i=nn.len-1;
       nn[i]++; ss[i]+=n0;
    };

    wbvector<unsigned> nn;
    wbvector<unsigned long> ss;

  private:

};

#ifdef WB_CLOCK
   STAT_DGEMM stat_dgemm;
#endif

#endif

