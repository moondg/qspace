#ifndef __WB_REC_HPSORT_HH__
#define __WB_REC_HPSORT_HH__

namespace Wb {

/* ------------------------------------------------------------------ */
// generalized matrix/record sort
// sorted by rows/records of length m, with actual leading length lda
/* ------------------------------------------------------------------ */

template<class T>
int hpsort(
    T *A,          
    size_t lda,    
    size_t N,      
    wbperm &P,
    char dir=+1,   
    char lex=+1    
);

template<class T>
int hpsort(
    wbMatrix<T> &ra,
    wbperm &P,
    char dir=+1,   
    char lex=+1    
);

template<class T>
void hpsort(
    wbvector<T> &ra,  
    wbperm &P,
    char dir=+1    
);

template<class T>
void hpsort(wbvector<T> &ra);

template <class T> inline
size_t matchSortedIdx_grp(
    const T *da, int lda, size_t na,
    const T *db, int ldb, size_t nb,
    wbindex &gra, wbindex &grb, size_t &mtot,
    int m,    
    char lex, 
    T eps=0   
);

template <class T> inline
size_t matchSortedIdx(
    const T *da, int lda, size_t na,
    const T *db, int ldb, size_t nb, wbindex &Ia, wbindex &Ib,
    int m, char lex,  
    widx_t *ma=NULL, 
    widx_t *mb=NULL,
    T eps=0           
);

template <class T> inline
size_t matchSortedIdx(
    const T *da, int lda, size_t na,
    const T *db, int ldb, size_t nb,
    wbindex &Ia, wbindex &Ib, wbvector<widx_t> &D, 
    int m, char lex, T eps=0  
);

}; 

#endif

