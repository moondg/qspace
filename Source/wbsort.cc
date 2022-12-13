#ifndef __WB_REC_HPSORT_CC__
#define __WB_REC_HPSORT_CC__

/* ------------------------------------------------------------------ *
 * Heap sort algorithm for vector array
 * Reference: Numerical Recipies C, p337f
 * AW (C) Jan 2006
 * ------------------------------------------------------------------ *
 * mergesort (ranks #2, vs. #3 for hpsort)
    + on average faster than hpsort
    + stable (does not reverse order within degenerate subspaces)
    - requires worst case an additional space of size N
    + in-place merge-sort: implemented in C++ STL (!) [stable_sort]
    + easily parallelizes!
    * while quicksort is faster in typical cases, it is
      (i) not stable, and (ii) has a worst-case performance of N^2
 * hpsort, for comparsion
    * is NOT stable (ie. may reverse order within degenerate subspaces)
    * slower on average, while in-place
    * still performs MEM_CPY even on already sorted array (!?)
    * reverse sort is actually faster than sorting a sorted array (!?)
 * Wb,Dec29,11
 * ------------------------------------------------------------------ */

// by including P, this ensures to keep initial order of degenerate
// subspaces, in that what is actually is sorted is [ra,(1:N)']

template<class T>
void Wb::hpsort(
    wbvector<T> &ra,    
    wbperm &P,
    char dir            
){
    size_t i,j,l,ir, N=ra.len; wperm_t ia,*p;
    T rra;

    char q, Q=(dir>0 ? +1 : -1);

    P.Index(N); p=P.data; if (N<2) return;

    l = (N>>1);
    ir= N-1;

    for (;;) {
        if(l>0) {                       
            rra=ra[--l];
            ia=p[l];
        } else {                        
            rra=ra[ir];                 
            ra[ir]=ra[0];               

            ia=p[ir]; p[ir]=p[0];

            if ((--ir)==0) {            
                ra[0]=rra;              
                p[0]=ia;
                break;
            }
        }

        i=l;

        j=l+l+1;
        while (j<=ir) {
            if (j<ir) {
               q=NUMCMP(ra[j],ra[j+1]); if (q==0) q=NUMCMP(p[j],p[j+1]);
               if (q!=Q) j++;       
            }

            q=NUMCMP(rra,ra[j]); if (q==0) q=NUMCMP(ia,p[j]);
            if (q!=Q) {             
                ra[i]=ra[j];
                p[i]=p[j];

                i = j;
                j = ((j+1)<<1) - 1;
            }
            else break;              
        }
        ra[i]=rra;                   
        p[i]=ia;
    }

    return;
};

template<class T>
int Wb::hpsort(
    T *A,        
    size_t lda,  
    size_t N,    
    wbperm &P,
    char dir,    
    char lex     
){
    size_t i,j,l, ia,ir, *p;

    char q, Q=(dir>0 ? +1 : -1); 

    P.Index(N); p=P.data; if (N<2) return 0;

    l=(N>>1); ir=N-1;
    T x[lda];

    for (;;) {
        if (l>0) { 
           ia=p[--l]; MEM_CPY<T>(x,lda,A+l*lda); 
        }
        else { 
            ia=p[ir];   MEM_CPY<T>(x,lda,A+ir*lda); 
            p[ir]=p[0]; MEM_CPY<T>(A+ir*lda,lda,A); 

            if ((--ir)==0) {
                p[0]=ia; MEM_CPY<T>(A,lda,x); 
                break;
            }
        }

        i=l; j=2*l+1;

        while (j<=ir) { T *aj=A+j*lda;
            if (j<ir) {
                q=Wb::recCompare(aj,aj+lda,lda,lex); 
                if (q==0) q=NUMCMP(p[j],p[j+1]); 
                if (q!=Q) { ++j; aj+=lda; } 
            }

            q=Wb::recCompare(x,aj,lda,lex); if (q==0) q=NUMCMP(ia,p[j]);
            if (q!=Q) {  
                MEM_CPY<T>(A+i*lda,lda,A+j*lda); 
                p[i]=p[j];

                i=j;
                j=((j+1)<<1)-1;
            }
            else break;              
        }
        p[i]=ia; MEM_CPY<T>(A+i*lda,lda,x);   
    }

    return 0;
};

template<class T>
int Wb::hpsort(
    wbMatrix<T> &ra, wbperm &P, char dir, char lex
 ){ return hpsort(ra.data, ra.dim2, ra.dim1, P, dir, lex); };

template<class T>
void Wb::hpsort(wbvector<T> &ra) { wbperm P; hpsort(ra,P); }

template <class T> inline
size_t Wb::matchSortedIdx_grp(
   const T *A, int lda, size_t na,
   const T *B, int ldb, size_t nb,
   wbindex &gra, wbindex &grb, size_t &mtot,
   int m,     
   char lex,  
   T eps      
){
   size_t ia,ib,iga=0,igb=0,ng=0,m1,m2,*p1,*p2;
   int nwrn=0; char c, cr=0;
   const T *a,*b;

   if (m<0) { m=lda;
      if (lda!=ldb) wblog(FL,
      "WRN %s() got different record length (%d/%d)",FCT,lda,ldb);
   }

   if (!m) wblog(FL,"ERR finite m required "
      "(m=%d; %dx%d <> %dx%d)",m,lda,na,lda,nb,ldb);
   if (m>lda || m>ldb || lda<=0 || ldb<=0) wblog(FL,
      "ERR m out of bounds (m=%d; %d/%d)",m,lda,ldb);
   if (!A || !B) wblog(FL,
      "ERR %s() got null data (0x%lX, 0x%lX)",FCT,A,B);

   gra.init(2*na); grb.init(2*nb); mtot=0;

   if (lex<=0) {
      A+=(lda-m);
      B+=(ldb-m);
   }

   a=A; b=B;
   for (ia=1; ia<na; ++ia, a+=lda) { 
      c=Wb::recCompare(a,a+lda,m,lex,eps,&nwrn);
      if (c) { cr=c; break; }
   }
   if (!cr) {
   for (ib=1; ib<nb; ++ib, b+=ldb) { 
      c=Wb::recCompare(b,b+ldb,m,lex,eps,&nwrn);
      if (c) { cr=c; break; }
   };
   if (!cr) { cr=-1; }} 

   for (a=A,b=B, ia=ib=0; ia<na && ib<nb;) {
       c=Wb::recCompare(a,b,m,lex,eps,&nwrn);

       if (c== cr) { ++ia; a+=lda; continue; } else
       if (c==-cr) { ++ib; b+=ldb; continue; } else
       if (c) wblog(FL,"ERR %s() c=%d/%d !?",FCT,c,cr);

       m1=m2=1; ++ng;
       p1=gra.data+(iga++); gra[iga++]=ia++; a+=lda;
       p2=grb.data+(igb++); grb[igb++]=ib++; b+=ldb;

       while (ia<na) {
           c=Wb::recCompare(a-lda,a,m,lex,eps,&nwrn);
           if (!c) { gra[iga++]=ia++; a+=lda; ++m1; }
           else {
              if (c!=cr) wblog(FL, 
                 "ERR %s() got unsorted recs (%d: %d/%d)",FCT,ia+1,c,cr);
              break;
           }
       }
       while (ib<nb) {
           c=Wb::recCompare(b-ldb,b,m,lex,eps,&nwrn);
           if (!c) { grb[igb++]=ib++; b+=ldb; ++m2; }
           else {
              if (c!=cr) wblog(FL, 
                 "ERR %s() got unsorted recs (%d: %d/%d)",FCT,ib+1,c,cr);
              break;
           }
       }
       mtot+=m1*m2; (*p1)=m1; (*p2)=m2;
   }

   gra.Shorten2(iga);
   grb.Shorten2(igb);

   return ng; 
};

template <class T> inline
size_t Wb::matchSortedIdx(
   const T *A, int lda, size_t na,
   const T *B, int ldb, size_t nb, wbindex &Ia, wbindex &Ib,
   int m, char lex, 
   widx_t *ma_,    
   widx_t *mb_,
   T eps            
){
   size_t i,j,k=0, ig=0,iga=0,igb=0, m1,m2,mt=0,ma=0,mb=0, ng,*p1,*p2;
   wbindex gra,grb;

   ng=Wb::matchSortedIdx_grp(A,lda,na,B,ldb,nb,gra,grb,mt,m,lex,eps);
   Ia.init(mt); Ib.init(mt); mt=0;

   for (; ig<ng; ++ig) { m1=gra[iga]; m2=grb[igb];
      if (m1==1 && m2==1) {
         Ia[k]=gra[iga+1]; iga+=2;
         Ib[k]=grb[igb+1]; igb+=2; ++k;
      }
      else {
         p1=gra.data+(iga+1); iga+=(m1+1);
         p2=grb.data+(igb+1); igb+=(m2+1);

         ma+=(m1-1);
         mb+=(m2-1); mt+=(m1*m2-1);

         for (i=0; i<m1; ++i)
         for (j=0; j<m2; ++j, ++k) { Ia[k]=p1[i]; Ib[k]=p2[j]; }
      }
   }

   if (ma_) { (*ma_)=ma; }
   if (mb_) { (*mb_)=mb; }

   return mt;
};

template <class T> inline
size_t Wb::matchSortedIdx(
   const T *A, int lda, size_t na,
   const T *B, int ldb, size_t nb,
   wbindex &Ia, wbindex &Ib, wbvector<widx_t> &D,
   int m, char lex, T eps 
){
   size_t i,j,k=0, ig=0,iga=0,igb=0, m1,m2,mt=0, ng,*p1,*p2;
   wbindex gra, grb;

   ng=Wb::matchSortedIdx_grp(A,lda,na,B,ldb,nb,gra,grb,mt,m,lex,eps);
   Ia.init(mt); Ib.init(mt); D.init(ng); mt=0;

   for (; ig<ng; ++ig) { m1=gra[iga]; m2=grb[igb];
      if (m1==1 && m2==1) {
         Ia[k]=gra[iga+1]; iga+=2;
         Ib[k]=grb[igb+1]; igb+=2; D[ig]=1; ++k;
      }
      else {
         p1=gra.data+(iga+1); iga+=(m1+1);
         p2=grb.data+(igb+1); igb+=(m2+1);

         D[ig]=m1*m2; mt+=(D[ig]-1);

         for (i=0; i<m1; ++i)
         for (j=0; j<m2; ++j, ++k) { Ia[k]=p1[i]; Ib[k]=p2[j]; }
      }
   }

   return mt;
};

#endif

