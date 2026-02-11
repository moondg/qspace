#ifndef __WB_XMAP_HCC__
#define __WB_XMAP_HCC__

/* ------------------------------------------------------------------ *
 * class XMAP - map to store contraction pattern
 * AW (C) Jun 2006
 * ------------------------------------------------------------------ */

// ------------------------------------------------------------------ //
// ------------------------------------------------------------------ //
// central element for contraction map

class XMAP_Node {

  public:

    XMAP_Node ()
    : Z(NULL), R(NULL), I(NULL),
      tflag(0), ic(0), i1(0), i2(0), isref(0), z(1.) {};

   ~XMAP_Node () {}; 

    free(){
       if (isref) { Z=R=NULL; }
       else {
          if (Z) { delete Z; Z=NULL; }
          if (R) { delete R; R=NULL; }
       }
    };

    void init(const mxArray *a);

    wbvector<unsigned> size() const {
       if (Z) return (*Z).SIZE; else
       if (R) return (*R).SIZE;
       else   return wbvector<unsigned>();
    };

    unsigned rank() const {
       if (Z) return (*Z).rank(); else
       if (R) return (*R).rank();
       else   return 0;
    };

    bool sameSize(const XMAP_Node &B) const {
       if ((!R ^ !Z) || (!B.R ^ !B.Z)) wblog(FL,
       "ERR Severe XMAP inconsistency.");

       return (
          (R   ? (  *R).SIZE : (  *Z).SIZE) !=
          (B.R ? (*B.R).SIZE : (*B.Z).SIZE)
       );
    };

    template<class T>
    void set(
       const wbArray<T> &A,
       unsigned ic, unsigned i1, unsigned i2, wbcomplex z, char ref=1
    );

    void set( 
       unsigned ic, unsigned i1, unsigned i2, wbcomplex z, char ref=1
    );

    wbArray<wbcomplex> *Z; 
    wbArray<double>    *R;

    char tflag;
    wbperm P;

    unsigned ic;    
    unsigned i1;    
    unsigned i2;    

    unsigned isref;

    wbcomplex z;    

  protected:
  private:

};

class XMAP : public wbvector<XMAP_Node> {
  public:
     XMAP(unsigned n) : isref(1) { init(n); }; 

    ~XMAP()  {}; 

     free()  {
        for (unsigned i=0; i<len; i++)
        data[i].free();
     };

     void init(unsigned n) { free(); init(n); };

     template<class T>
     void set(
        unsigned k, const wbArray<T> &A,
        unsigned ic, unsigned i1, unsigned i2, wbcomplex z
     ) { data[k].set(A,ic,i1,i2,z, isref); };

     void set(
        unsigned k,
        unsigned ic, unsigned i1, unsigned i2, wbcomplex z
     ) { data[k].set(ic,i1,i2,z, isref); };

     char isref;

  protected:
  private:
};

void XMAP_Node::init(const mxArray *a) {

   if (!a) wblog(FL,
   "ERR Input mxArray is NULL.");

   if (!mxIsDouble(a) || mxIsSparse(a)) wblog(FL,
   "ERR Need full double array upon input.");

   free();
   (*( mxIsComplex(a) ? Z : R )).init(a);
}

template<class T>
void XMAP_Node::set(
   unsigned k, const wbArray<T> &A,
   unsigned ic_, unsigned i1_, unsigned i2_, wbcomplex z_, char ref
){
   wblog(FL,"ERR XMAP::set() not available for data type %s",
   TSTR(T));
};

void XMAP_Node::set(
   unsigned k, const wbArray<wbcomplex> &A,
   unsigned ic_, unsigned i1_, unsigned i2_, wbcomplex z_, char ref
){
   ic=ic_; i1=i1_; i2=i2_; z=z_;

   if (A.rank()!=2) wblog(FL,
   "ERR expecting rank-2 tensors to contract.");

   if (ref) {
      Z=&A; isref=1;
   }
   else {
      if (!Z) Z = new wbArray<wbcomplex>;
      (*Z)=A; (*Z)*=z;
      z=1.;
   }
};

void XMAP_Node::set(
   const wbArray<double> &A,
   unsigned ic_, unsigned i1_, unsigned i2_, wbcomplex z_, char ref
){
   ic=ic_; i1=i1_; i2=i2_; z=z_;

   if (A.rank()!=2) wblog(FL,
   "ERR expecting rank-2 tensors to contract.");

   if (ref) {
      R=&A; isref=1;
   }
   else {
      if (z.i==0.) {
         if (!R) R = new wbArray<double>;
         (*R)=A; (*R)*=z.r;
      }
      else {
         if (!Z) Z = new wbArray<wbcomplex>;
         (*Z).set(A, wbArray<double>());
         (*Z)*=z;
      }
      z=1.;
   }
};

void XMAP_Node::set(
   unsigned ic_, unsigned i1_, unsigned i2_, wbcomplex z_
){
   ic=ic_; i1=i1_; i2=i2_; z=z_;

   if (Z && R) wblog(FL,
   "ERR Both Z and R are set (0x%lX, 0x%lX) ???", Z, R);

   if (Z) (*Z)*=z; else
   if (R) {
      if (z.i==0) (*R)*=z.r;
      else {
         Z = new wbArray<wbcomplex>;
         Z.set(R, wbArray<double>()); (*Z)*=z;
         delete R; R=NULL;
      }
   };

   z=1.;
};

void XMAP::init(const mxArray *a) {
   unsigned i,N;

   if (mxIsCell(a)) {
      N=mxGetNumberOfElements(a);
      init(N);

      for (i=0; i<N; i++)
      data[i].init(mxGetCell(a,i));
   }
   else
   if (mxIsDouble(a)) {
      N=1; init(N); data[0].init(a);
   }

   for (i=0; i<N; i++) {
      if (data[i].rank()!=2) wblog(FL,
      "ERR Expecting rank-2 tensors (%d)", data[i].rank());
   }
}

#endif

