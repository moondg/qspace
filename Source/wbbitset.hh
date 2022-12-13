#ifndef __WB_BITSET_HH__
#define __WB_BITSET_HH__

// ---------------------------------------------------------------- //
   namespace Wb {
// ---------------------------------------------------------------- //

// arbitrary length bitset
   class Bitset;

   template <class T> class bitset;

template <class T>
void SetBit(T &x, unsigned i) { 
   unsigned nbits=8*sizeof(T);
   if (i>=nbits) wblog(FL,
      "ERR %s() bit index out of range (%d/%d)",FCT,i,nbits);
   T r=1; r<<=i; x|=r;
};

template <class T> inline
T setbit(const T &x0, unsigned i) { T x=x0; SetBit(x,i); return x; };

template <class T>
void UnsetBit(T &x, unsigned i) { 
   unsigned nbits=8*sizeof(T);
   if (i>=nbits) wblog(FL,
      "ERR %s() bit-index out of range (%d/%d)",FCT,i,nbits);
   T r=1; r=(~(r<<i)); x&=r;
};

template <class T> inline
T unsetbit(const T &x0, unsigned i) { T x=x0; UnsetBit(x,i); return x; };

template <class T>
bool issetbit(const T &x, unsigned i) { 
   unsigned nbits=8*sizeof(T);
   if (i>=nbits) wblog(FL,
      "ERR %s() bit-index out of range (%d/%d)",FCT,i,nbits);
   T r=1; r<<=i;
   return (x&r ? 1 : 0);
};

template <class T>
wbstring bitstr(const T &x); 

class Bitset: public wbvector<char> {  

  public:

    Bitset(unsigned n=0   ) : nbits(0) { NEW_BITSET(n); };
    Bitset(const Bitset &b)
    : wbvector<char>(), nbits(0) { NEW_BITSET(b.nbits, b.data); };

    Bitset& init(unsigned n=0) { return NEW_BITSET(n); };
    Bitset& init(unsigned n, const char* b) { return NEW_BITSET(n,b); };

    Bitset& setall(bool t) {
       memset(data, t ? 255:0, len);
       if (t && len) { (data[len-1])>>=(8*len-nbits); }
       return *this;
    };

    Bitset& set(unsigned i) { 
       if (i>=nbits) wblog(FL,
          "ERR %s() bit-index out of range (%d/%d)",FCT,i,nbits);
       unsigned k=i/8; data[k] |= (1<<(i-8*k));
       return *this;
    };

    Bitset& unset(unsigned i) { 
       if (i>=nbits) wblog(FL,
          "ERR %s() bit-index out of range (%d/%d)",FCT,i,nbits);
       unsigned k=i/8; data[k] &= ~(1<<(i-8*k));
       return *this;
    };

    bool isset(unsigned i) { 
       if (i>=nbits) wblog(FL,
          "ERR %s() bit-index out of range (%d/%d)",FCT,i,nbits);
       unsigned k=i/8;
       return (data[k] & (1<<(i-8*k)) ? 1 : 0);
    };

    unsigned nnz() const {
       unsigned i,n=0;
       for (i=0; i<len; ++i) {
          const char &b=data[i];
          for (char r=1; r; r<<=1) { if (b&r) ++n; }
       }
       return n;
    };

    bool operator==(const Bitset &b) const {
       if (len!=b.len || nbits!=b.nbits) return 0;
       if (data!=b.data) {
          for (unsigned i=0; i<len; i++)
          if (data[i]!=b.data[i]) return 0;
       }
       return 1;
    };

    bool operator!=(const Bitset &b) const { return !((*this)==b); };

    wbstring toStr() const;

    mxArray* toMx() const;
    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    void put(const char *vname, const char *ws="base"
    ) const { put(0,0,vname,ws); };

    void put(const char *file, int line,
      const char *vname, const char *ws="base"
    ) const;

    unsigned nbits;

  protected: 

    Bitset& NEW_BITSET(unsigned n, const char* d=NULL) {
       nbits=n; if (n) RENEW(((n-1)/8)+1,d); 
       return *this;
    }

}; 

template <class T>
class bitset { 

  public:

    bitset() : data(0) {};
    bitset(const bitset<T> &b) : data(0) { data=b.data; };

    bitset& setall(bool t) {
       memset(&data, t ? 255:0, sizeof(T));
       return *this;
    };

    bitset& set(unsigned i) { 
       T r=1; r<<=i;
       if (r) { data |= r; } else wblog(FL,
         "WRN %s() bit-index out of range (%d/%d)",FCT,i,8*sizeof(T));
       return *this;
    };

    bitset& unset(unsigned i) { 
       T r=1; r<<=i;
       if (r) { data &=(~r); } else wblog(FL,
         "WRN %s() bit-index out of range (%d/%d)",FCT,i,8*sizeof(T));
       return *this;
    };

    bool isset(unsigned i) { 
       T r=1; r<<=i;
       if (!r) wblog(FL,
         "WRN %s() bit-index out of range (%d/%d)",FCT,i,8*sizeof(T));
       return (data&r ? 1 : 0);
    };

    unsigned nnz() const {
       unsigned n=0;
       for (T r=1; r; r<<=1) { if (data&r) ++n; }
       return n;
    };

    bool operator==(const bitset &b) const { return (data==b.data); };
    bool operator!=(const bitset &b) const { return (data!=b.data); };

    wbstring toStr() const;

    mxArray* toMx() const;
    mxArray* mxCreateStruct(unsigned m, unsigned n) const;
    void add2MxStruct(mxArray *S, unsigned i) const;

    void put(const char *vname, const char *ws="base"
    ) const { put(0,0,vname,ws); };

    void put(const char *file, int line,
      const char *vname, const char *ws="base"
    ) const;

    T data;

  protected: 

}; 

}; 

#endif

