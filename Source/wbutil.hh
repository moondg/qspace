/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace utility routines
 *
 * Copyright 2024 Andreas Weichselbaum
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * --------------------------------------------------------------------- */

#ifndef __WBUTIL_HCC__
#define __WBUTIL_HCC__

// NB! util.hh cannot use wbstring etc. as members of class
// definitions here since wbstring.hh etc. need to be fully
// defined (if an object is used as a member of a class,
// the compiler wants to know its full definition, so the
// forward declarations in wblib.h do not suffice)
// Wb,Sep03,16

namespace Wb {

class Path { 
  public:

    Path() {};
    Path(const char *s) { init(FL,s); };
    Path(const char *F, int L, const char *s) { init(F,L,s); };

    int init() { buf.init(); idx.init(); return 0; };

    int init( 
       const char *F, int L, const char *s, const char *delim=":",
       unsigned n_=-1
    );

    int append( 
       const char *F, int L, const char *s, const char *delim=":"
    );

    char* operator[] (unsigned i) {
       if (i+1>=idx.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,idx.len-1);
       return buf.data+idx[i];
    };
    const char* operator[] (unsigned i) const { 
       if (i+1>=idx.len) wblog(FL,
          "ERR %s() index out of bounds (%d/%d) !?",FCT,i,idx.len-1);
       return buf.data+idx[i];
    };

    const char* last() const {
       if (idx.len<2) {
          if (!idx.len)
               wblog(FL,"ERR %s() path.got uninitialized path !?",FCT);
          else wblog(FL,"ERR %s() path.got path with idx.len=%d !?",FCT,idx.len);
       }
       return (*this)[idx.len-2];
    };

    unsigned numel() const {
       if (idx.len<2) {
          if (!idx.len) wblog(FL,
             "ERR %s() got uninitialized Path (%d,%d)",FCT,idx.len,buf.len);
          else wblog(FL,"ERR %s() got Path with idx.len=%d",FCT,idx.len);
       }
       return (idx.len-1);
    };

    Wb::iterator<Path,char*> begin() {
       return Wb::iterator<Path,char*>(this); 
    };
    Wb::citerator<Path,char*> cbegin() const { 
       return Wb::citerator<Path,char*>(this);
    };

    void print(const char *istr="") const;

    wbstring buf;
    wbvector<unsigned> idx;

  protected:
  private:

};

void Path::print(const char *istr) const {
   unsigned i=0, n=numel();

   if (istr && istr[0]) wblog(FL,"==> Contents of %s%N",istr);
   for (; i<n; ++i) { PRINTF("  %4d. %s\n",i+1,(*this)[i]); }
   PRINTF("\n");
};

int Path::init(
   const char *F, int L, const char *s, const char *delim,
   unsigned n_
){

   if (!s || !s[0]) {
      if (int(n_)>0) buf.init(n_); else buf.init();
      idx.init(); return 0;
   }

   if (!delim || !delim[0]) wblog(F_L,
      "ERR %s() invalid delim '%s' !?",FCT,delim?delim:"(null)");

   unsigned i,m0,m=0, n=strlen(s); char *b,*b0;

   if (int(n_)<0) { buf=s; }
   else {
      if (n_<n) wblog(F_L,"ERR %s() "
         "input path length out of bounds %d/%d !?",FCT,n,n_);
      buf.init(n_); strcpy(buf.data,s);
   }

   b0=b=buf.data;

   while ((b0=strsep(&b,delim))) { ++m; }

   idx.init(m+1); m0=m; b=buf.data;

   for (m=i=0; i<=n; ++i) {
      idx[m++]=i;

      if (!b[i]) continue; 

      do { ++i; } while (b[i] && i<=n); 
   }

   if (m!=m0) wblog(FL,"ERR %s() m=%d/%d !?",FCT,m,m0);
   idx[m]=n;

   if (idx.last()>=buf.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d) !?",FCT,idx.last(),buf.len);
   return m;
};

int Path::append(
   const char *F, int L, const char *s, const char *delim
){
   if (!idx.len) { return init(F,L,s,delim); }

   if (!s || !s[0]) { buf.init(); idx.init(); return 0; }
   if (!delim || !delim[0]) wblog(FL,
      "ERR %s() invalid delim '%s' !?",FCT,delim?delim:"(null)");

   unsigned i,l=idx.len-1,m0,m=0, n=strlen(s);
   unsigned i0=idx[l]+1; 
   char *b=buf.data+i0, *b0=b;

   if (i0+n>=buf.len) wblog(FL,
      "ERR %s() index out of bounds (%d+%d/%d) !?",FCT,i0,n,buf.len);
   strcpy(b0,s);

   while ((b0=strsep(&b,delim))) { ++m; }

   idx.Resize(idx.len+m);
   m0=m; b=buf.data+i0;

   for (m=i=0; i<=n; ++i) {
      idx[l+(m++)]=i0+i;

      if (!b[i]) continue; 

      do { ++i; } while (b[i] && i<=n); 
   }

   if (m!=m0) wblog(FL,"ERR %s() m=%d/%d !?",FCT,m,m0);
   idx[l+m]=n;

   if (idx.last()>=buf.len) wblog(FL,
      "ERR %s() index out of bounds (%d/%d) !?",FCT,idx.last(),buf.len);

   return m;
};

double getTimeNow() {

   struct timespec now;
   int e=clock_gettime(CLOCK_REALTIME,&now);
   if (e<0) wblog(FL,"ERR %s() gettime returned error (e=%d)",FCT,e);

   return (now.tv_sec + (now.tv_nsec/1e9));
};

template<class T>
void getPerms(wbMatrix<T> &PP,unsigned r, bool l2r=0) {
   if (r<=2) {
      if (!r) { PP.init(); } else
      if (r==1) { PP.init(1,1); }
      else {
         const T pp[4] = {0, 1, 1, 0};
         PP.init(2,2,pp);
      }; return;
   }
   else if (r>12) {
      double N=1, i=0; for (; i<r; ++i) { N*=(i+1); }
      wblog(FL,"ERR %s() n=%d exceeds bounds (=> %.4g perms)",FCT,r,N);
   }

   unsigned i,j,l,n,p,k,k_,b=0; size_t N=1;

   for (i=0; i<r; ++i) {  N*=size_t(i+1); }
   PP.init(N,r);

   if (l2r) {
      for (i=0; i<r; ++i) { PP(0,i)=PP(1,i)=i; };
      SWAP(PP(1,0),PP(1,1)); l=2; 

      for (p=2; p<r; ++p) { n=l;  
         for (j=0; j<n; ++j, ++b) {
         for (k=0; k<p; ++k, ++l) {
            k_=( b%2 ? k: p-k-1 );

            if (l>=N) wblog(FL,"ERR %s() index out of bounds %d/%d",FCT,l,N);
            for (i=0; i<k_; ++i) { PP(l,i)=PP(j,i); }; PP(l,i)=p;
            for (++i; i<=p; ++i) { PP(l,i)=PP(j,i-1); }
            for (   ; i<r;  ++i) { PP(l,i)=i; }
         }}
      }
   }
   else {
      for (i=0; i<r; ++i) { PP(0,i)=PP(1,i)=i; };
      SWAP(PP(1,i-2),PP(1,i-1));
      l=2; 

      for (p=r-3; p<r; --p) { n=l;  
         for (j=0;   j<n; ++j, ++b) {
         for (k=p+1; k<r; ++k, ++l) {
            k_=( b%2 ? p+r-k : k);

            if (l>=N) wblog(FL,"ERR %s() index out of bounds %d/%d",FCT,l,N);
            for (i=0; i<p;  ++i) { PP(l,i)=i; }
            for (   ; i<k_; ++i) { PP(l,i)=PP(j,i+1); }; PP(l,i)=p;
            for (++i; i<r;  ++i) { PP(l,i)=PP(j,i); }
         }}
      }
   }
   if (l!=N) wblog(FL,"WRN %s() index mismatch %d/%d",FCT,l,N);
};

class counter { 
  public:
    counter() : n(0) {};

    counter(const char *F, int L, const char *fmt_) : n(0) {
       file=F; line=L; fmt=fmt_;
    };

   ~counter() { if (n) {
       if (!fmt) wblog(FL,"ERR %s() got undefined fmt",FCT);
       unsigned n=128; char s[n];
       snprintf(s,n,fmt.data,n);
       if (file)
            { wblog(file.data, line,"--> %s",s); }
       else { wblog(FL,             "--> %s",s); }
    }};

    counter& operator++() { ++n; return *this; }

    long n;
    int line;
    wbstring file, fmt;

  protected:
  private:
};

}; 

#endif
