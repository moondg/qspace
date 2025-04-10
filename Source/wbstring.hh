/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : wbstring
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

#ifndef __WB_STRING_HH__
#define __WB_STRING_HH__

namespace Wb {

   char* memsize2Str(double x, char *s, unsigned l);

   wbstring repHome(const wbstring &file);

   int findRegEx(const char *s, const char *pat, char icase=0);

   size_t string_hash(const char *s, unsigned n=-1, size_t hid=0);
};

class wbstring : public wbvector<char> { 

  public:

    wbstring (unsigned n=0) : wbvector<char>(n) {};

    wbstring (unsigned n, const char c)
     : wbvector<char>(n+1) { set(c); data[len-1]=0; };

    wbstring (const char* s,unsigned n) { init(s,n); }

    wbstring (const wbstring &s) : wbvector<char>() { init(s.data); };
    wbstring (const char* s1) : wbvector<char>() { init(s1); };

    wbstring (
       const char* s1, const char* s2,
       const char* s3=NULL, const char* s4=NULL
    ) : wbvector<char>() {
       unsigned
          l1 = (s1? strlen(s1) : 0),
          l2 = (s2? strlen(s2) : 0),
          l3 = (s3? strlen(s3) : 0),
          l4 = (s4? strlen(s4) : 0);

       RENEW(l1+l2+l3+l4+1);
       if (l1) { strcpy(data,         s1); }
       if (l2) { strcpy(data+l1,      s2); }
       if (l3) { strcpy(data+l1+l2,   s3); }
       if (l4) { strcpy(data+l1+l2+l3,s4); }
       data[len-1]=0; 
    };

#ifndef NOMEX
    wbstring(const mxArray *a) : wbvector<char>() { init(0,0,a); };
    wbstring& init(const mxArray *a) { return init(0,0,a); }
#endif

    wbstring (const char *F, int L, const mxArray *a)
     : wbvector<char>() { init(F,L,a); };

#if 0 
    ~wbstring() {
       if (len) { 
          if (data[len-1]) {
             unsigned i=0; for (; i<len && data[i]; ++i) { };
             if (i>=len) { wblog(FL,
                "WRN %s() string out of bounds? (len=%d/%d)",FCT,i,len);
                fprintf(stdout,"\n");
                for (i=0; i<len && data[i]; ++i) fputc(data[i],stdout);
                fprintf(stdout,"\n\n");
             }
          }
       }
    };
#endif

    wbstring& init(unsigned  n=0) {  RENEW(n); return *this; };

    wbstring& init(const char *s) {
       if (s)  {
          if (s[0])
               { RENEW(strlen(s)+1, s); data[len-1]=0; }
          else { RENEW(1); data[0]=0; } 
       } else  { RENEW(0); } 
       return *this;
    };

    wbstring& init(const char* s,unsigned n) {
       unsigned l=(s? strlen(s) : -1);
       RENEW(n+1,s,l); data[n]=0; 
       return *this;
    };

    wbstring& init(const char *F, int L, const mxArray *a) { 
       int i,n;
       if (!a || !(n=mxGetNumberOfElements(a))) {
          RENEW(0); return *this;
       };

       if (!mxIsChar(a)) wblog(F_L,
          "ERR initializing string with type '%s'",mxGetClassName(a));

       RENEW(++n); 
       i=mxGetString(a,data,n);

       if (i) wblog(F_L,"ERR failed to read string >%s<",data);
       return *this;
    };

    wbstring& initrep(const char* s, unsigned n) { 
       if (s && s[0] && n) {
         unsigned l=strlen(s); init(l*n+1);
          if (l==1) { memset(data,s[0],n); data[len-1]=0; }
          else {
             unsigned i=0; char *d=data;
             for (; i<n; ++i, d+=l) { strcpy(d,s); }
          }
       }
       else { init(); }

       return *this;
    };

    wbstring& initrep(char c, unsigned n) {
       if (n) {
          if (!isprint(c)) wblog(FL,
             "WRN %s() got non-print char=%c<%d>",FCT,c,c);
          init(n+1); data[len-1]=0; memset(data,c,n);
       }
       else { init(); }
       return *this;
    };

    void set(char c)  {
       if (len) {
       memset(data, c, len-1); data[len-1]=0; }
    };

    wbstring& operator= (const char* s) { return init(s); }

    bool gotTerm0() {
        for (unsigned i=0; i<len; ++i) { if (data[i]==0) return 1; }
        return 0;
    };

    bool strchr(const char *s) {
        if (data && data[0]) {
           if (!s || !s[0]) wblog(FL,"WRN %s() got empty char set !?",FCT);
           else { for (; *s; ++s) { if (::strchr(data,*s)) return 1; }}
        }
        return 0;
    };

    wbstring& push(const char *F, int L, const char* s, char extend=1);
    wbstring& push(const char* s, char extend=1) {
        return push(0,0,s,extend); };

    wbstring& pushf(const char *F, int L, const char *fmt, ...);
    wbstring& cpy(const char *F, int L, const char *s); 

    int RegEx_replace(const char *F, int L, const char *pat, const char *rep,
       char icase=0, char gflag=1);

    int printf(const char *F, int L, const char *fmt, ...);

    wbstring& time() { return time(::time(NULL)); };
    wbstring& time(const time_t &t) { RENEW(32);
      #pragma omp critical (got_CPTR_TIME)
       { ctime_r(&t,data); } 
       return *this;
    };

    wbstring& time_sys() { return time_sys(::time(NULL)); };
    wbstring& time_sys(const time_t &t) { RENEW(32);
      #pragma omp critical (got_CPTR_TIME)
       { strftime(data,31,"%a %b %d %H:%M:%S %Z %Y",localtime(&t)); }
       return *this;
    };

    wbstring& Upper() {
       for (unsigned n=strlen(data), i=0; i<n; i++)
       data[i]=std::toupper(data[i]);
       return *this;
    };

    wbstring& Lower() {
       for (unsigned n=strlen(data), i=0; i<n; i++)
       data[i]=std::tolower(data[i]);
       return *this;
    };

    wbstring toupper() const {
    wbstring sout(*this); sout.Upper(); return sout; };

    wbstring tolower() const {
    wbstring sout(*this); sout.Lower(); return sout; };

    wbstring& Chomp() { 
       if (data && data[0]) { int i=strlen(data);
          for (; i>=0; --i) {
             if (data[i]>20 && data[i]!=' ') { data[i+1]=0; break; }
          }
          if (i<0) data[0]=0;
       }
       return *this;
    };

    char getopt(char c) {
       for (unsigned i=0; data[i]; i++) {
          if (data[i]==c) {
             data[i]=-1; 
             return 1;
          }
       }; return 0;
    };

    char getopt() {
       for (unsigned i=0; data[i]; i++) {
          if (data[i]>0) { return data[i]; } 
       }; return 0;
    }

    const char* basename(const char c='/') const {
       if (!data) return data;
       return Wb::basename(data,c);
    };

    mxArray* toMx() const {
       return mxCreateString(data ? data : "(null)");
    };

    size_t toHash(unsigned offset=0) const;

    void add2MxStruct(mxArray *S, unsigned i) const {
       mxSetCell(S,i,toMx());
    };

    template <class T>
    char* init2Fmt(const T& x __attribute__ ((unused)), int n=-1, int p=-1){
       (*this)=wbstring(Wb::num2Fmt<T>(FL,n,p));
       return data;
    };

    template <class T> wbstring& operator<< (const T* x); 
    template <class T> wbstring& operator<< (const T& x);

    wbstring operator+ (const wbstring &s) const { return (*this)+s.data; };
    wbstring operator+ (const char *s) const {
       wbstring sout; 
       unsigned n1=data? strlen(data):0, n2=s && s[0]? strlen(s):0;
       sout.init(n1+n2);
       if (n1) { strcpy(sout.data, data); }
       if (n2) { strcpy(sout.data+n1, s); }
       return sout;
    };

    wbstring& operator+=(const char *s) {
       if (data && s && s[0]) {
          unsigned n1=strlen(data), n2=strlen(s);
          if (n2) { Resize(n1+n2); strcpy(data+n1, s); }
       }
       else if (s) { (*this)=s; }
       return *this;
    };

    wbstring& operator+=(const wbstring &s) {
       if (s.data) { (*this)+=s.data; }
       return *this;
    };

    unsigned cat(unsigned &k, const char *s) {
       if (s) { unsigned i=0;
          for (; s[i]!=0 && k<len; ++i) { data[k++]=s[i]; }
          if (s[i]) wblog(FL,"ERR wbstring::cat out of bounds (%d/%d)",k,len);
       }
       data[k]=0; return k;
    };

    wbstring toStr(char level=0) const { 
       if (data) {
          if (isPrint()>level) 
               { return *this; }
          else { return wbvector<char>::toStr(); }
       }
       else { return wbstring("(null)"); } 
    };

    char& operator[] (unsigned i) const { return data[i]; };

    unsigned strlen_()  const { 
       return (data && data[0] ? strlen(data) : 0); };

    bool isEmpty() const { return (data && data[0] ? 0 : 1); };

    char isPrint() const {  
        if (data) { unsigned i=0, n=0;
           for (; i<len && data[i]; ++i) {
              if (!isprint(data[i])) { return 0; }}
           for (; i<len; ++i) { if (data[i]) {
               if (isprint(data[i])) { ++n; } else { return 0; }
           }}
           return (n ? 2 : 1);
        }
        return 0;
    };

    bool operator! () const { return (data && data[0] ? 0 : 1); }
    explicit operator bool() const { return (data && data[0] ? 1 : 0); };

    unsigned isName(int L=-1) const {
       unsigned i, n=strlen(data);
       if (L>=0 && L<(int)n) return 0;
       for (i=0; i<n; i++) if (!isalnum(data[i])) return 0;
       return n;
    };

    int operator== (const char *s) const {
        return (strcmp(data, s)==0);
    };
    int operator!= (const char *s) const {
        return (strcmp(data, s)!=0);
    };

    int operator== (const wbstring &s) const {
        return strcmp(data, s.data)==0;
    };
    int operator!= (const wbstring &s) const {
        return strcmp(data, s.data)!=0;
    };
    int operator<  (const wbstring &s) const {
        return strcmp(data, s.data)<0;
    };
    int operator<= (const wbstring &s) const {
        return strcmp(data, s.data)<=0;
    };
    int operator>  (const wbstring &s) const {
        return strcmp(data, s.data)>0;
    };
    int operator>= (const wbstring &s) const {
        return strcmp(data, s.data)>=0;
    };

    int firstLocation (const wbstring &, unsigned k=0);
    int firstLocation (const char*, unsigned offset=0);

    std::istream & getline (std::istream &istr, char c='\n');

  protected:
  private:
};

std::istream & wbstring::getline (std::istream &istr, char c) {
    const unsigned maxlen=256;
    char tmp[maxlen];
    istr.getline(tmp, maxlen-1, c);

    init(strlen(tmp)+1); strcpy (data, tmp);

    return istr;
}

int wbstring::firstLocation (const wbstring &s, unsigned k) {
    return firstLocation(s.data,k);
}

int wbstring::firstLocation (const char *s, unsigned offset) {
    int i, j, n=strlen(data), m=strlen(s);

    for (j=0, i=offset; i<n && j<m; i++, j++) {
        if (data[i] != s[j]) {
            i-=j; 
            j=-1;
        }
    }

    if (j==m) return i-m; 
    else return -1;
}

namespace Wb {
   template<class T> 
   wbstring num2Str(const T &x, const char *fmt=0, unsigned l=0);

   wbstring int2Str(long x, const char *fmt=0, unsigned l=0);
   wbstring char2Str(int x, unsigned l=0);

wbstring repStr(const char* s_, unsigned n) { 
   wbstring s; 
   return s.initrep(s_,n);
};

class String : public wbvector<char> { 
  public:

    String(unsigned n)
     : wbvector<char>(n), ipos(0), fline(NULL) { if (n) { data[n]=0; }};

    String(const char *F, int L, unsigned n)
     : wbvector<char>(n), ipos(0), fline(NULL) { if (n) { data[n]=0; }
       if (F && F[0]) { fline=shortFL(F,L); }
    };

    int pushf(const char *fmt, ...);

  protected:
  private:

    unsigned ipos;
    const char *fline;
};

}; 

#include <queue> 

class wblogbuf_struct {

 public:

   wblogbuf_struct() : file(NULL), line(0), text(NULL) {};

   wblogbuf_struct (const char *s) : file(NULL), line(0), text(NULL) {
      if (s) {
         unsigned n=strlen(s);
         WB_NEW(text,n+1); strcpy(text,s);
      }
   };

   wblogbuf_struct (const char *F, int L, const char *s)
    : file(NULL), line(L), text(NULL) { unsigned i=0,k=0,n;
      if (F) {
         for (; F[i]; i++) { if (F[i]=='/' || F[i]=='\\') k=i+1; }
         n=strlen(F+k); WB_NEW(file,n+1); strcpy(file,F+k);
      }
      if (s) { n=strlen(s); WB_NEW(text,n+1); strcpy(text,s); }
   };

   wblogbuf_struct(const wblogbuf_struct &S)
    : file(S.file), line(S.line), text(S.text
   ){
   };

  ~wblogbuf_struct (){
      if (file) { WB_DELETE(file); }
      if (text) { WB_DELETE(text); }
      line=0;
   };

   void unset() { file=NULL; line=0; text=NULL; };

   char *file; int line; char *text;

 protected:
 private:
};

class wblogbuffer {

  public:

   void push(const char *F, int L, const char *s) {
      wblogbuf_struct S(F,L,s); 
      buf.push(S); S.unset();
   };

   void flush() {
      while (!buf.empty()) {
         wblogbuf_struct S=buf.front();
         wblog(
            S.file ? S.file : __FILE__,
            S.line ? S.line : __LINE__, S.text
         );
         S.unset();
         buf.pop(); 
      }
      doflush(); 
   };

 protected:
 private:

   queue<wblogbuf_struct> buf;
};

  wblogbuffer wblogBuf;

#endif

