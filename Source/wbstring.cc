/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wbstring
 *
 * Copyright 2022 Andreas Weichselbaum
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

#ifndef __WB_STRING_CC__
#define __WB_STRING_CC__

// ----------------------------------------------------------------- //
// tags: filesize fsize kB MB GB

char* Wb::memsize2Str(double x, char *s, unsigned l) {
    if (int(l)<0) {
       if (x<5E3) sprintf(s,"%g bytes",x);        else
       if (x<1E6) sprintf(s,"%.3g kB",x/(1<<10)); else
       if (x<1E9) sprintf(s,"%.3gM",x/(1<<20));   else
                  sprintf(s,"%.3gG",x/(1<<30));
    }
    else {
       unsigned i=0; if (l) {
       if (x<5E3) i=snprintf(s,l,"%g bytes",x);        else
       if (x<1E6) i=snprintf(s,l,"%.3g kB",x/(1<<10)); else
       if (x<1E9) i=snprintf(s,l,"%.3gM",x/(1<<20));   else
                  i=snprintf(s,l,"%.3gG",x/(1<<30));
       }
       if (i>=l) wblog(FL,
         "WRN %s() string out of bounds (%s; %d/%d)",
          FCT,s,strlen(s),l
       );
    }
    return s;
};

wbstring Wb::repHome(const wbstring &file) {
   return Wb::repHome(file.data);
};

namespace Wb {

template<class T>
wbstring num2Str(const T &x, const char *f, unsigned l) {
    wbstring s(int(l)<16 ? 16 : l), fmt; 
    if (!f || f[0]) {
       fmt.init2Fmt(x,-1,4); f=fmt.data; 
    }
    l=snprintf(s.data,s.len,f,x);
    if (l>=s.len) wblog(FL,"ERR %s() string out of bounds "
       "(%d/%d)\n%s '%s' -> '%s'",FCT,l,s.len,TSTR(T),f,s.data);
    return s;
};

template<>
wbstring num2Str(const wbcomplex &x, const char *f,
    unsigned l __attribute__ ((unused))
){ return x.toStr(f && f[0] ? f : "%.4g"); };

#ifdef QS_USING_MPFR
template<>
wbstring num2Str(const Wb::quad &x, const char *f, unsigned l
){ return num2Str(double(x),f,l); }; 
#endif

template<>
wbstring num2Str(const int &x, const char *f, unsigned l) {
    return int2Str(x,f,l);
};
template<>
wbstring num2Str(const unsigned &x, const char *f, unsigned l) {
    return int2Str(x,f,l);
};
template<>
wbstring num2Str(const long &x, const char *f, unsigned l) {
    return int2Str(x,f,l);
};
template<>
wbstring num2Str(const unsigned long &x, const char *f, unsigned l) {
    return int2Str(x,f,l);
};

template<>
wbstring num2Str(const char &x,
    const char *f __attribute__ ((unused)), unsigned l) {
    return char2Str(x,l);
};
template<>
wbstring num2Str(const unsigned char &x,
    const char *f __attribute__ ((unused)), unsigned l) {
    return char2Str(x,l);
};

wbstring int2Str(long x, const char *f, unsigned l) { 
   wbstring s_(l?l:32); 
   unsigned i=0, j,k=0, m=0, n=0; char *s=s_.data;

   l=snprintf(s,s_.len,f && f[0]? f:"%ld",x);
   for (; s[i]; ++i) {
       if (!isdigit(s[i])) { if (s[i]!='-') ++m; else ++k; }
       else {
          for (++n, ++i; s[i]; ++i) { if (isdigit(s[i])) ++n; else break; }
          break;
       }
   }
   if (s[i] || i!=l || !n || k>1) wblog(FL,
      "ERR %s() invalid int `%s' !?",FCT,s);

   i=(n-1)/3; 
   j=l; if (i>m) { j+=(i-m); }
   if (j>=s_.len) wblog(FL,
      "ERR %s() string out of bounds (%ld / %d)",FCT,x,s_.len);

   if (n>3) { s[j]=0; --n;
      for (i=l; i<=l; --n) { if (i && j) {
         s[--j]=s[--i];
         if (j>1 && n<l && (l-i)%3==0) { s[--j]=','; }
      } else { break; }}
      if (j<i) wblog(FL,
         "ERR %s() invalid int %ld `%s' (%d/%d/%d) !?",FCT,x,s,i,j,n);
   }
   return s_;
};

wbstring char2Str(int x, unsigned l) {
   wbstring s(l?l:16); 

   if (isalnum(x)) l=snprintf(s.data,s.len,"'%c'",x); else
   if (isprint(x)) l=snprintf(s.data,s.len,"'%c'",x); 
   else l=snprintf(s.data,s.len,"%d",x);

   if (l>=s.len) wblog(FL,
      "WRN %s() string out of bounds `%s' (%d/%d)",FCT,s.data,l,s.len);
   return s;
};

}; 

int Wb::findRegEx(const char *str, const char *pat, char icase) {

   int q=0; regex_t r;

   if (!str) wblog(FL,"ERR %s() got null input !?",FCT);
   if (!pat || !pat[0]) wblog(FL,
      "ERR %s() got invalid regexp (empty)",FCT);

   if (strstr(pat,"\\d")) wblog(FL,
      "ERR %s() does not support \\d for digits (use [0-9], instead)",FCT);

   q=(REG_EXTENDED | REG_NOSUB); if (icase) { q|=REG_ICASE; }
   if ((q=regcomp(&r,pat, q))) wblog(FL,
      "ERR %s() invalid regexp '%s' (e=%d)",FCT,pat,q);
   if (!regexec(&r,str,0,NULL,0)) { q=1; } 

   regfree(&r); 
   return q;
};

template<class T>
wbstring& wbstring::operator<< (const T& x) {
   wblog(FL,"ERR %s() type `%s' not specialized",FCT,TSTR(x));
   return *this;
};

template<class T>
wbstring& wbstring::operator<< (const T* x) {
   wblog(FL,"ERR %s() type `%s' not specialized",FCT,TSTR(*x));
   return *this;
};

template<>
wbstring& wbstring::operator<< (const char* s) {
   return (*this)+=s;
};
template<>
wbstring& wbstring::operator<< (const char& x) {
   return (*this)+=Wb::char2Str(x);
};

template<>
wbstring& wbstring::operator<< (const int& x) {
   return (*this)+=Wb::int2Str(x);
};
template<>
wbstring& wbstring::operator<< (const unsigned& x) {
   return (*this)+=Wb::int2Str(x);
};
template<>
wbstring& wbstring::operator<< (const long& x) {
   return (*this)+=Wb::int2Str(x);
};
template<>
wbstring& wbstring::operator<< (const unsigned long& x) {
   return (*this)+=Wb::int2Str(x);
};

template<>
wbstring& wbstring::operator<< (const double& x) {
   return (*this)+=Wb::num2Str(x);
};
template<>
wbstring& wbstring::operator<< (const float& x) {
   return (*this)+=Wb::num2Str(x);
};
template<>
wbstring& wbstring::operator<< (const wbcomplex& x) {
   return (*this)+=x.toStr();
};

size_t wbstring::toHash(unsigned offset) const {

   size_t h=0; 

   if (offset) {
      if (!data) wblog(FL,
         "ERR %s() got offset=%d for null string",FCT,offset);
      if (offset>=len) { wblog(FL,
         "WRN %s() offset out of bounds (%d/%d)",FCT,offset,len);
         return h;
      }
      for (unsigned i=0; i<offset; ++i) { if (!data[i]) wblog(FL,
         "ERR %s() string terminating prior to offset "
         "(i=%d/%d, len=%d)\n%s",FCT,i+1,offset,len,data);
      }
   }

   if (data) {
      h=Wb::string_hash(data+offset,len-offset); }
   return h;
};

size_t Wb::string_hash(const char *s, unsigned n,
   size_t hid 
) {
   if (s) { unsigned i=0;
      if (int(n)<0) { 
         for (; s[i]; ++i) { hid = ((hid<<5) + hid) ^ s[i]; }
      }                
      else {
         for (; i<n && s[i]; ++i) { hid = ((hid<<5) + hid) ^ s[i]; }
         if (s[i]) { wblog(FL,
            "WRN %s() string not zero-terminated (n=%d)",FCT,n);
         }
      }
   }
   else if (n) wblog(FL,"ERR %s() got n=%d for null string",FCT,n);
   return hid;
};

inline wbstring& wbstring::push(
    const char *F, int L, const char* s, char extend) {

    if (!s || !s[0]) { return *this; }

    unsigned n1=(data ? strlen(data) : 0), n2=strlen(s);

    if (n1+n2>=len && extend) { 
       unsigned l=n1+n2+1, m=(l%32); if (m) { l=l+(32-m); }
       if (WBLOG_RESIZE) wblog(F_L,
          "TST %s() resizing %d -> %d",FCT,len,l);
       Resize(l);
    }

    if (n1+n2<len) { strcpy(data+n1,s); }
    else if (!len) wblog(F_L,
       "ERR %s() uninitialized string (%d; %p)",FCT,len,data);
    else {
       if (n1+1<len) { strncpy(data+n1,s,len-n1); data[len-1]=0; }
       wblog(F_L, 
         "WRN %s() string out of bounds (%d+%d/%d)%N%N> %s%N",
         FCT,n1,n2,len,data);
    }

    return *this;
};

wbstring& wbstring::pushf(const char *F, int L, const char *fmt, ...) {

    if (fmt && fmt[0]) {
       size_t l, n1=strlen(data), n2=len-n1;
       if (len<=1) wblog(F_L,"ERR %s() got empty string (%d)",FCT,len);
       if (n1>=len) wblog(F_L,
          "ERR %s() got string without null-terminator (%d)!",FCT,len);

       for (int retry=1; retry<=2; ++retry) { 
         #pragma omp critical (got_va_list) 
          { va_list args;
            va_start(args,fmt);
            l=vsnprintf(data+n1,n2,fmt,args);
            va_end(args);
          }
          if (l<n2) { break; }
          if (retry==1) {
             l+=(n1+1); l+=(-l % (l>128 ? 256 : 64));
             if (WBLOG_RESIZE) wblog(F_L,
                "TST %s() resizing %d -> %d",FCT,len,l);
             Resize(n1+l+1); n2=len-n1; 
          }
          else { data[len-1]=0; wblog(F_L,
             "ERR %s() string out of bounds (%d->%d/%d)%N%N> %s%N",
              FCT,n1,n1+l,len,data);
          }
       }
    }
    return *this;
};

int Wb::String::pushf(const char *fmt, ...) {

    if (!fmt || !fmt[0]) { return -1; }
    if (len<=1) wblog(FL,
       "ERR %s got empty string (%d)",fline? fline:FCT, len);
    if (ipos+strlen(fmt)>=len) wblog(FL,
       "ERR %s string out of bounds (%d+%d/%d)",fline? fline:FCT,
       ipos,strlen(fmt),len);

    unsigned ip0=ipos;

   #pragma omp critical (got_va_list) 
    { va_list args;
      va_start(args,fmt);
      ipos+=vsnprintf(data+ipos,len-ipos,fmt,args);
      va_end(args);
    }

    if (ipos>=len) { data[len-1]=0; wblog(FL,
       "ERR %s() string out of bounds (%d->%d/%d)%N%N> %s%N",
       fline? fline:FCT, ip0,ipos,len,data);
    }

    return (ipos-ip0);
};

inline int wbstring::printf(const char *F, int L,
    const char *fmt, ...
 ){
    int l=-1;
    if (fmt && fmt[0]) {
       if (len<=1) wblog(F_L,"ERR %s() got empty string (%d)",FCT,len);

       va_list args;
       va_start(args,fmt);
       l=vsnprintf(data,len,fmt,args);
       va_end(args);

       if (unsigned(l)>=len) { data[len-1]=0; wblog(F_L,
          "ERR %s() string out of bounds (%d/%d)%N%N> %s%N",FCT,l,len,data);
       }
    }
    return l;
};

inline wbstring& wbstring::cpy(const char *F, int L,
    const char *s
 ){
    if (s) {
       size_t l=strlen(s);
       if (l+1>=len) wblog(F_L,
          "ERR %s() string out of bounds (%s; %d/%d)",FCT,s,l,len);
       if (s[0]) strcat(data,s);
       else data[0]=0;
    }
    else wblog(F_L,"ERR %s() got NULL string",FCT);

    return *this;
};

#endif

