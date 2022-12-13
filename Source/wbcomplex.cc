#ifndef __WB_COMPLEX_CC__
#define __WB_COMPLEX_CC__

/* ---------------------------------------------------------------- */
/* ---------------------------------------------------------------- */
// complex number class
// Wb,Dec20,01

wbstring wbcomplex::toStr(const char *fmt) const {
    wbstring zstr(32), zfmt(16);
    unsigned k,n;

    double in=abs(), rn=0; 
    if (in) { rn=r/in; in=i/in;
       if (std::fabs(rn)<1E-12) rn=0;
       if (std::fabs(in)<1E-12) in=0;
    }
    else { in=i; rn=i; }

    if (!rn) {
       if (in) {
          zfmt.printf(FL,"%si",fmt);
          zstr.printf(FL,zfmt.data,i); return zstr;
       }
       else {
          zstr.printf(FL,fmt,r); return zstr;
       }
    }
    else if (!in) { zstr.printf(FL,fmt,r); return zstr; }

    for (n=strlen(fmt), k=0; k<n; ++k) { if (fmt[k]=='+') break; }
    if (k<n) {
       zfmt.printf(FL,"%s%si",fmt,fmt);
    }
    else {
       zfmt.printf(FL,"%s+%si",fmt,fmt);
       zfmt[n]='%'; zfmt[n+1]='+'; 
    }

    zstr.printf(FL,zfmt.data,r,i);
    return zstr;
};

inline wbstring toStr(const wbcomplex &z, const char* fmt) {
   return z.toStr(fmt);
}

inline wbstring toStr(const double &d, const char* fmt) {

   size_t l=0, n=32; char s[n];
   l=snprintf(s,n,fmt,d);
   if (l>=n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,l,n);
   return s;
};

#endif

