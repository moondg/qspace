/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace MPFR routines
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

#ifndef __WB_MPFR_CC__
#define __WB_MPFR_CC__

/**********************************************************************/
   namespace Wb {
/**********************************************************************/

template <unsigned P>
mpfr__<P>::operator wbcomplex() const {
   return wbcomplex(mpfr_get_d(f,DEF_RND),0);
};

template <unsigned P>
mpfr__<P>& mpfr__<P>::init_s(const char *F, int L,
   const char *s, unsigned base
){

   if (!s) wblog(F_L,"ERR mpfr::%s() got null string !?",FCT);
   if (!s[0]) wblog(F_L,
      "ERR mpfr::%s() got empty string (zero data) !?",FCT);

   int e=mpfr_set_str(f,s,base,GMP_RNDD);
   if (e) wblog(F_L, 
      "ERR '%s'\n=> mpfr_set_str() returned %d @ base=%d !?",s,e,base);

   return *this;
};

template <unsigned P>
mpfr__<P>& mpfr__<P>::init_d(
   const char *F, int L, double val, char dcheck
){
   if (::isnan(val) || ::isinf(val))
      wblog(FL,"WRN %s() got val=%g",FCT,val);

   init(F_L,val,'l');

   if (val && dcheck) {
      double x=getval_d(), e=std::fabs((x-val)/val);
      if (e>1e-14) { 
         fprintf(stderr,"\n  %30.20g\n  %30.20g\n\n",val,x);
         wblog(F_L,"ERR %s(%g) conversion error @ %.3g /1e-14",FCT,val,e);
      }
      else if (F && (e || std::fabs(x-int(x)))) {
         wblog(F,L,"TST %s(%g) conversion error @ %.3g",FCT,val,e);
      }
   }
   return *this;
};

template <unsigned P>
mpfr__<P>& mpfr__<P>::init(const char *F, int L, double val,
   char lflag 
){
   int e=0;

   if (::round(val)==val && std::fabs(val)<1E16) {
      mpfr_set_si(f,::round(val),DEF_RND);
      return *this;
   }

   if (std::fabs(val)>1e-4) {
      long p,q; unsigned niter; double r,v0=val; char rflag=0;
      e=Rational(val,p,q,&r,&niter,NULL,6,1000000,1e-10,1e-12);

      if (e || std::fabs((r*q)/p)>1e-14) { val=v0*v0; rflag=1;
         e=Rational(val,p,q,&r,&niter,NULL,6,1000,1e-10,1e-12);
         if (!e && std::fabs((r*q)/p)>1e-14) { e=2; }
      }
      if (!e) {
         mpfr_t b; mpfr_init2(b,P);
         mpfr_set_si(f,p,DEF_RND);
         mpfr_set_si(b,q,DEF_RND); mpfr_div(f,f,b,DEF_RND);
         mpfr_clear(b); 

         if (rflag) {
            mpfr_sqrt(f,f,DEF_RND);
            if (v0<0) mpfr_setsign(f,f,1,DEF_RND);
         }
         return *this;
      }
      val=v0;
   }

   double q=::exp10(6-std::floor(std::log10(std::fabs(val)))), x=val*q;
   if (std::fabs(::round(x)-x)<1e-7) {
      char s[20], l=snprintf(s,20,"%.12g",val);
         if (l>14) wblog(FL,"WRN %s() s=%s",FCT,s);
      init_s(FL,s); return *this;
   }

   if (lflag) {
      e=mpfr_set_d(f,val,DEF_RND);
      if (e) wblog(F_L,"WRN %s() got e=%d for val=%.12g",FCT,e,val);
   }
   else { wblog(F_L,"ERR %s() only accepts int or simple fractions\n"
      "got x=%.16g", FCT,val);
   }

   return *this;
};

template <unsigned P>
wbstring mpfr__<P>::to_str(int base) const {
   wbstring s(prec(base)+7); 
   to_str(s.data,s.len,base); 
   return s;
};

template <unsigned P>
size_t mpfr__<P>::to_str(char *s, size_t n, int base) const {

   size_t l=prec(base)+7;
   mp_exp_t e; 
   if (l>n) wblog(FL,"ERR %s() string out of bounds (%d/%d)",FCT,n,l);

   mpfr_get_str(s,&e,base,0,f,GMP_RNDU); 

   l=strlen(s); if (l<n) {
   l+=snprintf(s+l,n-l," @%d",int(e)); }

   if (l>=n) wblog(FL,
      "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);

   return l;
};

template <unsigned P>
wbstring& mpfr__<P>::toStr(
   wbstring &sout, int base, int nd, char align) const {

   unsigned n=prec(base)+10;
   if (sout.len!=n) sout.init(n); 

   toStr(sout.data,sout.len,base,nd,align);

   return sout;
};

template <unsigned P>
void mpfr__<P>::toStr(
   char *s, unsigned n, int base, int nd, char align) const {

   unsigned i=0, l=3;  char *q, xflag=0;
   mp_exp_t e; 

   s[0]=s[1]=' '; if (nd<0) { xflag=1; nd=0; }

   q=mpfr_get_str(s+l,&e,base,nd,f,GMP_RNDU); 
   if (q!=s+l) wblog(FL,"ERR mpfr_get_str() returned %p",q);

   if (s[l]=='@') { e=0; } 
   else if (!e) {          
      if (s[l]=='-') { s[l-2]=s[l]; l-=1; } else { l-=2; }
      memcpy(s+l,"0.",2); --l;
   }
   else if (e>1 && e<4) {  
       unsigned m=l+e; if (s[l]=='-') { ++m; };
       for (i=l; i<m && s[i]; ++i) { s[i-1]=s[i]; }
       if (i>=m) { s[i-1]='.'; l=m-e-2; e=0; }
       else wblog(FL,"ERR `%s' e%d (%d/%d) !?",FCT,s,e,i,m);
   }
   else { i=(--l); --e; 
      s[i]=s[i+1]; if (s[i]=='-') { ++i;
      s[i]=s[i+1]; } else { --l; }
      s[i+1]='.';
   }

   if (!align) {
      for (l=0; l<n; ++l) { if (s[l]!=' ') { break; }}
   }
   if (l) { 
      for (i=l; i<n; ++i) { s[i-l]=s[i]; if (!s[i]) break; }
      l=i-l;
   }
   else if (e) { l=strlen(s); }

   if (e) { if (l<n) {
      l+=snprintf(s+l,n-l,"%c%d",base<=10 ? 'e':'@', int(e)); }
      if (l>=n) wblog(FL,
         "ERR %s() string out of bounds (%d/%d)\n`%s'",FCT,l,n,s);
   }

   if (xflag) { mpfr__<P> x;
      try { x.init_s(FL,s,base); }
      catch (...) { wblog(FL,"ERR %s() `%s'\n%s",FCT,s,to_str().data); }

      if (mpfr_cmp(x.f,f)) {
         wblog(FL,"%NTST `%s'\n%s\nf=%.16g",s,to_str().data,double(*this));
         fprintf(stderr,"\n   "); mpfr_out_str(stderr,base,0,  f,GMP_RNDD);
         fprintf(stderr,"\n   "); mpfr_out_str(stderr,base,0,x.f,GMP_RNDD);
         fprintf(stderr,"\n");

         wblog(FL,"ERR %s() got base-%d difference",FCT,base);
      }
   }
};

template <unsigned P>
void mpfr__<P>::Cleanup(const char *F, int L) const {

#ifdef LD_CLEBSCH_QS
   if (CG_VERBOSE>5) {
   wblog(F_L," *  free mpfr cache (%s)",myname); }
#endif

   mpfr_free_cache();
};

}; 

#endif

