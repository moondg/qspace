/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace clock routines
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

#ifndef __WB_CLOCK_HH__
#define __WB_CLOCK_HH__

#include <set>

namespace Wb {
   void pause(const char *F, int L, double tsec); 
   void pause(double tsec) { pause(0,0,tsec); }
   wbstring sec2Str(double t);
};

#define CHRONO_NOW  std::chrono::high_resolution_clock::now

#define CHRONO_CLOCK std::chrono::time_point<std::chrono::high_resolution_clock>
#define CHRONO_DIFF2SEC(A) (1E-9 * \
double(std::chrono::duration_cast<std::chrono::nanoseconds>(A).count()))

#define WBCLK_USE_TAG 1

namespace Wb {

void get_Clock_name(wbstring &name, const char *s, char use_tag);

class ClockSet { 

  public:

   ~ClockSet() { init(); }

    void init(char vflag=1);  

    int info_u (unsigned u, char rflag=1);
    int reset_u(unsigned u, char iflag=0);

    int reset(const char *istr, char uset_tag=1);

    unsigned info(char vflag=0) const; 

    Clock* use(const char *istr, char mode); 

    Clock* use(const char *F, int L, char mode) {
       return use(shortFL(F,L),mode); 
    };

    Clock* get(const char *s, char use_tag=1);

    Clock* insert(Wb::Clock* clk, char lflag=0);

    int erase(Wb::Clock* clk); 

    int add2Mx(MXPut &Iout) const; 

    map<std::string, Wb::Clock*> buf; 

  protected:
  private:
};

   ClockSet Clocks; 

class Clock { 

  public:

    Clock()
     : ttot(0), ctot(0), tneg(0), tref(), cref(0), nz(0), ncall(0),
       clk_flags(0), user(0), gcs(NULL) {};

    Clock(const char *F, int L, char use_tag=1, unsigned u=0)
     : Clock( L>0? shortFL(F,L) : F, L>0 ? use_tag : L, u, &Clocks) { };

    Clock(const char *istr, char use_tag=0)
     : Clock(istr,use_tag,0,&Clocks) { };

    Clock(const char *istr, char use_tag, unsigned u,
       ClockSet *gcs_=NULL, 
       Clock **g=NULL,      
       char rflag=1
     ) : Clock() { user=u;
       if (istr) {
          if (use_tag) { clk_flags|=WBCLK_USE_TAG; }

          init(istr,use_tag); if (rflag) { resume(); }
          gcs=gcs_;   if (!gcs && g) { gcs=&Wb::Clocks; }
          if (gcs) {
             Clock *x=gcs->use(name.data,0); 
             if (g) { (*g)=x; }
          }
       }
       else if (g || gcs) wblog(FL,"ERR %s() invalid usage",FCT);
    };

   ~Clock() {  done(); } 

    void done() {
       if (ncall) { stop(); } 
       if (gcs) {
          Clock *g = gcs->use(name.data,0);
          g->Add(*this); 
       }
       else if (WBLOG_CLK) { info(); }
       init();
    }

    void init() {
       ttot=ctot=tneg=0; tref=CHRONO_CLOCK();
       nz=ncall=clk_flags=user=0; cref=0; gcs=NULL;
       name.init();
    };

    void reset() {
       ttot=ctot=tneg=0; tref=CHRONO_CLOCK();
       nz=ncall=0; cref=0;
    };

    void init(const char *s, char use_tag=1) {
       reset();
       get_Clock_name(name,s,use_tag);
    };

    void start() {
       cref=clock(); ctot=ttot=tneg=0; nz=0; ncall=1;
       tref=CHRONO_NOW();
    };

    void resume() {
       if (!ncall) { start(); } else
       if ((ncall%2)==0) { 
          cref=clock(); ++ncall;
          tref=CHRONO_NOW();
       }
    };

    int stop(const char *F=NULL, int L=0) {
       if (ncall%2) { double tw,tc; ++ncall; 
          ttot+=(tw=get_dt_wall());
          ctot+=(tc=double(clock()-cref)/CLOCKS_PER_SEC);

          if (tw<=0 || tc<=0) { ++nz; if (tw) tneg-=tw;  }
          return 1;
       }
       else {
          if (!ncall) wblog(F_L,
             "WRN clock `%s' not yet started",name.data? name.data:"");
          return 0;
       }
    };

    double gettime(char cflag=0) const; 

    void Add(const Clock &t) { 
       if (ncall%2) { wblog(FL,"WRN %s() clock still running",FCT); }
       ttot += t.ttot;
       ctot += t.ctot; ncall+=2;
    };

    int Switch(const char *istr, Clock **g=NULL, char rflag=1) {
       if (!name.data || !gcs) wblog(FL,
          "ERR %s() invalid usage (not yet initialized)",FCT);
       int q=stop(); Clock *x;

       reset();

       init(istr,clk_flags&1); if (rflag) { resume(); }

       x=gcs->use(name.data,0); 
       if (g) { (*g)=x; }

       return q;
    };

    mxArray* toMx (const char *F=NULL, int L=0) const;
    mxArray* toMxS(const char *F=NULL, int L=0) const;

    void info(const char *istr=NULL, char vflag=0) const;

    double ttot, ctot, tneg;

    CHRONO_CLOCK tref;
    clock_t  cref;

    size_t nz;  

    size_t ncall;

    wbstring name;

    unsigned clk_flags;

    unsigned user;

    ClockSet *gcs;

  protected:
  private:

    double get_dt_wall() const {
       return ( CHRONO_DIFF2SEC( CHRONO_NOW() - tref) );
    };

    double getTickFreq() const;
};

}; 

double Wb::Clock::gettime(char cflag) const {

   double dt=0; 

   if (cflag) {
      dt=ctot; if (ncall%2) { 
      dt+=double(clock()-cref)/CLOCKS_PER_SEC; }
   }
   else {
      dt=ttot; if (ncall%2) {
      dt+=get_dt_wall(); }
   }

   return dt;
};

void Wb::Clock::info(const char *istr, char vflag) const {

   static int first_call=1;

   if (vflag>8) {
      if (vflag=='v') { vflag=1; } else
      if (vflag=='V') { vflag=2; } else
      wblog(FL,"WRN %s() unexpected vflag=%s",FCT,cSTR(vflag));
   }
   if (!ncall && vflag<2) { return; }

   unsigned l=0, n=16; char s[n];
   double tc=gettime('c'), tw=gettime();

   if (nz) {
      double p0=1-nz/(double)(ncall/2); 
      l=snprintf(s,16,"/%5.1f%%", 100*p0);
   }  else s[0]=0;

   if (clk_flags && l<n) {
   l+=snprintf(s+l,n-l,"/%d", clk_flags); } 

   if (first_call) { first_call=0;
      PRINTF("\n   %-28s Count   CPU-time  Wall-time     / call\n",
      "Wb::Clocks"); 
   }

   PRINTF(" %c %-26s%8ld %10s %10s  %9.3g  %s\n",
      (ncall%2) ? '*':' ', istr? istr:name.data, ncall/2, 
      sec2Str(tc).data, sec2Str(tw).data,  
      tw/(ncall/2), s
   );
};

double Wb::Clock::getTickFreq() const {

   unsigned i,j,m=512,n=32; 
   clock_t c1,c2,nc=0;

   CHRONO_CLOCK t1,t2;
   double dt=0;

   for (i=j=0; i<n; ++i) {
      t1=CHRONO_NOW(); c1=clock(); while (c1==clock()) { if (++j>m) break; };
      t2=CHRONO_NOW(); c2=clock();

      if (i) {
         nc+=(c2-c1);
         dt+=CHRONO_DIFF2SEC( t2-t1 );
      }
   }

   return (dt/nc); 
};

namespace Wb {
   class Profile;

   class line_stat__ {
     public:
        line_stat__() : dt(0), aux(0), ncall(0) {};

        void add(double dt_, double dx) { dt+=dt_; aux+=dx; ++ncall; };

        double dt, aux;
        size_t ncall;

     protected:
     private:
   };

   void save_and_clear_Profiling();
};

   map< std::string, map<int, Wb::line_stat__ > > gwb_Profs;

namespace Wb {
class Profile { 
  public:

    Profile(const char *F, int L, const char *istr=0) {
       unsigned i=0, l=0, n=32; char fl[32];
       if (!F || !F[0] || !L) wblog(FL,"ERR %s() invalid empty F !?",FCT);

       for (; F[i]; ++i) { if (F[i]=='/') { l=i+1; }}
       for (i=0; F[l] && i<n; ++l, ++i) {
          if (isalpha(F[l])) { fl[i]=F[l]; } else
          if ((F[l]=='.' || F[l]=='_') && i) { fl[i]='_'; } else { fl[i]='X'; }
       }; l=i;

       if (l<n) { l+=snprintf(fl+l,n-l,"_%03d",L); }
       if (l<n && istr && istr[0]) { l+=snprintf(fl+l,n-l,":%s",istr); }
       if (l>=n) { fl[n-1]=0; wblog(FL,
          "ERR %s() string out of bounds (%d/%d)\n%s",FCT,l,n,fl);
       }
       std::string fl_=fl;

      #pragma omp critical (gwb_Profs__)
       { auto im=gwb_Profs.find(fl_);
         if (im==gwb_Profs.end()) {
            gwb_Profs[fl_]=map<int, Wb::line_stat__ >();
            if ((im=gwb_Profs.find(fl_))==gwb_Profs.end()) 
            { wblog(FL,"ERR %s() map entry failed !?",FCT); }
         }
         M=&(im->second);
       }
       tic();
    };

    void tic() { tref=CHRONO_NOW(); };
    void toc(int L, double x=0) {
       double dt=CHRONO_DIFF2SEC( CHRONO_NOW() - tref );

      #pragma omp critical (gwb_ProfL__)
       { auto im=M->find(L);
         if (im==M->end()) { (*M)[L]=Wb::line_stat__();
            if ((im=M->find(L))==M->end()) 
            { wblog(FL,"ERR %s() map entry failed !?",FCT); }
         }
         im->second.add(dt,x);
       }
    };

    CHRONO_CLOCK tref;

  protected:
  private:

    map<int, Wb::line_stat__> *M; 
};

}; 

#endif

