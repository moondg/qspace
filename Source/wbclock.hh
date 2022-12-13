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

namespace  Wb {

   class Clock;

   void pause(const char *F, int L, double tsec); 
   void pause(double tsec) { pause(0,0,tsec); }
   wbstring sec2Str(double t);

};

set<Wb::Clock*> gwb_Clocks;

#define CHRONO_NOW  std::chrono::high_resolution_clock::now

#define CHRONO_CLOCK std::chrono::time_point<std::chrono::high_resolution_clock>
#define CHRONO_DIFF2SEC(A) (1E-9 * \
double(std::chrono::duration_cast<std::chrono::nanoseconds>(A).count()))

namespace  Wb {

class Clock { 

  public:

    Clock(const char *s=0)
     : ttot(0), ctot(0), tneg(0), tref(), cref(0), nz(0), ncall(0)
     { init(s); gwb_Clocks.insert(this); };

   ~Clock() {
       if (ncall) stop(); 
       if (WBLOG_CLCK) info(); 
       name.init(); gwb_Clocks.erase(this);
    };

    void Init() { info(); name.init(); }

    void init(const char *s=NULL) { char s_[32]; reset();
       #ifdef PROG_TAG
         if (s) { snprintf(s_,32,"%.16s::%s",PROG_TAG,s); name=s_; }
         else { name = PROG_TAG; }
       #elif defined(myname)
         if (s) { snprintf(s_,32,"%.16s::%s",myname,s); name=s_; }
         else { name = myname; }
       #else
         name = s? s:"";
       #endif
    };

    void reset() { ttot=ctot=tneg=0; cref=0; nz=ncall=0; } 

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

    mxArray* toMx (const char *F=NULL, int L=0);
    mxArray* toMxS(const char *F=NULL, int L=0);

    void info(const char *istr=NULL) const;

    double ttot, ctot, tneg;

    CHRONO_CLOCK tref;
    clock_t  cref;

    size_t nz, ncall;

    wbstring name;

  protected:
  private:

    double get_dt_wall() const {
       return ( CHRONO_DIFF2SEC( CHRONO_NOW() - tref) );
    };

    double getTickFreq() const;
};

class Clock_resume { 
  public:

    Clock_resume(Clock *t_) : t(t_) { t->resume(); };

   ~Clock_resume() { if (t) t->stop(); };

    int stop() {
       int r=0; if (t) { r=t->stop(); t=NULL; }
       return r;
    };

    int Switch(Clock *t_) {
       int r=0; if (t) r=t->stop(); t=t_; t->resume();
       return r;
    };

    Clock *t;

  protected:
  private:
};

}; 

unsigned showAllClocks() {
   set<Wb::Clock*>::iterator i=gwb_Clocks.begin();
   for (; i!=gwb_Clocks.end(); ++i) { (*i)->info(); }
   return gwb_Clocks.size();
};

unsigned initAllClocks() {
   set<Wb::Clock*>::iterator i=gwb_Clocks.begin();
   for (; i!=gwb_Clocks.end(); ++i) { (*i)->info(); (*i)->reset(); }
   return gwb_Clocks.size();
};

namespace Wb {

class UseClock { 
  public:

    UseClock(Wb::Clock *c) : C(c) { if (C) C->resume(); };
   ~UseClock() { if (C) C->stop(); };

    void done() { if (C) { C->stop(); C=0; }};
    void stop() { if (C) C->stop(); };
    void resume() { if (C) C->resume(); };

  protected:
  private:

    Wb::Clock *C;
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

void Wb::Clock::info(const char *istr) const {

   static int first_call=1;

   if (!ncall) return;

   char s[16];
   double tc=gettime('c'), tw=gettime();

   if (nz) {
      double p0=1-nz/(double)(ncall/2);
      snprintf(s,16,"/%5.1f%%", 100*p0);
   }  else s[0]=0;

   if (first_call) { first_call=0;
      PRINTF("\n  %-30s Count   CPU-time  Wall-time     / call\n",
      "Wb::Clocks"); 
   }

   PRINTF("%c %-26s%10ld %10s %10s  %9.3g  %s\n",
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

