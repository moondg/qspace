/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : wblog (logging routines)
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

#ifndef __WB_WBIO_H__
#define __WB_WBIO_H__

/* -------------------------------------------------------------------- */
// tags: tmpfile, tempfile, tempmat, TMPDIR // Wb,Jun06,19

namespace Wb {

wbstring tmpfile(const char *F, int L,
   const char *t=NULL, const char *x=NULL, const char *f=NULL);

class tmpFile { 
  public:
     tmpFile() : fid(NULL), ncall(-1) { };

     tmpFile(const char *F, int L,
        const char *f=NULL, const char *x=NULL, const char *mode=NULL)
      : fid(NULL), ncall(-1) { init(F,L,f,x,mode); };

    ~tmpFile() { if (fid) {
        wblog(PFL,
           "I/O closing tmpfile `%s'",fout.data? fout.data:"(null!?)");
        fprintf(fid,"\n%-18s %5s %s close %s\n\n", 
           shortFL(PFL),"", Wb::TimeStamp('T').data, Wb::TimeStamp('D').data);
        fclose(fid); fid=NULL;
     }};

     void init(const char *F, int L,
        const char *f=NULL, const char *x=NULL, const char *mode=NULL);

     int blogf(const char *F, int L, const char *fmt, ...);

     int sepline(const char *w="-", unsigned len=-1);

     void flush() { if (fid) fflush(fid); };

   FILE *fid;
   int ncall; 
   wbstring fout;

  protected:
  private:
};

}; 

template<class T>
int Str2Idx(       
   const char *F, int L,
   const char* s, wbvector<T> &idx, T offset=0, char cmpct=0
);

template<class T>
int Str2Idx(
   const char* s, wbvector<T> &idx, T offset=0, char cmpct=0
);

template<class T> int Str2Idx(
   const char* s, wbvector<T> &idx, unsigned offset=0, char cmpct=0
){ return Str2Idx(s,idx,T(offset),cmpct); };

wbvector<unsigned> Str2Idx(
   const char *F, int L, const char* s, unsigned offset=0
){
   wbvector<unsigned> I; 
   int e=Str2Idx(F_L,s,I,offset);
   if (e<0) wblog(FL, 
      "ERR %s() invalid '%s' (%d;%d)",FCT,s?s:"null",-e,offset);
   return I;
};

wbvector<unsigned> Str2Idx(const char* s, unsigned offset=0) {
   return Str2Idx(FL,s,offset);
};

template<class T>
int Str2Idx(
   const char* F, int L, 
   wbvector<T> &I, const mxArray *a, T offset=0, char cmpct=0
);

template<class T>
int Str2Idx(
   wbvector<T> &I, const mxArray *a, T offset=0, char cmpct=0
){ return Str2Idx(FL,I,a,offset,cmpct); };

void WbPrintMatrixC (
    const wbMatrix<wbcomplex> &M, const char *istr="", int space=10);

#endif
