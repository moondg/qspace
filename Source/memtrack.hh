/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0 pre-release)
 * Class   : QSpace memory routines
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

#ifndef __WB_MEM_TRACK_HH__
#define __WB_MEM_TRACK_HH__

// Wb,Jul27,12

namespace Wb {
   void MemStat(const char *F=0, int L=0, char lflag=0);
};

#if !( defined __WBDEBUG__ || defined __WB_MEM_CHECK__ )

#define WB_NEW(p,n) Wb::NEW(p,n)
#define WB_NEW_1(p) Wb::NEW_1(p);

#define WB_NEW_2(p,x) { p = new x; \
  if (!p) wblog(__FILE__,__LINE__,"ERR failed to allocate instance"); \
}

#define WB_DELETE(p)   if (p) { delete [] p; p=NULL; }
#define WB_DELETE_1(p) if (p) { delete p; p=NULL; }

namespace Wb {

  template<class T>
  inline void NEW(T* &p, const size_t &n) {

     if (n) {
      { try { p = new T[n]; } catch (...) { p=NULL; } }

        if (!p) { printf("\n");
            MemStat(FL); 
            wblog(FL,"ERR out of memory (%ld * %d = %.1fG) !?",
               n, sizeof(T), n*sizeof(T)/double(1<<30)
            );
        }
     }
     else { p=NULL; }
  };

  template<class T>
  inline void NEW_1(T* &p) { 
     try { p = new T; } catch (...) { p=NULL; }

     if (!p) { printf("\n");
         MemStat(FL); 
         wblog(FL,"ERR failed to allocate instance !?");
     }
  };

}; 

#else

#define WB_NEW(p,n) Wb::NEW(__FILE__,__LINE__,p,n);
#define WB_NEW_1(p) Wb::NEW_1(__FILE__,__LINE__,p);

#define WB_NEW_2(p,x) { p = new x; \
  if (!p) wblog(__FILE__,__LINE__,"ERR failed to allocate instance"); \
}

#define WB_DELETE(p) \
   if (p) { Wb::gML.rm_ptr(__FILE__,__LINE__,p); delete [] p; p=NULL; }

#define WB_DELETE_1(p) \
   if (p) { Wb::gML.rm_ptr(__FILE__,__LINE__,p); delete p; p=NULL; }

#include "memtrack.cc"

#endif 

#endif 

