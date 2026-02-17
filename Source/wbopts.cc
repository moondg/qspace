/* ---------------------------------------------------------------------
 * Project : QSpace tensor library (v4.0)
 * Class   : QSpace MEX routine options class
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

#ifndef __WB_OPTS_CC__
#define __WB_OPTS_CC__

/* ------------------------------------------------------------- */

void OPTS::init(const mxArray** ain, const unsigned l) {

   wbvec<char> s(128);
   if (MAT.mfp) wblog(FL,
      "ERR OPTS already initialized by file\n%s",MAT.fname.data);

   aa.init(l, (mxArray**)ain);
   name.init(l);

   for (unsigned i=0; i<aa.len; ++i) {
      if (!mxIsChar(aa[i])) { continue; } else
      if (mxGetString(aa[i],s.data,s.len)) { wblog(FL,
         "WRN failed to read OPTS string (len<%d)",s.len);
         continue;
      }

      name[i].init(s.data);
      aa[i]=NULL; 
   }
};

#endif
