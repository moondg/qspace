/* LICENSE/license.txt ; Class: QSpace MEX routine options class */

#ifndef __WB_OPTS_CC__
#define __WB_OPTS_CC__

/* ------------------------------------------------------------- */

void OPTS::init(const mxArray** ain, const unsigned l) {

   wbvec<char> s(128);
   if (MAT.mfp) wblog(FL,
      "ERR OPTS already initialized by file\n%s",MAT.fname.data);

   aa.init(l, (mxArray**)ain);
   name.init(l);

/* extract all strings and store them into name
 * (wheather they are a name or not) - name=NULL otherwise */
   for (unsigned i=0; i<aa.len; ++i) {
      if (!mxIsChar(aa[i])) { continue; } else
      if (mxGetString(aa[i],s.data,s.len)) { wblog(FL,
         "WRN failed to read OPTS string (len<%d)",s.len);
         continue;
      }

      name[i].init(s.data);
      aa[i]=NULL; // unset since "used" as name
   }
};

/* ------------------------------------------------------------- */

#endif
