
/* COMMENTS / CHANGE LOG ============================================= *
   e.g. that also allows to use backspace on command line output
   for i=1:999, printfc([zeros(1,10)+8, sprintf('i=%03d',i)]); end

   NB! it seems disp(str) does the same thing
   e.g. it shows color escape sequences!! // Wb,Feb18,19
   NB! same for fprintf(1,'%s',str)
   where str contains char(27) = \e // Wb,Aug02,22

 * =================================================================== */

char USAGE[]=""; // outsourced to printfc.m // Wb,Feb14,19

#include "wblib.h"

void mexFunction(
    int nargout, mxArray *argout[],
    int nargin, const mxArray *argin[]
){ Wb::CleanUp aclu; unsigned w=0; try { 

    MX_CHECK_HELPER_NARGS(0,-1,-1);
    if (nargout) usage(FLINE,"no output arguments available!");
    if (!nargin) return;

    unsigned i=0,j,l,m=0; char kflag=0, tag[4]="";
    mxArray *a;

    str[0]=str[STRLEN]=0; 

    for (; int(i)<nargin; ++i) { if (!mxIsChar(argin[i])) break; }

    if (i) {
       if (mxGetString(argin[0],str,STRLEN)) wblog(FL,
          "ERR %s() fmt string out of bounds (%d) !?",FCT,STRLEN);

       if (!strcmp(str,"-k")) { kflag=1; --i; --nargin; ++argin;
          if (mxGetString(argin[0],str,STRLEN)) wblog(FL,
          "ERR %s() fmt string out of bounds (%d) !?",FCT,STRLEN);
       }
    }

    l=strlen(str);
    if (i>1 && l==3)  { strcpy(tag,str);
       if (!strcmp(str,"ERR")) { w|=(1<<0); } else
       if (!strcmp(str,"WRN")) { w|=(1<<1); } else
       if (!strcmp(str,"TST")) { w|=(1<<2); } else
       if (!strcmp(str,"*  ")) { w|=(1<<3); } else
       if (!strcmp(str," * ")) { w|=(1<<4); } else
       if (!strcmp(str,"  *")) { w|=(1<<5); }
       else {
          for (i=0; i<l; ++i) { if (str[i]=='\\' || str[i]=='%') break; }
          if (i==l) { w|=(1<<6);  }
       }
       if (w) { --i; --nargin; ++argin;
          if (mxGetString(argin[0],str,STRLEN)) wblog(FL,
             "ERR %s() fmt string out of bounds (%d) !?",FCT,STRLEN);
          l=strlen(str);
       }
    }

    if (!i) wblog(FL,
       "ERR %s() invalid usage (missing fmt string)",FCT);
    if (l>511) wblog(FL,
       "ERR %s() fmt string out of bounds (%d/%d) !?",FCT,l,STRLEN);

    for (i=0; i<l; ++i) { if (str[i]=='\\') { ++m; }}
    if (m) {
       mxArray* args[nargin];
       for (j=i+m; i<=l; --i, --j) {
          str[j]=str[i]; if (str[i]=='\\') { str[--j]='\\'; }
       }
       memcpy(args,argin,nargin*sizeof(mxArray*));
       args[0]=mxCreateString(str);
       if (!args[0]) { i=-1; }
       else {
          i=mexCallMATLAB(1,&a,nargin,args,"sprintf");
          mxDestroyArray(args[0]);
       }
    }
    else {
       i=mexCallMATLAB(1,&a,nargin,(mxArray**)argin,"sprintf");
    }

    if (i) {
       wblog(FL,"ERR %s() printf() returned error (e=%d)",FCT,i);
       return;
    }

    if ((i=mxGetString(a,str,STRLEN))) wblog(FL,
      "ERR %s() fmt string out of bounds (%d/%d; e=%d) !?",
       FCT,strlen(str),STRLEN,i);
    l=strlen(str);

       for (j=i=0; i<l; ++i, ++j) {
          if (str[i]!='\\') { str[j]=str[i]; continue; }
          switch (str[++i]) {
             case 'n': str[j]='\n'; break;
             case 't': str[j]='\t'; break;
             case 'r': str[j]='\r'; break;
             case 'b': str[j]='\b'; break;
             case 'e': str[j]='\e'; break;

             default : wblog(FL,
             "WRN unknown escape sequence %c%c<%d>",str[i-1],str[i],str[i]);
          }
       }
       str[j]=0; l=j;

    if (w) { char t[32];
       if ((w&3) && !kflag) {
          if (w&1) sprintf(t,"   \e[31m%s\e[0m ",tag); 
          else     sprintf(t,"   \e[35m%s\e[0m ",tag); 
       }
       else { sprintf(t,"   %s ",tag); }

       if (w&1) { printf("\n"); }

       for (j=i=0; i<l; ++i) {
          if (str[i]=='\n') { j=i+1; } else if (str[i]!=' ') break;
       }; if (j) { str[j-1]=0; printf("%s\n",str); }

       for (i=j; i<l; ++i) { if (str[i]=='\n') {
           str[i]=0; printf("%s%s\n",t,str+j);
           j=i+1;
       }}
       if (j<l) { printf("%s%s",t,str+j); }

       if (w) { printf("\n"); }
    }
    else { printf("%s",str); } 

}  catch (Wb::LogException &e) { ExitMsg(e.istr); }
   catch (...) { ExitMsg("caught exception in printfc()"); }

   if (w&1) {
      printf("\n   Issuing ..."); doflush();
      mexErrMsgIdAndTxt("Wb:MEX:printfc","");
   }
   aclu.Check();
};

