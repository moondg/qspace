# load as follows: using QSModules
# Wb,Sep14,16

module QSModules

using Printf
using Dates

export repHome, cto, fsize2str, sec2str, time2str, 
    wbnow, wbstamp, wblog, @sayhello
# , readmat, @matread

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

macro sayhello()
   return :( println("Hello, world!") )
end

# macro matread(F,v)
#    return quote
#       local f=matopen(F); #> ERROR: UndefVarError: matopen not defined
#       local f=MAT.matopen(F); #> ERROR: UndefVarError: MAT not defined
#       local a=read(f,v); close(f);
#       a
#    end
# end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

"""
cto() - change to directory based on cto.pl
"""
function cto(s="")

   d=match(r"^\s*$",s);
   if d==nothing
        d=readall(`cto.pl -d $s`); chomp(d);
   else d=homedir();
   end
   cd(d); 

   s=repHome(d); h=gethostname();
   s=string("\033]0;Julia[",base(36,getpid()),"] $h:$s\007");
   print(STDOUT,s);
 # s=nothing
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# Wb,Sep19,16

# # methods(MAT.matread) is in namespace MAT [after using MAT]
# # methods(readmat)
# 
# """
# readmat(F,vname) - read single variable from mat file (requires using MAT)
# """
# function readmat(F::AbstractString,v::AbstractString)
#    f=matopen(F); a=read(f,v);
#    close(f); a;
# end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

"""
repHome() - shorten path strings
"""

function repHome(P...)

 # unpack tuple into array (otherwise can't assign to P[i])
   P=[P...];

   H=ENV["HOME"];

   for i in 1:length(P)
      P[i] = replace(P[i],Regex(H) => "~");
   end

   if length(P)==1
        return P[1];
   else return P;
   end
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# Wb,Sep29,16

"""
fsize2str() convert given number bytes to string
"""
function fsize2str(s)
   if     s<(1<<10); s=@sprintf("%g",s);
   elseif s<(1<<20); s=@sprintf("%.3gk",s/(1<<10));
   elseif s<(1<<30); s=@sprintf("%.3gM",s/(1<<20));
   elseif s<(1<<40); s=@sprintf("%.3gG",s/(1<<30));
   else;             s=@sprintf("%.3gT",s/(1<<40));
   end
   s;
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# Wb,Sep14,16

"""
sec2str() convert given number seconds to string
"""
function sec2str(t)
   if t<100; s=@sprintf("%.3g",t);
   elseif t<3600
      s=t%60; m=(t-s)/60;
      if s==round(s)
           s=@sprintf("%02d:%02d",  m,s);
      else s=@sprintf("%02d:%02.3f",m,s);
      end
   else
      s=t%60; t=(t-s)/60;
      m=t%60; t=(t-m)/60;
      h=t%24; d=(t-h)/24;

      s=@sprintf("%02d:%02d:%02d",h,m,s);
      if d!=0; s=@sprintf("%d-",d)*s; end
   end
   s;
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# Wb,Sep14,16

"""
time2str() like Libc.strftime, but add split-second resolution
"""
function time2str(t=time())
   s=Libc.strftime("%D %T",t);
   f=@sprintf("%.3f",t-floor(t)); f=replace(f,r"^0\." => ".");
   s*f;
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #
# Wb,Sep14,16

"""
wbnow() - current time in string format "Wed, 14 Sep 2016 21:24:07"
"""
function wbnow(t=now())
 # @printf("%s",Dates.format(t,Dates.RFC1123Format));
                Dates.format(t,Dates.RFC1123Format);
end

# -------------------------------------------------------------------- #

"""
wbstamp() - current time stamp in string format "Wb,Sep14,16"
"""
function wbstamp(t=now())
   Dates.format(t,"Wb,ud,yy")
end

# -------------------------------------------------------------------- #
# Wb,Aug21,24

"""
wblog() 
   general log output 
   with leading file:line segment of location in source caller
"""
function wblog(fmt::String,args...)
   S=stacktrace()[2];
   FL=@sprintf("%s:%d",replace(String(S.file),r".*\/" => ""),S.line);

   Fmt="\e[38;5;243m%-20s\e[0m "*fmt*"\n";
   if match(r"(to be cont[^\s]*d)",fmt) !== nothing
      Fmt="\n$Fmt\n";
   end

   @eval @printf($Fmt,$FL,$(args...))
end

# -------------------------------------------------------------------- #
# -------------------------------------------------------------------- #

end

