# Wb,Sep19,16

module RCStore

using MAT  # eg. required for getID() below
using Printf
using QSModules

export cstat2str, fname2nby, getID, getOM, qset2fname

# -------------------------------------------------------------------- #
# Wb,Sep19,16

"""
cstat2str(q.cid) - convert cgd.cid(1:4) into more verbose string
"""
function cstat2str(q,w=3)

 # first two bits of w set which output is included
 # 3rd bit uses more space (trying to align data from multiplet calls)
   if (w&1!=0)
      if w&4!=0 # see i jul regarding non-const format string!
         s=@sprintf("%20s (%10s) %-6s",
         time2str(q[1]), sec2str(q[2]-q[1]), @sprintf("#%04x",q[3]));
      else
         s=@sprintf("%20s (%s) %s",
         time2str(q[1]), sec2str(q[2]-q[1]), @sprintf("#%04x",q[3]));
      end
   end

   if (w&2!=0)
      istr=[ "default / empty"                # CGD_DEFAULT,  //  0
             "abelian"                        # CGD_ABELIAN,  //  1  
             "identity"                       # CGD_IDENTITY, //  2 
             "initialized to ref"             # CGD_REF_INIT, //  3
             "initialized to bare size"       # CGD_BSZ_INIT, //  4
             "computed via cgw*CRef"          # CGD_FROM_CGR, //  5
             "computed via contraction"       # CGD_FROM_CTR, //  6
             "computed via decomposition"     # CGD_FROM_DEC, //  7
             "computed via tensorprod (std3)" # CGD_STD3,     //  8
             "1j symbol via std3"             # CGD_1JSY_ST3, //  9 
             "1j symbol via get1J_gen"        # CGD_1JSY_GEN, // 10
             "rank3 obtained via 1j+perm"     # CGD_STD3_X,   // 11
             "obtained via 1j+perm"           # CGD_GEN_X,    // 12
             "explicitely generatec (vac)"    # CGD_EXPLICIT  // 13
      ];

      k=convert(Int,q[4]);
      if k<0 || k>=length(istr) || q[4]!=round(q[4])
         error(@sprintf("ERR invalid cstat type %g !?",q[4]));
      end

      if (w&1!=0); s*"   "*istr[k+1]; 
      else istr[k+1]; end
   end
end

# -------------------------------------------------------------------- #
# Wb,Sep28,16

"""
fname2nby() - get number of boxes from qlables in given filename assuming SU(N)
"""
function fname2nby(f)

    if match(r"^\.",f) !== nothing
       @printf(STDERR,"\n  ERR invalid file '%s'\n",repHome(f));
       return [];
    end

    f=replace(replace(f,r".*\/" => ""),r"\.[\w\d_-]+$" => "");
    f=replace(f,r"_[\d]+\**\s*" => " "); # relevant for X-files
    f=replace(f, r"[(),;\*\s]+" => " ");
    f=replace(replace(f,r"^\s+" => ""),r"\s+$" => "");

    qq=split(f,r"\s+"); m=-1; nb=zeros(size(qq));
    for k = 1:length(qq)
       q=split(qq[k],""); x=zeros(size(q))
       if m>=0
          if m!=length(q); error(@sprintf(
             "\n  inconsistent length of qset (%s %d/%d) !?\n  %s",
             "$q",length(q),m,"$f"));
          end
       else m=length(q);
       end
       for i=1:length(q) # COMPACT_QLABELS
           x[i]=q[i][1]-'0'; if x[i]>9; x[i]-=7; end
       end
       q = x'*(1:length(x));
       nb[k]=q[1];
    end
    nb;
end

# -------------------------------------------------------------------- #
# Wb,Sep26,16

"""
getID(fname) - get ID from mat file for given fname
"""
function getID(f)
   if !contains(f,"cgd")
      error("\n  ERR got invalid CData file '$f'\n\n");
   end
   fid=matopen(f);
   R=read(fid,"CRef"); close(fid);
   R["cid"][3];
end

# -------------------------------------------------------------------- #
# Wb,Sep20,16

"""
getOM(c) - get dimension of outer multiplicity
"""
function getOM(c)
   r=length(c["qdir"]);
   l=length(c["cgd"]["S"]); d=1;
   if l==r+1
      d=c["cgd"]["S"][l];
   elseif l!=r
    # NB! l=r=0 assumes scalard => d=1.
      error(@sprintf("\n  Wb:ERR invalid CData or CRef (r=%d/%d) !?",r,l));
   end
   convert(Int,d);
end

# -------------------------------------------------------------------- #

"""
qset2fname() - convert CRef into more verbose string
"""
function qset2fname(Q,rsym)
    q=copy(Q["qset"]); d=Q["qdir"]; r=length(d); 

    if !isa(q,Array) # q=[q];
     # NB! may have q=[0] for |vac> // Wb,Aug21,24
       if (q!=0); show(Q); flush(stdout);
         error(@sprintf("ERR invalid CRef"));
       end
       return "vac";
    end

    for i=1:length(q) # COMPACT_QLABELS
       q[i] = (q[i]>9 ? q[i]+7 : q[i]) + Int('0'); # '9' + 8 = 'A'
    end
    ss=[];

    q=reshape(convert(Array{Char},q),rsym,r);

    for i=1:r
       push!(ss,join(q[:,i]));
       if i<r; push!(ss, d[i+1]!=d[i] ? ";" : ","); end
    end

    d * "/(" * join(ss) * ").cgd";
end

# -------------------------------------------------------------------- #

end # of Module RCStore

