#!/usr/local/bin/julia -p 12
#
#    /software/opt/xenial/x86_64/julia/0.4.5/bin/julia -p 12
#    check consistency of XStore data.
#
# Wb,Sep14,16 ; Wb,Aug21,24

# Matlab equivalent $SYM/private/check_XStore.m is much slower than this!
# julia --startup-file=no -p 3 ./check_XStore.jl $RCS SU2

# Examples sequence
# > check_XStore.jl --sym SU5 > check_XStore.log 2>&1 &
# > check_XStore.pl check_XStore.log
# > egrep TASK check_XStore.log | transXStore.pl --sym SU5 > check_XStore_SU5.sge

# using Base
# using Statistics
# import REPL

  using Printf  # required by @printf, etc.
  using Dates   # now(), etc

  me=realpath(@__FILE__);
  # @__FILE__ works with @everywhere, PROGRAM_FILE does not!

  @eval @everywhere push!(LOAD_PATH,replace($me,r"\/+[^\/]*$" => ""));

  me=replace(PROGRAM_FILE, r".*\/" => "");
  @eval @everywhere me=$me;

# dump(LOAD_PATH)

@everywhere using QSModules
@everywhere using RCStore
@everywhere using MAT  # required by all threads

# NB! `@everywhere' is only defined when julia is called with option -p ##
  @everywhere tflag = false;

  if length(ARGS)!=2; error("\n\n"*
     "  ERR 2 arguments requires (got "*repr(length(ARGS))*"):\n"*
     "  ERR usage: $me [RCStore] [RCsym]\n");
  end

  RC_=ARGS[1];
  sym=ARGS[2]; RCS="$RC_/$sym";

  if isdir(RC_) == false
     error("\n\n  ERR $me: invalid directory '$RC_'\n"); end

  m=match(r"^S[UOp]\d+$",sym);
  if m === nothing
     error("\n\n  ERR $me: invalid symmetry '$sym'\n"); end

  if isdir(RCS) == false
     istr=["","invalid RC_STORE directory '"*repHome(RCS)*"'"];
     if match(r"RCStore",RCS) === nothing
     istr[3]="need to include symmetry"; end
     error(join(istr,"\n   ERR "));
  end

  @printf("\n  %s\n  - checking %s\n",me,repHome(RCS));

  m = SubString(sym,1,2);
  rsym = parse(Int,SubString(sym,3));

      if m == "SU"; rsym=rsym-1;
  elseif m == "Sp"; rsym=rsym/2;
  elseif m == "SO"; rsym=floor(rsym/2);
  else error("invalid symmetry $sym"); end

  @eval @everywhere rsym=$rsym;
  @eval @everywhere sym=$sym;
  @eval @everywhere RCS=$RCS;

  @printf("  - having %d workers (%d procs, this being %d)\n",
  nworkers(), nprocs(), Threads.threadid());

# wblog("to be cont'd"); exit(1);

# ------------------------------------------------------------- #
# if not empty, search for these IDs, like ["#d1f9","#997e",...]
  @everywhere xID=[];

  @everywhere ID=zeros(Int,1,length(xID));
  @everywhere for i=1:length(xID);
     ID[i]=parse(Int,replace(xID[i],r"^#",""),16);
  end
# ------------------------------------------------------------- #

# @everywhere Xdir=["++","+-","++-","+--","scalar"];
  @everywhere Xdir=readdir(RCS*"/XStore");

# -------------------------------------------------------------------- #
  tstart=now(); odir=pwd();

  for X in Xdir
    DIR="$RCS/XStore/$X";

    if (!isdir(DIR))
       @printf(stderr,"\n  ERR invalid directory %s\n",repHome(DIR));
       continue
    end

    FX=readdir(DIR);
     # skip debug, test files, etc
       FX = filter(a -> match(r"\.x3d\s*$",a) !== nothing, FX);
    
    nx=length(FX);
    @printf("\n  %s (%g files)\n",repHome(DIR),nx);

    cd(DIR); t1=time(); nq=[0 0 0];

    nq = @distributed (+) for fx in FX
  # for fx in FX
      F=fetch(DIR)*"/"*fx; s=stat(F); n1=[1 0 0]; f="";
      try
        q=matread(F); q=q["x3m"];
        abc=["a","b"]; # including "c" below (be aware of scalars!)

        if isempty(q["c"])
           if isempty(q["x3"]["data"])
             r=q["rtype"]; n1[2]=1;
             if tflag != true
                @printf("=> CTR_ZERO %s\n",fx); # flush(stdout);
             end
           end
        elseif !isempty(q["c"]["qdir"]); push!(abc,"c"); end

      # NB! may have relevant(!) zero-contractions to be updated
      # e.g. having invalid A.ID or B.id // Wb,Sep28,16

        cx=[]; c0=[]; e_=0; err=0; i=0;

        for x in abc; i+=1;
          f=qset2fname(q[x],rsym); if f=="vac"; continue; end
          f="$RCS/CStore/"*f;
          # @printf(" a: %s (%d)\n",f,isfile(f));
          fp=matopen(f); c=read(fp,"CRef"); close(fp);

          cx=q[x]["cid"];
          c0=c["cid"];

          if any(ID.==cx[3])
             @printf("==> FOUND_ID #%04x in %s [%s]\n",cx[3],F,x);
          end
          if any(ID.==c0[3])
             @printf("==> FOUND_ID #%04x in %s\n",c0[3],f);
          end

          if (cx[[1,3]]!=c0[[1,3]])
             err|=(1<<(i-1)); # severe (ctime or ID) mismatch
          end
          if (cx[2]!=c0[2])
             if (cx[2]>c0[2])
              # @printf("  TASK CHECK ERR got newer cstat in X-data !?\n");
                @printf("  TASK RECALC! %s\n",F); # n1[3]=1;
                err|=(1<<(i-1));
             end
             err|=(1<<(2+i)); # only mtime mismatch
          end

          if e_!=err; e_=err;
            if (err&7 != 0)
               @printf("%s\n%3s: [%02d] %s\n  <> [%02d] %s\n",
               repHome(f),x, getOM(c),cstat2str(c0),
               getOM(q[x]),cstat2str(cx));
            end
          end
        end

        if (err&7 != 0) # first 3 bits // err!=0
          es=@sprintf("err = %d '%s'",err,bitstring(convert(Int8,err)));
        # if (err&7 != 0) # first 3 bits
               @printf(" WRN severe cstat mismatch: %s !?\n",es);
        # else @printf(" --> mtime cstat mismatch: %s\n",es);
        # end

          if (contains(F,".x3d")) # Wb,Sep18,16
            f=replace(F,r".*\/" => ""); f=replace(f,r".x3d" => "");
            f=replace(f,r"\(" => " ("); f=replace(f,r"_" => " ");
            f=split(strip(f),r"\s+");
            if length(f)==4
              @printf(" ==> contractRC %s %s #> %s\n", sym,join(f," "),
                 isempty(q["c"]) ? "(empty)" : qset2fname(q["c"],rsym));
            else @printf(stderr,
             "\n  ERR got %s\n-> %s (%g args) !?\n\n",F,f,numel(f));
            end
          end

        # error(es);
          @printf("TASK RECALC%s %s\n\n",(err&7)!=0 ? "!" : " ",F);
          n1[3]=1;
        end

      catch ex

        estr="$ex"; n1[3]=1;
        @printf("ERR %s\nF = '%s'\nf = '%s'\n",
           estr,repHome(F),repHome(f));
        showerror(stdout,ex,catch_backtrace());

        if isa(ex,InterruptException)

           @printf("\n==> got Ctrl-C interrupt (line %d)\n",@__LINE__);
           rethrow(ex);

        elseif isa(ex,ErrorException);

           if match(r"File .*\.cgd.* does not exist",estr) != nothing
            # NB! remove x3d data (pointless to keep when *cgd is empty)
              @printf("TASK REM %s\n",F);
           elseif match(r"\.cgd.* is not a MAT file",estr) != nothing
            # NB! remove x3d data (pointless to keep when *cgd is corrupted)
              @printf("TASK MAT %s\n",F);
           else 
              @printf("TASK CHECK %s\n    error %s\n",F,estr);
           end

        else
           @printf("TASK CHECK %s\n",F);
         # rethrow(ex);
        end

        flush(stdout);
      # continue
      end

      n1 # added to nq(!)
    end
    t2=time();

    if (nq[2]!=0) 
         @printf("     %d/%d zero-files\n",nq[2],nq[1]);
    else @printf("     no zero-files\n"); end

    if (nq[3]!=0)
         @printf("     %d problematic files\n",nq[3]);
    else @printf("     no problematic files\n"); end

    @printf("     at %.3g sec / file and worker\n",(t2-t1)/nx*nworkers());
    flush(stdout);
  end

  cd(odir);

  tfin=now();
  @printf("\n  Started : %s\n",wbnow(tstart));
  @printf(  "  Finished: %s\n",wbnow(tfin));

