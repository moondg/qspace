function [Iout,mat_]=lma(varargin)
% function [Iout,mat]=lma([file-pattern, opts][, {dir2() opts} ])
%
%    load mat-file from local matlab directory (LMA)
%    in case of file not found, the returned mat_ is empty.
%
% Options
%
%   '-s',..  sort by file name (showing all then)
%   '-x',..  exclude pattern
%   '-a'     show all mat files
%   '-n',..  number of records on display
%   '-d',..  specifies directory tag as recognized by cto (default directory: lma)
%   '.'      choose local directory (same as -d .)
%
% Examples
%
%    lma -n 12 . kro*    % search current directory for at max 12 kro*.mat files
%
% Wb,Oct25,05  Wb,May11,07

  pat='*.mat';
  getopt('INIT',varargin);
     nk   =getopt('-n',4);
     aflag=getopt('-a');  if aflag, nk=Inf; end
     lma_=getopt('-d','lma'); % since `lma' is name of this!
     xpat =getopt('-x','');
     sflag=getopt('-s');
     use_current=getopt('.');

     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0;
     elseif getopt('-Q'), vflag=-1; else vflag=1; end

  args=getopt('get_remaining');
  if ~isempty(args) && ~isempty(regexp(args{1},'^-\d'))
       i=str2num(args{1}(2:end)); args(1)=[];
  else i=-1; end

  if ~isempty(args) && iscell(args{end})
       odir=args{end}; args(end)=[];
  else odir={}; end

  if vflag<=1, oc={'-q'}; else oc={}; end

  nargs=numel(args); mat_='';

  prevDir=pwd;

  try if ~use_current, cto(lma_,oc{:}); end

    if nargs, pat=args{1}; args(1)=[];
       if ~isempty(regexp(pat,'^\w+$')) % \w = [a-z_A-Z0-9] (`word')
          pat=[ '*' pat '*.mat' ]; if sflag, sflag=sflag+2; end
       else
          if isempty(regexpi(pat,'\.mat$'))
             pat=[ pat '*.mat' ]; if sflag, sflag=sflag+1; end
          end
          if isempty(regexp(pat,'\/')) && isempty(regexp(pat,'^[A-Z]'))
             pat=[ '*' pat ];     if sflag, sflag=sflag+1; end
          end
       end
       if isempty(regexp(pat,'[\\\[\]\(\)]'))
          pat=regexprep(pat,'\.\.','*');
          if sflag, sflag=sflag+1; end
       end
    end

    ff=dir2(pat,odir{:}); nf=numel(ff);
    if nf
       if vflag>0 && nf>1, wblog(' * ','%2g files found in lma',nf); end
    elseif ~use_current
       [i,p]=system('cto.pl -q -d uma');
       if ~i && ~isempty(p), cd(p);
          ff=dir2(pat,odir{:}); nf=numel(ff);
          if nf, if vflag>=0, 
             if isdesktop>1
                  wesc_={[e '[30;1m' ],[e '[0m']};
             else wesc_={'',''};
             end
             fprintf(1,[ wesc_{1} ...
                '>> NB! using UMA data (no files found in LMA) ' ...
                wesc_{2} '\n']);
             end
          elseif vflag>1, disp(args)
             fprintf(1,'  pwd: %s\n',pwd);
             fprintf(1,'  pat: %s\n\n',pat);
          end
       end
    end

    if isempty(ff)
       if use_current, fprintf(1,'\n'); end
       if vflag>0 || ~nargout
          fprintf(1,'   No files "%s" found (%s)\n\n',pat,repHome(pwd)); end
       cd(prevDir); if nargout, Iout=[]; end
       return
    elseif ~isempty(xpat)
       j=find(arrayfun(@(x) isempty(regexp(x.name,xpat)), ff));
       ff=ff(j); nf=numel(ff);
    end

    if i>0
       if i>nf, wberr('file ID out of bounds (%d/%d)',i,nf); end
       l=nf-i+1;
    else
       l=max([1,nf-nk]):nf;
       if sflag, q={ff.name}; l=1:nf;
          if sflag>1
             p=regexprep(pat,'\*(.*)\*.mat$','$1');
             if ~isempty(regexp(p,'^\d+$')), p=[p '\d*(\.\d+|)']; end
             p=['.*(' p ').*'];
             for i=1:numel(q), q{i}=regexprep(q{i},p,'$1'); end
          end
          [~,is]=sort(q);
       else
          if isfield(ff(1),'datenum') % in this case 'date' *IS* a string!
               [~,is]=sort([ff.datenum]);
          else [~,is]=sort({ff.date}); end
       end

       ff=ff(is(l)); nf=numel(ff);

       for i=1:nf
          ff(i).datestr=datestr(ff(i).date,'mm/dd/yy HH:MM');
       end

       if vflag>0, fprintf(1,'\n    DATE %13s SIZE  FILENAME\n','');
          for i=1:nf % s=num2str2(ff(i).bytes,'--bytes');
             fprintf(1,'%2d  %s %7.1fM  %s\n', nf-i+1, ...
             ff(i).datestr, ff(i).bytes/(2^20), basename(ff(i).name));
          end
          fprintf(1,'\n');
       end

       if nf==1, i=1; else
          i=input(sprintf('>> Select file <#[1]> (0=none) ... '));
          if isempty(i), i=1; end
       end

       if i<1
          fprintf(1,'\n'); cd(prevDir); if nargout, Iout=[]; end
          return
       else i=nf-i+1;
       end
    end

    mat=ff(i).name; mat_=[pwd '/' mat];
    if vflag>0, fprintf(1,'>> load %s\n',mat); end

    if nargout
       Iout=load(mat_,args{:});
       if numel(args)==1, ff=fieldnames(Iout);
          if numel(ff)==1
             Iout=getfield(Iout,ff{1});
          end
       end
       cd(prevDir); return
    end

    cmd=sprintf(',''%s''',mat_,args{:});
    cmd=['load(' cmd(2:end) ');'];

    v='varx__'; setuser(0,v,'');

    evalin('caller', cmd);

    evalin('caller',['if exist(''mat''), setuser(0,''' v ''',mat); end']);
    m2=getuser(0,v,'-rm'); if ~ischar(m2), m2=''; end

    if ~isempty(m2)
       m0=regexprep(regexprep(mat,'.*\/',''),'\.mat$','');
       m2=regexprep(m2, '\.mat$','');
       if isempty(regexp(m2,m0)), q=[m0 ' -> ' m2];
          assignin('caller','mat',q);
          fprintf(1,['\n   NB! %s.m => variable ''mat'' differs' ... 
          ' from mat-file name\n   mat ''%s''\n\n'],mfilename,q);
       end
    end

  catch l

    cd(prevDir); disp(l.stack(1));
    if length(l.stack)>1, disp(l.stack(end)); end
    rethrow(l)

  end

  cd(prevDir);

end

