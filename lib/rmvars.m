function rmvars(fls,varargin)
% function rmvars(mat,var1,var2,...)
% 
%    Remove variables from specified mat file(s).
%    In case of multiple files, fls is assumed to be a structure
%    with field names 'name' and 'folder', e.g. as returned by dir().
%
%    The input names var1, var2, ... are assumed to use
%    regexp() type regular expressions. Vars that start (end) with
%    a letter, are assumed to indeed start (end) with that letter.
%
% Options
%
%    -q     quiet mode
%    -v     more verbose mode
%    -t     test mode
%    ~p     do not preserve time stamp (by default, time stamp is preserved)
%
% Wb,Mar09,22

% tags: rm_variables, rm_vars, save2

  if nargin<2, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  if ischar(fls)
     fls=struct('name',fls,'folder','');
  elseif ~isstruct(fls) || ~isfield(fls,'name')
     wbdie('invalid usage');
  end

  nf=numel(fls); nerr=0;
  if ~nf, fprintf(1,'\n  No files specified.\n\n'); return; end
  for k=1:nf, f=fls(k).name;
     if isempty(regexpi(f,'\.mat$'))
        f=[f '.mat']; fls(k).name=f;
     end
     F=fls(k).folder; if isempty(F), F=f; else F=[F '/' f]; end
     if exist(F)~=2
        nerr=nerr+1; if nerr==1, fprintf(1,'\n'); end
        fprintf(1,'  ERR invalid file %s\n',F);
     end
  end
  if nerr, wbdie('file(s) not found'); end

  getopt('INIT',varargin);
     if     getopt('-v'), vflag=2;
     elseif getopt('-V'), vflag=3;
     elseif getopt('-q'), vflag=0; else vflag=1; end

     if     getopt('-t'), tflag=1;
     elseif getopt('-T'), tflag=2; else tflag=0; end

     pflag=~getopt('~p');
  args=getopt('get_remaining'); nargs=numel(args);

  for i=1:numel(args)
     if regexp(args{i},'^[-~!]'), wbdie('invalid option %s',args{i}); end
     if regexp(args{i},'^[A-Za-z]'), args{i}=['\<' args{i}]; end
     if regexp(args{i},'[A-Za-z]$'), args{i}=[args{i} '\>']; end
  end
  vpat=strjoin(args,'|');

  WSZ='MATLAB:save:sizeTooBigForMATFile';
  wsz=regexprep(WSZ,'.*:','');
  v73='-v7.3';

  if vflag>2, fprintf(1,'\n  using regex pattern: %s\n\n',vpat);
  elseif vflag, fprintf(1,'\n'); end

  tstr='...'; if tflag, tstr='TST'; end

  for k=1:nf, f=fls(k).name;
    F=fls(k).folder; if isempty(F), F=f; else F=[F '/' f]; end
    v=whos('-file',F);
    ik=find( cellfun(@(x) isempty(regexp(x,vpat)), {v.name}));
    ix=1:numel(v); ix(ik)=[]; nx=numel(ix);

    if vflag
       s=fls(k).folder; if ~isempty(s), s=['[' repHome(s) '/] ']; end
       s={s,''}; if ~nx, s{2}=' (no matching variables found)'; end
       fprintf(1,'  %2d/%d %s%s %s\n',k,nf,s{1},f,s{2});
    end
    if ~nx, continue; end

    if vflag
       s={'',sprintf(' %s',v(ix).name)}; if nx~=1, s{1}='s'; end
       b=[ sum([v(ix).bytes]), sum([v.bytes]) ];
       fprintf('%6s removing %d var%s:%s @ %s / %s = %.1f%%\n',tstr,...
       nx, s{:}, num2str2(b(1),'-b'), num2str2(b(2),'-b'),...
       100*b(1)/b(2) );
    end
    if tflag==1, fprintf(1,'\n'); continue; else lastwarn(''); end

    X=load(F); F2=[F '~'];

    warning('off',WSZ);
    save(F2,'-struct','X',v(ik).name);
    warning('on',WSZ);

    [msg,id]=lastwarn;
    if ~isempty(findstr(id,wsz))
        if vflag, fprintf(1,...
        '%6s large variables encountered; using %s\n','...',v73); end
        save(F2,'-struct','X',v(ik).name,v73);
    end
    if pflag
       cmd=sprintf('touch -r ''%s'' ''%s''',F,F2);
       system(cmd);
    end
    if tflag
       fprintf(1,'%6s see %s\n',tstr,F2);
    else
       cmd='/bin/mv -f'; if vflag>2, cmd=[cmd 'v']; end
       system(sprintf('%s ''%s'' ''%s''',cmd,F2,F));
    end
  end

  if vflag && tflag~=1, fprintf(1,'\n'); end

end

