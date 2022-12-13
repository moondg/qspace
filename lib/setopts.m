function opts=setopts(varargin)
% Function: [opts=]setopts(opts, opt1, 'opt2', ...) 
%
%    Adds options to given option list `opts'; if no output
%    variable is specified and opts is specified on input,
%    output directly replaces opts.
%
% Adding options with values
%
%    opt1, ...     adds opt1 if set in calling routine
%   'opt1:', val   adds opt1 with value specified
%   {'opt1',val}   adds opt1 with value specified
%   'var1?'        adds option var1 with value of variable var1 in caller space
%   {'opt1'}       removes options opt1 from opts if it exists
%
% Adding flags without values
%
%    '--flag?'      adds flag '-flag' to opts iff variable flag exists and is ~=0
%    '--flag'       adds flag '-flag' to opts 
%    '~-flag'       removes flag '-flag' from opts
%
% Other options
%
%    '-v'           verbose mode
%    '-q'           quiet mode
%
% Wb,Oct23,06  Wb,Jun14,07

% NB! see also MatLab's join() routine!

  evalin('caller','global gopt_val__ gopt_flag__');
  global gopt_val__ gopt_flag__

  if ~nargin || iscell(varargin{1}) && nargin<2
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  vn=cell(1,nargin);
  for i=1:nargin, vn{i}=inputname(i); end

  aflag=1;
  vflag=-1; mark=zeros(1,nargin);
  for i=1:nargin, q=varargin{i};
     if ischar(q), q=lower(q);
        if vflag<0
           if     isequal(q,'-q'), vflag=0; mark(i)=1; 
           elseif isequal(q,'-v'), vflag=2; mark(i)=1; end
        end
        if aflag && ~mark(i)
           if isequal(q,'-p'), aflag=0; mark(i)=2; end
        end
     end
  end
  if vflag<0, vflag=1; end

  i=find(mark); if ~isempty(i)
     varargin(i)=[]; vn(i)=[];
  end
  qflag_=2-vflag;

  i=0; noassign=1;  opts={};
  if isequal(varargin{1},'-'), i=1;
  elseif iscell(varargin{1}) && ~isempty(vn{1})
     opts=varargin{1}; noassign=0; i=1;
  end

  n=length(varargin);
  while i<n, i=i+1;

     nm=varargin{i}; gotval=0; qflag=qflag_;
     if iscell(nm) && isempty(vn{i}), l=numel(nm);

        if l==2,     val=nm{2}; nm=nm{1}; gotval=1;
        elseif l==1, val=[];    nm=nm{1}; gotval=-1;
        else wberr('invalid usage {''name'',val }'); end

     elseif i<n && ~iscell(varargin{i+1}) && ...
       ~ischar(varargin{i+1}) && isempty(vn{i+1})

        i=i+1; val=varargin{i}; gotval=1;

     else

        if ~isempty(vn{i}), val=nm; nm=vn{i}; gotval=1; end
     end

     if ~ischar(nm) || isempty(nm) wberr(...
        'invalid usage (option name not of type char)'); end

     if ~gotval
        if nm(end)==':', i=i+1;
           if i>n, wberr('missing value for `%s'')',nm); end
           nm=nm(1:end-1); val=varargin{i}; gotval=1;
        elseif nm(1)=='-' || nm(1)=='~'
           gotval=nm(1); nm=nm(2:end);
           if nm(end)=='?', cflag=1; nm=nm(1:end-1); else cflag=0; end
        elseif isempty(regexp(nm,'^[a-zA-Z]\w*[!?]?$'))
           wberr('invalid usage (invalid option name ''%s'')',nm);
        end
     end

     if ~gotval
        switch nm(end)
           case '?', qflag=1; nm=nm(1:end-1);
           case '!', qflag=2; nm=nm(1:end-1);
           otherwise qflag=0;
        end

        evalin('caller',sprintf([
          'gopt_flag__=0; if exist(''%s'',''var'')==1, ' ...
          'gopt_val__=%s; gopt_flag__=1; end'],nm,nm));
        if gopt_flag__, gotval=1; val=gopt_val__; end

        if ~gotval && qflag+1<qflag_, continue; end
        qflag=max(qflag,qflag_);
     end

     k=findcell(opts,nm);

     if length(k)>1, wblog('WRN',...
     'option `%s'' exists %d times !?',nm,length(k)); end

     if gotval=='-', gopt_flag__=1;
        if cflag, s=nm; s=regexprep(s,'^[-]+','');
           cmd=sprintf(['if ~exist(''%s'',''var'') || ' ...
             'numel(%s)~=1 || ~%s, gopt_flag__=0; end'],s,s,s);
           evalin('caller',cmd);
        end
        if isempty(k) && gopt_flag__
           if aflag
                opts{end+1}=nm;
           else opts={nm, opts{:}};
           end
        end
        continue
     elseif gotval=='~'
        if ~isempty(k)
           if ~qflag, wblog(' * ','removing option `%s''',nm); end
           opts(k)=[];
        end
        continue
     end

     if gotval<=0, if ~gotval && qflag~=1
        wberr('missing value for option ''%s''',nm); end

        if ~isempty(k), opts([k,k+1])=[]; k=[]; end
        continue
     elseif length(k)>1
        opts([k(2:end), k(2:end)+1])=[]; k=k(1);
     end

     if ~isempty(k)
        if vflag>1 && ~isequal(opts{k+1},val)
           wblog('<o>','overwriting option ''%s''',nm);

           if (ischar(opts{k+1}) || numel(opts{k+1})<12) && ...
              (ischar(val) || numel(val)<12)
              s=evalc('disp(opts{k+1})');
              fprintf(1,'  %s\n->%s', s(1:end-2), evalc('disp(val)'))
           end
        end
        opts(k:k+1)={nm,val};
     elseif aflag
        k=length(opts)+1;
        opts(k:k+1)={nm,val};
     else
        opts={nm,val,opts{:}};
     end
  end

  if ~nargout && ~noassign
     if isempty(vn{1}), wberr('invalid usage'); end
     assignin('caller',vn{1},opts); clear opts
  end

  clear global gopt_val__ gopt_flag__

end

