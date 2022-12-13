function rval = getopt(varargin)
% Usage #1: rval = getopt(args)
%
%    options can be grouped (i.e. multiple options implying the same)
%    getopt({opt1,opt2} [,default_value])
%
% Usage #2: rval = getopt(I,'field1','subfield1',...,default_value)
%
% Wb,Nov05 ; Wb,Jan09,08

  persistent args errcount chkcase

  i=find(cellfun(@ischar,varargin)); if ~isempty(i)
    if ~isempty(regexp([varargin{i}],'dbstop_if')), error('Wb:ERR','ERR'); end
  end
  if iscell(args)
     i=find(cellfun(@ischar,args)); if ~isempty(i)
        if ~isempty(regexp([args{i}],'dbstop_if')), error('Wb:ERR','ERR'); end
     end
  end

  if nargin<1, eval(['help ' mfilename]); return; end
  if nargout, rval=[]; end

  if isnumeric(varargin{1})
       lflag=varargin{1}; varargin(1)=[];
  else lflag=0; end

  if isstruct(varargin{1})
   % auxilliary usage to extract option from structure
   % usage: q=getopt(I,'field1','subfield1',...,default_value)
   % NB! if a field has a trailing '?', it is considered
   % optional (e.g. relevant for compatibility of altered
   % info structures) // Wb,Feb22,17
     I=varargin{1}; n=numel(varargin)-1;
     for i=2:n, f=varargin{i};
        if ~ischar(f) || isempty(f), varargin, wberr('invalid usage'); end
        if f(end)=='?', f=f(1:end-1);
           if isfield(I,f), I=getfield(I,f); end
        elseif isfield(I,varargin{i})
           I=getfield(I,varargin{i});
        else n=n+1; break; end
     end
     if i==n, rval=I; else rval=varargin{end}; end
     return
  end

  if ~iscell(args)
     if isempty(args), args={};
     else args = {args}; end
  end

  nargs=numel(varargin); o=varargin{1};
  if nargs>2
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if nargin==2
     if isequal(o,'init')
        if ~isempty(args)
           wblog('WRN','overwriting current other getopt scan !?'); 
           q=dbstack(); dispstack(q(2:end)); disp(args)
        end
        args=varargin{2}; errcount=0; chkcase=0;
        return

     elseif isequal(o,'INIT')
        if numel(varargin)>2, wberr('invalid usage'); end
        args=varargin{2}; errcount=0; chkcase=1;
        return
     end
  end

  if nargs==1
     if strcmpi(o,'check_error')
        if length(args)>0
           wblog(2,'ERR','\Ninvalid option(s)\N');
           disp(args)
           errcount = errcount + 1;

           fprintf(1,'Stack: '); s='-> ';

           [stack, index] = dbstack; m=length(stack);
           for i=m:-1:1
              if i==1, s=''; end
              fprintf(1,'%s::%d %s', stack(i).name, stack(i).line, s);
           end

           if length(stack), fprintf(1,'\n\n\n'); end

           error(' ');
        end

        if nargout, rval=errcount; end
        args={}; errcount=0;
        return

     elseif strcmpi(o,'get_remaining')
        rval=args; args={}; errcount=0;
        return

     elseif strcmpi(o,'get_status')
        fflag=numel(find(o<'Z'));
        if fflag
           if ~isempty(args)
              u=getuser(0,'getopt_dbg');
              if isequal(u,args), fflag=-fflag; else
                 setuser(0,'getopt_dbg',args);
              end
           end
        end
        if ~nargout || fflag, S='getopt()';
           if fflag, S=dbstack;
              if numel(S)>1, S=S(2);
                 S=sprintf('%s:%g',regexprep(S.file,'.*\/',''),S.line);
              end
           end
           if ~isempty(args), s={'',S,numel(args),'','',''};
              if s{3}==1, s{4}='entry'; else s{4}='entries'; end
              if numel(errcount), s{5}=num2str(errcount); end
              if numel(chkcase),  s{6}=num2str(chkcase);  end
              if ~fflag, s{1}=char(10); end
              fprintf(1,'%s   %s got %g %s still (e=%s, c=%s)',s{:});
              disp(args)
              if fflag>1, dispstack(dbstack); end
           else
              if fflag, s=''; else s=char(10); end
              fprintf(1,'%s   %s no current entries%s\n',s,S,s);
           end
        end
        return
     end
  end

  if strcmpi(o,'get_last'), n=numel(args);

     if n==1
        rval=args{1}; args={};
     elseif n, getopt('check_error');
     elseif nargs>1, rval=varargin{2};
     else wberr('invalid usage (missing default value)');  end

     return

  else
     NAME=o; narg=length(args); found=0;

     if ~iscell(NAME), NAME={NAME}; end
     no=length(NAME);

     for io=1:no
        if chkcase, vname=NAME{io}; else vname=lower(NAME{io}); end
        for i=1:narg, if ~ischar(args{i}), continue; end
            if chkcase, ai=args{i}; else ai=lower(args{i}); end
            if isequal(vname,ai)
               if chkcase || isequal(NAME{io},args{i})
                    found=1;
               else found=2; end
               break
            end
        end
        if found, break, end
     end

     if nargs<2
        if found, if lflag, wblog(1,' * ','%s',NAME{io}); end
             rval=found;
             args(i)=[];
             if i<=length(args) && isnumeric(args{i}) && NAME{io}(1)~='-'
                if isscalar(args{i}), wblog('WRN',...
                  'asking for optional flag, yet value %g specified !?',args{i});
                   rval=(args{i}~=0); args(i)=[];
                else wberr(...
                'asking for optional flag, yet array specified !?'); end
             end
        else rval=0; end
     else
        if found
           if i==narg
              wblog('ERR - value expected for option %s', vname);
              errcount=errcount+1;
              args(narg)=[];
           else
              rval=args{i+1}; args(i:i+1)=[];
              if lflag
                 if isnumeric(rval)
                      wblog(1,' * ','%-8s: %g',vname,rval);
                 else wblog(1,' * ','%-8s: %s',vname,rval); end
              end
           end
        else rval=varargin{2}; end
     end
  end

end

