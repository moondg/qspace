function [rval,found] = getopt(varargin)
% function [rval,found] = getopt(..)
%
%    General routine to deal with input options
%    either specified as cell (usage #1) or structure (uage #2).
%    The second return argument specifies whether a particular
%    option was found in the input.
%
% Usage #1
%
%    Analyse options provided by a cell vector such as varargin
%    This consists of three steps:
%
%    1) Initialization
%       getopt('init',args)
%       getopt('INIT',args)
%
%    with args = {'opt1',val1,'opt2',val2,...,'-flag',...}.
%    the first usage is case-insensitive on the name of the options (optN)
%    where the second usage with capital 'INIT' is casesensitive.
%
%    2) Queries
%
%       valI=getval('optI',defI)  % with defI a default value for option I, or
%       flagJ=getval('-flagJ')    % checking flags (boolean)
%
%       where flags typically start with '-' or '--').
%       Options can also be grouped (i.e. multiple options referring
%       to the same thing): getopt({opt1,opt2} [,default_value])
%
%    3) Clean up by any of the following:
%
%       getopt('check_error')
%          in case that any position in args went unchecked, 
%          this throws an error.
%
%       a=getopt('get_last',a_def);
%          args may include one value that does not come
%          with a paired name, and that typically also does not
%          represent a flag either, but a value of some sort.
%          This value, irrespective of the place where it occured
%          in between options in args, is returned by this call.
%          If more than one option in args is still left unchecked,
%          this throws an error. If all options have been checked
%          already, this returns the default a_def.
%
%       args=getopt('get_remaining');
%          get all remaining options in argsthat have not yet been checked.
%
% Usage #2 - based on input structure I
%
%    rval = getopt(I,'field1','subfield1',...,default_value)
%
%    This looks for a sequence of fields specified by strings;
%    if the final object exists, its value is returned,
%    otherwise getopt() returns the specified default value.
%
% Wb,Nov05 ; Wb,Jan09,08

  persistent args errcount chkcase

  i=find(cellfun(@ischar,varargin)); if ~isempty(i)
    if ~isempty(regexp([varargin{i}],'dbstop_if')), wbdie(''); end
  end
  if iscell(args)
     i=find(cellfun(@ischar,args)); if ~isempty(i)
        if ~isempty(regexp([args{i}],'dbstop_if')), wbdie(''); end
     end
  end

  if nargin<1, eval(['help ' mfilename]); return; end
  if nargout, rval=[]; end

  if isnumeric(varargin{1})
       lflag=varargin{1}; varargin(1)=[];
  else lflag=0; end

  found=0;

  if isstruct(varargin{1})
   % auxilliary usage to extract option from structure
   % usage: q=getopt(I,'field1','subfield1',...,default_value)
   % NB! if a field has a trailing '?', it is considered
   % optional (e.g. relevant for compatibility of altered
   % info structures) // Wb,Feb22,17
     I=varargin{1}; n=numel(varargin)-1;
     for i=2:n, f=varargin{i};
        if ~ischar(f) || isempty(f), varargin, wbdie('invalid usage'); end
        if f(end)=='?', f=f(1:end-1);
           if isfield(I,f), I=getfield(I,f); end
        elseif isfield(I,varargin{i})
           I=getfield(I,varargin{i});
        else n=n+1; break; end
     end
     if i==n
          rval=I; found=1;
     else rval=varargin{end}; end
     return
  end

  if ~iscell(args)
     if isempty(args), args={};
     else args = {args}; end
  end

  nargs=numel(varargin); o=varargin{1};
  if nargs>2
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  str1=(ischar(o) && ~isempty(o) && o(1)~='-');

  if nargin==2 && str1
     if isequal(o,'init')
        if ~isempty(args)
           wblog('WRN','interfering other concurrent getopt scan !?'); 
           dispstack(1); disp(args)
        end
        args=varargin{2}; errcount=0; chkcase=0;
        return

     elseif isequal(o,'INIT')
        if numel(varargin)>2, wbdie('invalid usage'); end
        args=varargin{2}; errcount=0; chkcase=1;
        return
     end
  end

  if nargs==1 && str1
     if strcmpi(o,'check_error')
        if length(args)>0
           fprintf(1,'\n'); disp(args)
           errcount = errcount + 1;
           wbdie(1,'invalid option(s)');
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
              if fflag>1, dispstack(); end
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
     else wbdie('invalid usage (missing default value)');  end

     return

  else
     NAME=o; narg=length(args);

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
             if i<=length(args) && ...
                isnumeric(args{i}) && isempty(NAME{io}(1)=='-~:!')
                if isscalar(args{i}), wblog('WRN',...
                  'asking for optional flag, yet value %g specified !?',args{i});
                   rval=(args{i}~=0); args(i)=[];
                else wbdie(...
                'asking for optional flag, yet array specified !?'); end
             end
        else rval=0; end
     else
        if found
           if i==narg
              wblog('ERR','value expected for option %s', vname);
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

