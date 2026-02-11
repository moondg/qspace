function rval=setdef(varargin)
% Function: val=setdef(name1, val1, name2, val2, ...)
%
%    If variable <name_i> is not set yet in workspace
%    it is set to default value <val_i>.
%
%    alternatively the input list also can be specified
%    as structure with s.name_i = val_i.
%
%    if return argument is requested with two input arguments
%    a=setdef('b',val); % returns exists('b') ? b : val
%    if more arguments are specified, input arguments are
%    returned as structure.
%
% Options
%
%    -v   verbose mode
%    -e   do not overwrite empty values
%         (by default, empty values are considered uninitialized)
%
% Examples
%
% setdef('U',0.12, 'epsd','= -U/3', 'Gamma',0.01);
%    since U may be defined through first arguemnt, the expression '= ...'
%    evaluates the string (if variable is not defined yet).
%
% Wb,Jun15,07  Wb,Jun16,25

% [Wb-06/16/25] added option -e

  evalin('caller','global sd_val__ sd_gotval__');
  global sd_val__ sd_gotval__

  getopt('init',varargin);
     vflag=getopt('-v');
     eflag=getopt('-e');
  args=getopt('get_remaining'); nargs=length(args);

  if ~nargs, return; end
  if nargs==1 && isstruct(args{1}), S=args{1};
     ff=fieldnames(S); n=length(ff); args=cell(1,n);
     for i=1:n, args{i}={ff{i},getfield(S,ff{i})}; end
     args=cat(2,args{:}); nargs=length(args);
  end

  if nargs<2 || mod(nargs,2)
     eval(['help ' mfilename]); args
     if nargin, wbdie('invalid usage'); else return; end
  end

  if ~nargin || mod(nargs,2), wbdie('invalid usage'); end

  for i=1:2:nargs
     vn=args{i}; v=args{i+1}; sd_gotval__=0;
     if ~ischar(vn), wbdie(...
       'invalid usage (variables must be specified by name)'); end

     evalin('caller',sprintf( ...
       'if exist(''%s'',''var''), sd_gotval__=1; sd_val__=%s; end',vn,vn));

     if sd_gotval__ && (eflag || ~isempty(sd_val__))
        args{i+1}=sd_val__;
        continue
     end

     if ischar(v) && ~isempty(v) && v(1)=='='
        evalin('caller',sprintf('eval(''sd_val__=%s;'');',v(2:end)));
        v=sd_val__;
     end

     if ~nargout
        if vflag, wblog(' * ','%-10s = %g',vn,v); end
        assignin('caller',vn,v);
     end
  end

  if nargout
     if nargin==2, rval=args{2};
     else rval=struct;
        for i=1:2:nargs
        rval=setfield(rval,args{i},args{i+1}); end
     end
  end

  clear global sd_val__ sd_gotval__

end

