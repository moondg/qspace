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
% Examples
%
% setdef('U',0.12, 'epsd','= -U/3', 'Gamma',0.01);
%    since U may be defined through first arguemnt, the expression '= ...'
%    will evaluate the string if variable is not defined yet.
%
% Wb,Jun15,07  Wb,Mar01,08

  evalin('caller','global sdval__ sdflag__');
  global sdval__ sdflag__

  getopt('init',varargin);
     verb=getopt('-v');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if ~narg, return; end
  if narg==1 && isstruct(varargin{1}), S=varargin{1};
     ff=fieldnames(S); n=length(ff); varargin=cell(1,n);
     for i=1:n, varargin{i}={ff{i},getfield(S,ff{i})}; end
     varargin=cat(2,varargin{:}); narg=length(varargin);
  end

  if narg<2 || mod(narg,2)
     eval(['help ' mfilename]); varargin
     if nargin, wbdie('invalid usage'); else return; end
  end

  if nargout && narg<=2
     if narg~=2, wbdie('invalid usage'); end
     n=varargin{1}; rval=varargin{2}; sdflag__=0;

     evalin('caller',sprintf(...
     'if exist(''%s'',''var''), sdval__=%s; sdflag__=1; end',n,n));

     if ~sdflag__ && ischar(rval) && ~isempty(rval) && rval(1)=='='
        evalin('caller',sprintf(...
        'eval(''sdval__=%s;''); sdflag__=1;',rval(2:end)));
     end
     if sdflag__, rval=sdval__;
        if verb, wblog(' * ','%-10s = %g',n,rval); end
     end
  else
     for i=1:2:narg
        n=varargin{i}; v=varargin{i+1}; sdflag__=0;

        evalin('caller',sprintf( ...
        'if exist(''%s'',''var''), sdflag__=1; sdval__=%s; end',n,n));
        if sdflag__, varargin{i+1}=sdval__; continue; end

        if ischar(v) && ~isempty(v) && v(1)=='='
           evalin('caller',sprintf('eval(''sdval__=%s;'');',v(2:end)));
           v=sdval__;
        end

        assignin('caller',n,v);

        if verb, wblog(' * ','%-10s = %g',n,v); end
     end

     if nargout, rval=struct;
        for i=1:2:narg
        rval=setfield(rval,varargin{i},varargin{i+1}); end
     end
  end

  clear global sdval__ sdflag__

end

