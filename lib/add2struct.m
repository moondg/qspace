function Sout=add2struct(varargin)
% Function: Sout=add2struct([Sin,]var1,var2,..)
%
%   Adds variables to given structure with the fieldname
%   being the variable name (therefore no expressions are
%   allowed as field values.
%
%   NB! input structure Sin may also be an emtpy object,
%   or equivalently '-', in which case a new structure is
%   created (this is required if the first object to be
%   added to Sout is a structure itself).
%
% Variable name may be specified explicitely as
%
%   'var1'   variable with name `var1' must exist
%   'var1?'  variable with name `var1' (only if it exists)
%   'name1:var1[?]'  using `name1' as fieldname for variable `var1'
%
% Wb,Jul09,07

  global flag__ val__

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  Sout=varargin{1}; k=1; v1=inputname(1);
  if ~isstruct(Sout)
     if ~isempty(v1) || ~isequal(Sout,'-'), k=0; end
     Sout=struct; 
  end

  for i=k+1:length(varargin), n=inputname(i);
      if ~isempty(n)
         Sout=setfield(Sout,n,varargin{i});
      elseif ischar(varargin{i}), n=varargin{i};
         if n(end)~='?', opt=0; else n=n(1:end-1); opt=1; end

         l=find(n==':' | n=='=');
         if ~isempty(l), l=l(1);
            nx=n(l+1:end); n=n(1:l-1);
            if n(end)=='?', n=n(1:end-1); opt=1; end
         else nx=n; end

         cmd=sprintf([...
           'global flag__ val__; flag__=0; ' ...
           'try, val__=%s; flag__=1; catch; end '], nx);
         evalin('caller',cmd);

         if flag__
            Sout=setfield(Sout,n,val__);
         elseif ~opt
            wblog('ERR','invalid variable or expression ''%s''',nx);
            l=lasterror; inl(1), disp(l.message), inl(1)
         end
      else 
      wberr('failed to assign data (arg #%d)',i+1); end
  end

  clear global flag__ val__

  if k && ~nargout && ~isempty(v1)
     assignin('caller',v1,Sout);
     clear Sout
  end

end

