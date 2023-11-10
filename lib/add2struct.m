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
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  Sout=varargin{1}; k=1; v1=inputname(1);
  if ~isstruct(Sout)
     if ~isempty(v1) || ~isequal(Sout,'-'), k=0; end
     Sout=struct; 
  end

  for i=k+1:length(varargin), n=inputname(i);
      if ~isempty(n)
         Sout=setfield(Sout,n,varargin{i});
      elseif ischar(varargin{i}), n=varargin{i}; nx='';
         opt=0; flag__=0; val__=[];

         l=find(n==':' | n=='=',1);
         if ~isempty(l)
            if l<2 || isempty(regexp(n(1:l-1),'^[A-Za-z][\w_]*$'))
               wbdie('invalid field name ''%s''',n);
            end
            nx=n(l+1:end); n=n(1:l-1);

         elseif n(end)=='?', n=n(1:end-1); opt=1;
         elseif ~isempty(regexp(n,'^[A-Za-z][\w_]*$')), opt=2;
            nx=n;
         else wbdie('invalid expression ''%s''',n);
         end

         cmd=['global flag__ val__; val__=' nx '; flag__=1;'];
         try evalin('caller',cmd); catch; end

         if flag__ || opt==2
            Sout=setfield(Sout,n,val__);
         elseif ~opt
            l=lasterror; inl(1), disp(l.message), inl(1)
            wbdie('invalid expression ''%s''',nx); % wblog('ERR',..)
         end
      else 
         wbdie('failed to assign data (got %s as arg #%d)',...
         class(varargin{i}),i+1);
      end
  end

  clear global flag__ val__

  if k && ~nargout && ~isempty(v1)
     assignin('caller',v1,Sout);
     clear Sout
  end

end

