function q=getfields(S,varargin)
% function S=getfields(S,  field1, field2, ... )
% function S=getfields(S, {field1, field2, ...})
%
%    get substructure of S in terms of the specified fields.
%
% Wb,Apr10,14

  if nargin<2 || ~isstruct(S)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end
  if numel(varargin)==1 && iscell(varargin{1})
       ff=varargin{1};
  else ff=varargin; end

  nf=numel(ff);
  for i=1:nf
     if ~ischar(ff{i}), wbdie(...
        'invalid usage (expecting strings for field names [#%g])',i+1);
     end
  end

  f0=fieldnames(S); e={};
  for i=1:nf
     if ~isfield(S,ff{i})
         q=regexp(f0,ff{i}); k=[];
         for j=1:numel(q), if ~isempty(q{j}), k(end+1)=j; end, end
         if ~isempty(k)
            ff{i}=f0(k);
         else e(end+1)=ff(i); ff(i)=[];
         end
     else ff{i}={ff{i}};
     end
  end

  ff=cat(1,ff{:})';
  f2=setdiff(ff,f0);

  if ~isempty(f2)  || ~isempty(e)
     s=sprintf(', ''%s''',e{:},f2{:}); wbdie(...
    'asking for invalid field%s {%s}',iff(numel(f2)>1,'s',''),s(3:end));
  end

  fx=setdiff(f0,ff);
  if ~isempty(fx)
       q=rmfield(S,fx);
  else q=S; end

  q=orderfields(q,ff);

end

