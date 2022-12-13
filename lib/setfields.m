function s=setfields(s,varargin)
% function s=setfields(s,'field1',val1,'field2',val2,...)
%
%    set several fields in structure s
%
% Alternative usage
% function s=setfields(s,sref)
%
%    set empty yet *existing* fields in s to value found
%    in structure sref (if field exists).
%
% Options
%
%    '-k'  keep current structure, i.e. do not add extra fields
%
% Wb,Apr16,11 ; Wb,Sep08,15

   e=0; if ~isstruct(s)
     if isempty(s) || isequal(s,'-'), s=struct;
     else e=1; end
   end

   if e || nargin<2
      helpthis, if nargin || nargout, wberr('invalid usage'), end
      return
   end

   if isstruct(varargin{1})
      s=set_fields_based_on_reference(s,varargin{:});
      return
   end

   n=numel(varargin); if ~n, return; end

   if mod(n,2), e=2; else
      for i=1:2:n, if ~ischar(varargin{i}), e=3; break; end, end
   end
   if e
      helpthis, if nargin || nargout, wberr('invalid usage'), end
      return
   end

   for i=1:2:n
      s=setfield(s,varargin{i},varargin{i+1});
   end

   if ~nargout
      n=inputname(1);
      if ~isempty(n), assignin('caller',n,s); clear s; end
   end

end

% -------------------------------------------------------------------- %
function S=set_fields_based_on_reference(S,S0,kflag)

   if numel(S)~=1 && numel(S0)~=1, wberr(...
     'invalid usage (single data set required)'); end

   if nargin==3
      if ~isequal(kflag,'-k'), wberr('invalid kflag'); end

      ff=fieldnames(S);
      for i=1:numel(ff), f=ff{i};
         if isempty(getfield(S,f)) && isfield(S0,f)
            S=setfield(S,f,getfield(S0,f));
         end
      end
   else
      ff=fieldnames(S0);
      for i=1:numel(ff), f=ff{i};
          S=setfield(S,f,getfield(S0,f));
      end
   end

end

% -------------------------------------------------------------------- %

