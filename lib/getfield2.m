function [A,nd] = getfield2(S,varargin)
% function A = getfield2(S,fieldname1,subfieldname2,... [, {def_value}])
%
%   Similar to getfield but for structure arrays
%   This returns resulting data reshaped into original form.
%
%   If the last argument is a cell with a single entry, it is
%   taken (and returned) as the default value if the specified
%   field does not exist. The field must exist otherwise.
%
% Wb,Nov28,09 ; Wb,Jul31,17

% [Wb,Jul31,17] extended to subfields
% [Wb,May27,25] options '-q', '--def' => replace by trailing {[]} etc.

  got_dval=0; nd=0;

% default value is specified last (for the case that field does not exist)
  if ~isempty(varargin) && iscell(varargin{end}) && numel(varargin{end})==1
     dval=varargin{end}{1}; varargin(end)=[];
     got_dval=1;
  end

  nf=numel(varargin);
  if ~nf
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end
  for i=1:nf
     if ~ischar(varargin{i}), wbdie('invalid usage'); end
  end

  A=repmat({[]},size(S)); nA=numel(A);

  if ~isfield(S,varargin{1})
     if got_dval, nd=nd+1; 
     else wbdie('invalid field ''%s''',varargin{1}); end
  elseif nA==1
     q=getfield(S,varargin{1});
     fld=struct('type','.','subs',varargin);
     for i=2:nf
        if ~isfield(q,fld(i).subs)
           if got_dval, nd=nd+1; q=[]; 
           else wbdie('invalid field ''%s''',sprintf('.%s',fld(1:i).subs)); end
           break
        else
           q=subsref(S,fld(1:i));
        end
     end
     A={q};
  elseif nA==0
  elseif nf==1
     for i=nA:-1:1
        A{i}=getfield(S(i),fld);
     end
  elseif got_dval
     for i=nA:-1:1
        q=getfield(S(i),varargin{1});
        for l=2:nf
           if isfield(q,varargin{l})
                q=getfield(q,varargin{l});
           else nd=nd+1; q=[]; break; end
        end
        A{i}=q;
     end
  else
     fld=struct('type','.','subs',varargin);
     for i=nA:-1:1
        A{i}=subsref(S(i),fld);
     end
  end

  s=ones(nA,1);
  for i=1:nA, s(i)=numel(A{i}); end

  if nA==1
     if ~nd, A=A{1}; else A=dval; end
  else 
     iz=find(s==0); i1=all(s<=1);
     if ~isempty(iz)
        if i1
           if isempty(dval), dval=nan;
           elseif numel(dval)>1, whos dval, wbdie('invalid usage'); end
        end
        A(iz)={dval};
     end
     if i1
        A=reshape([A{:}],size(A));
     end
  end

end

