function A = getfield2(S,varargin)
% function A = getfield2(S,fieldname1,subfieldname2,...)
%
%   similar to getfield
%   but for structure arrays returns resulting data
%   reshaped into original form.
%
% Wb,Nov28,09 ; Wb,Jul31,17

% Wb,Jul31,17 : extended to multiple input subfields
% e.g. perBC=getfield2(HAM,'info','param','perBC','--def',0);

% default value if field does not exist
  x__='__undef__'; xdef=x__; qflag=0;

  if ~isempty(varargin) && iscell(varargin{end}) && numel(varargin{end})==1
     xdef=varargin{end}{1}; varargin(end)=[];
  end

  getopt('init',varargin);
     if getopt('-q'), qflag=1; end
     if isequal(xdef,x__)
        xdef=getopt('--def',x__); if isequal(xdef,x__)
        xdef=getopt('--default',x__); end
     end
  varargin=getopt('get_remaining');
  if ~qflag && ~isequal(xdef,x__), qflag=1; end

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
     if qflag, qflag=qflag+1;
     else wbdie('invalid field ''%s'' !?',varargin{1}); end
  elseif nA==1
     q=getfield(S,varargin{1});
     fld=struct('type','.','subs',varargin);
     for i=2:nf
        if ~isfield(q,fld(i).subs)
           if qflag, qflag=qflag+1; q=[];
           else wbdie(...
             'invalid field ''%s'' !?',sprintf('.%s',fld(1:i).subs));
           end
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
  elseif qflag
     for i=nA:-1:1
        q=getfield(S(i),varargin{1});
        for l=2:nf
           if isfield(q,varargin{l})
              q=getfield(q,varargin{l});
           else qflag=qflag+1; q=[]; break; end
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
     if qflag<=1, A=A{1}; else A=xdef; end
  else 
     iz=find(s==0); i1=all(s<=1);
     if ~isempty(iz)
        if i1
           if isempty(xdef), xdef=nan;
           elseif numel(xdef)>1, whos xdef, wbdie('invalid usage'); end
        end
        A(iz)={xdef};
     end
     if i1
        A=reshape([A{:}],size(A));
     end
  end

end

