function FF=mfull(FF,varargin)
% function FF=mfull(FF [,opts])
% Wb,Jul11,11

  if nargin<1 || nargout>1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     flds=getopt('fields',{});
  getopt('check_error');

  if isempty(FF), return; end

  if iscell(FF)

     if isnumeric(FF{1}), nF=numel(FF);
        for i=1:nF
           FF{i}=full(FF{i});
        end
     else wbdie('invalid usage'); end

  elseif isstruct(FF)

     nf=numel(flds); nF=numel(FF);
     if ~nf, wbdie('invalid usage'); end
     for i=1:nf
        if ~isfield(FF,flds{i}), wbdie('invalid field ''%s''',flds{i}); end
     end

     for i=1:nF
        for j=1:nf, fj=flds{j}; q=getfield(FF(i),fj);
           if iscell(q), nq=numel(q);
              for k=1:nq, q{k}=full(q{k}); end
           else q=full(q); end
           FF(i)=setfield(FF(i),fj,q);
        end
     end

  else wbdie('invalid usage'); end

end

