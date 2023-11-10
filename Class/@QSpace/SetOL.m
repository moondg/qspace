function A = SetOL(A,varargin)
% function A=SetOL(A,X,k [,otag])   % usage #1
% function A=SetOL(A,itag[,otag])   % usage #2
%
%    Set operator labels of the type { itag, itag* [,'op*'] }.
%    default operator tag (otag) for rank-3 tensors is 'op'.
%
% Usage #1
%
%    set index labels for operator A based on the label found
%    in tensor X at dimension k.
%
% Wb,Aug04,12

  if nargin<2, wblog('ERR','%s() invalid usage',mfilename); end

  if ~ischar(varargin{1})
     if nargin<3 || nargin>4, wbdie('invalid usage #1'); end

     X=varargin{1}; k=varargin{2};
     if nargin==4, otag=varargin{3}; else otag='op'; end

     if ~isa(X,'QSpace')
        if ~isfield(X,'Q') || ~isfield(X,'info')
           wbdie('invalid 2nd argument (QSpace expected)'); end
        X=QSpace(X);
     elseif ~isfield(struct(X).info,'itags') || isempty(X.info.itags)
        wbdie('invalid 2nd argument (QSpace got no itags)');
     end

     [r,t]=gotITags(X);

     if ~r, wblog('ERR','got no itags in reference QSpace (X)');
     elseif k>r, wblog('ERR','index out of bounds (%d/%d)',k,r);
     else t=t{k}; end

  else
     if nargin<2 || nargin>3, wbdie('invalid usage #2'); end

     t=varargin{1};
     if nargin==3, otag=varargin{2}; else otag='op'; end
  end

 if ~ischar(otag), otag, wbdie('invalid otag'); end

  c='*';
  if t(end)==c
     t=t(1:end-1);
     if isempty(t), wblog('ERR','invalid itag (empty)'); end
  end

  t0={ t, [t c]};

  for k=1:numel(A), t=t0;
     n=numel(A(k).Q);
     if n==3, t{3}=[otag c];
     elseif n<2 || n>3, wbdie('invalid operator of rank-%d',n);
     end

     A(k).info.itags=t;
  end

% DON'T! may use SetOP() on command line for testing!
% if ~nargout
%    n=inputname(1); if ~isempty(n)
%       assignin('caller',n,A);
%       clear A
%    end
% end

end

