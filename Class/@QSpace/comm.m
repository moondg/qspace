function C=comm(A,varargin)
% function C=comm(A [,B,'conj','+'])
%
%    Compute commutator (or anticommutator, if '+' is specified).
%    if B is ommitted, this assumes B=A && 'conj'.
%    Here 'conj' is interpreted as conjA, i.e. C = {A',B}
%    For anticommutator, use comm(A [, B, cflag], '+');
%
%    The operators A (and B) must be rank-2 or rank-3, where the latter
%    case assumes operator index order (ss',o), i.e. the irop (=spinor)
%    index comes last.
%
% Examples
%
%    comm(F,'+') computes {F',F} // same as acomm(F)
%
% Wb,Feb06,16

  if nargin>1 && ~ischar(varargin{1})
       B=varargin{1}; varargin(1)=[]; conjA=0;
  else B=A; conjA=1; end

  getopt('init',varargin);
     if getopt('conj'); conjA=1; end
  wc=getopt('get_last',-1);

  if ischar(wc)
     if     isequal(wc,'+'), wc=+1;
     elseif isequal(wc,'-'), wc=-1; end
  end
  if ~isequal(wc,-1) && ~isequal(wc,+1), disp(wc)
     error('Wb:ERR','\n   ERR invalid comm/acomm specs');
  end

  n=[numel(A), numel(B)]; C=QSpace(n);
  for i=1:n(1)
  for j=1:n(2), C(i,j)=commAB_1(A(i),B(j),conjA,wc);
  end, end

  if     n(1)==1, C=reshape(C,size(B));
  elseif n(2)==1, C=reshape(C,size(A)); end

end

% -------------------------------------------------------------------- %

function C=commAB_1(A,B,conjA,wc)

  ra=numel(A.Q);
  rb=numel(B.Q);

  if ra<2 || ra>3, error('Wb:ERR',...
     '\n   ERR invalid usage (got A of rank %g !?)',ra); end
  if rb<2 || rb>3, error('Wb:ERR',...
     '\n   ERR invalid usage (got B of rank %g !?)',rb); end

  if ra==rb
     if conjA
        if ra==2
             C=[contractQS(A,'1*', B,'1' ), contractQS(B,'2', A,'2*' )];
        else C=[contractQS(A,'13*',B,'13'), contractQS(B,'23',A,'23*')];
        end
     else
        if ra==2
             C=[contractQS(A,'2', B,'1' ), contractQS(B,'2', A,'1' )];
        else C=[contractQS(A,'23',B,'13'), contractQS(B,'23',A,'13')];
        end
     end
  else
     if conjA
          C=[contractQS(A,'1*',B,'1'), contractQS(B,'2',A,'2*')];
     else C=[contractQS(A,'2', B,'1'), contractQS(B,'2',A,'1' )];
     end
     if ra==2, i=2; else i=1; end
     C(i)=permuteQS(C(i),'132');
  end

  if wc<0
       C=QSpace(C(1))-C(2);
  else C=QSpace(C(1))+C(2);
  end

  if ra~=rb
     if ra>2, q=A.info.otype; else q=B.info.otype; end
     C.info.otype=q;
  end

end

% -------------------------------------------------------------------- %

