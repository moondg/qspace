function [P,iP]=initperm(r,varargin)
% function [P,iP]=initperm(r [,opts])
%
%    Initialize permutation of rank r, where the instructions
%    via the options are applied in the sequence they appear.
%    Second return argument, if requested, is the corresponding
%    inverse permutation.
%
% Options
%
%    '--2end',i      move index / indizes i to the end
%    '--2front',i    move index / indizes i to the front
%    '--rot',i       rotate indices by i to the right for i>0 (left for i<0)
%    '--2pos',{i,j}  move index i -> j
%
%    '--last2',i     move last index to position i (same as: '--2pos',{-0,i})
%    '--first2',i    move first index to position i (same as: '--2pos',{1,i})
%
%    '--t2', [i j]   transposition of the two indices i,j
%    '--r3', [i j k] rotate the three indices forward by one
%    '--r3i',[i j k] rotate the three indices backward by one
%
% Wb,Oct07,19

  nargs=numel(varargin);
  if mod(nargs,2), wberr('invalid usage'); end

  P=1:r; e=0;

  for k=1:2:nargs
     q=varargin{k}; j=varargin{k+1};
     if ~ischar(q) || ~isequal(q(1),'-'), wberr('invalid usage'); end

     if     isequal(q,'--2end'  ), i=1:r; i(j)=[]; P=P([i,j]);
     elseif isequal(q,'--2front'), i=1:r; i(j)=[]; P=P([j,i]); 
     elseif isequal(q,'--last2' ), P=[1:j-1,r,j:r-1];
     elseif isequal(q,'--first2'), P=[2:j,1,j+1:r];
     elseif isequal(q,'--rot') && numel(j)==1, j=mod(j,r);
        if     j>0, P=P([r-j+1:r, 1:r-j]); % rotate `to the right'
        elseif j<0, P=P([1-j:r,1:-j]);     % rotate `to the left'
        end
     elseif isequal(q,'--2pos') && numel(j)==2
        if iscell(j)
             l=j{2}; j=j{1}; if numel(l)~=1, e=e+10; continue; end
        else l=j(2); j=j(1);
        end

        iend=find(j<=0); if ~isempty(iend), j(iend)=r-j(iend); end

        i=1:r; i(j)=[]; P=P([i(1:l-1), j, i(l:end)]);

     elseif isequal(q,'--t2')  && numel(j)==2, P(j)=P(j([2 1  ]));
     elseif isequal(q,'--r3')  && numel(j)==3, P(j)=P(j([3 1 2]));
     elseif isequal(q,'--r3i') && numel(j)==3, P(j)=P(j([2 3 1]));
     else e=e+1; end
  end

  if e, wberr('invalid usage (k=%g)',k); end

  if ~isequal(sort(P),1:r)
     wberr('invalid permutation'); end
  if nargout>1, iP(P)=1:r; end

end

