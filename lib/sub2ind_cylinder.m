function k=sub2ind_cylinder(I,J,L,W,varargin)
% function k=sub2ind_cylinder(I,J,L,W [,opts])
%
%    Convert regular lattice indices (i,j) to k in zigzag order
%    where indices outside LxW square patch are mapped back
%    into the lattice assuming periodic boundary conditions.
%
% Options
%
%    --snake   map indices insite LxW square patch to snake order
%    -pX       use periodic BC in x-direction (default: open BC)
%    -oY       use open BC in y-direction (default: periodic BC, i.e. cylinder)
%
% Wb,Feb21,21

% outsourced / adapted from setup_HeisenbergSL_chiral.m

  getopt('init',varargin);
     snake  = getopt('--snake');
     perBCx = getopt('-pX');
     perBCy =~getopt('-oY');
  getopt('check_error');

  n=[numel(I), numel(J)];
  if any(~n), wberr('invalid usage (empty index)'); end

  if any(n>1)
     if     n(1)==1, I=repmat(I,size(J));
     elseif n(2)==1, J=repmat(J,size(I)); end
  end
  sI=size(I); if ~isequal(sI,size(J)), wberr('invalid usage'); end

  I=reshape(I,[],1);
  if perBCx, I=mod(I-1,L)+1;
  elseif any(I<1 | I>L)
     wberr('index I out of bounds ([%g .. %g]/%g)',min(I),max(I),L);
  end

  J=reshape(J,[],1);
  if perBCy, J=mod(J-1,W)+1;
  elseif any(J<1 | J>W)
     wberr('index J out of bounds ([%g .. %g]/%g)',min(J),max(J),W);
  end

  if snake
     i=find(mod(I,2)==0);
     if ~isempty(i), J(i)=W+1-J(i); end
  end

  k=reshape(sub2ind([W L],J,I),sI);

end

