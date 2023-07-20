function [HAM]=setup_tb_ladder(varargin)
% function [HAM]=setup_tb_ladder(varargin)
% Wb,Aug27,15

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     NC   = getopt('NC',  1 );
     L    = getopt('L',   []);
     t    = getopt('t',  -1 );
     epsd = getopt('epsd',0.);
     perBC= getopt('-perBC');
     rflag= getopt('-r');
     tflag= getopt('-t');
  if isempty(L)
       L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), wbdie('length L not specified'); end

  if NC==1
       sym='Acharge';
  else sym='Acharge,SUNchannel';
  end

  param=add2struct('-',L,t,epsd,NC,perBC);

  [F,Z,IS]=getLocalSpace('Fermion',sym,'NC',NC,'-v');
  if ~isfield(IS,'istr')
     IS.istr=sprintf('spinless fermions (%s, NC=%g)',sym,NC);
  end

  HAM=struct(Hamilton1D);
  HAM.info.istr=sprintf('tight-binding ladder of %s',IS.istr);
  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.oez=[
     init_ops(IS.E,'local identity operator (E)')
     init_ops(Z,'fermionic parity (Z)')
  ];

  HAM.ops=init_ops(F,'annihilation operator (F)','-ferm');

% DMRG site order -------------------------------------------------- %
%
%   1 --- 3 --- 5 --- 7 ---  ...
%   |     |     |     |
%   |     |     |     |
%   2 --- 4 --- 6 --- 8 ---  ...
%
% ------------------------------------------------------------------ %

% XY=[
%   0 1    \   0 0 1 0
%   0 0    /  
%   1 1    \   1 1 1 0
%   1 0    /  
%   ...        ...
% ];

  XY=[ repmat((0:(L-1))',1,2), repmat([1 0],L,1) ];
  XY=reshape(permute(reshape(XY,[L,2,2]),[3 2 1]),2,[])';
  HAM.info.XY=XY;;

  HH=zeros(2*(L-1)+L,5); l=1;

  for i=2:2:2*(L-1)
     HH(l,:)=[ [i-1, 1], [i+1, 1], t]; l=l+1;
     HH(l,:)=[ [i  , 1], [i+2, 1], t]; l=l+1;
     HH(l,:)=[ [i-1, 1], [i  , 1], t]; l=l+1;
  end

  i=2*L;

  if perBC
     HH(l,:)=[ [1  , 1], [i-1, 1], t]; l=l+1;
     HH(l,:)=[ [2  , 1], [i  , 1], t]; l=l+1;
  end
     HH(l,:)=[ [i-1, 1], [i  , 1], t]; l=l+1;

  if rflag
     wblog('NB!','randomizing hopping matrix elements');
     HH(:,end)=t*randn(size(HH,1),1);
  end

  HAM=setup_mpo(HAM,HH);

  HAM.store='DMRG_tbl';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

