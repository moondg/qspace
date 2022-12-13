function [HAM]=setup_tightbinding(varargin)
% function [HAM]=setup_tightbinding(varargin)
% Wb,Jan15,15 ; Wb,Aug26,15

  if nargin<1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     L    = getopt('L',   []);
     t    = getopt('t',  -1 );
     epsd = getopt('epsd',0.);
     perBC= getopt('-perBC');
     tflag= getopt('-t');
     rflag= getopt('-rand');
  if isempty(L)
       L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), error('Wb:ERR',...
    '\n   ERR length L not specified'); end

  sym='Acharge'; NC=1;
  [F,Z,IS]=getLocalSpace('Fermion',sym,'NC',NC,'-v');

  if ~isfield(IS,'istr')
     IS.istr=sprintf('spinless fermions (%s, NC=%g)',sym,NC);
  end

  HAM=struct(Hamilton1D);
  HAM.info.istr=sprintf('tight-binding chain of %s (with site permute)',IS.istr);
  HAM.info.IS=IS;

  HAM.oez=[
     init_ops(IS.E,'local identity operator (E)')
     init_ops(Z,'fermionic parity (Z)')
  ];

  i='annihilation operator (F)';
  HAM.ops=init_ops(F,i,'-ferm');

  HH=zeros(L-1,5);

  for i=2:L
     HH(i-1,:)=[ [i-1, 1], [i, 1], t];
  end

  if perBC
     HH(L,:)=[[1, 1], [L 1], t];
  end

  if rflag
     wblog('NB!','randomizing hopping matrix elements');
     HH(:,end)=t*randn(size(HH,1),1);
  end

  if 1
   % permute site order
   % P=1:L;
   % P=L:-1:1;
   % l=fix(L/2); P=[2:l, 1, l+1:L]; % move 1 site to center
   % l=fix(L/2); P=[l, 1:l-1, l+1:L]; % move center site to front

     P=randperm(L);

     wblog('TST','applying random permutation to site order!'); 
     HAM.info.P=P;
  end

  HH(:,1)=P(HH(:,1));
  HH(:,3)=P(HH(:,3));

  HAM=setup_mpo(HAM,HH);

  HAM.store='DMRG_tb';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

