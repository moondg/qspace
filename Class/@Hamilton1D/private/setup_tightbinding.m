function [HAM]=setup_tightbinding(varargin)
% function [HAM]=setup_tightbinding(varargin)
% Wb,Jan15,15 ; Wb,Aug26,15

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
     tflag= getopt('-t');
     kferm= getopt('-kf');
     rflag= getopt('-rand');
     use_mpo = getopt('--mpo',{'not_specified'});
  if isempty(L)
       L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), wbdie('length L not specified'); end
  if ~iscell(use_mpo)
     wbdie('invalid usage (--mpo requires cell)');
  end

  if NC==1
       sym='Acharge';
  else sym='Acharge,SUNchannel';
  end

  [F,Z,IS]=getLocalSpace('Fermion',sym,'NC',NC,'-v');
  if ~isfield(IS,'istr')
     IS.istr=sprintf('spinless fermions (%s, NC=%g)',sym,NC);
  end

  param=add2struct('-',t,epsd,L,NC,perBC);

  HAM=struct(Hamilton1D);
  HAM.info.istr=sprintf('tight-binding chain of %s',IS.istr);

  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.info.XY=[ (1:L)', zeros(L,1) ];
  if perBC
     HAM.info.XY(2:2:end,2)=1;
  end

  HAM.oez=[
     init_ops(IS.E,'local identity operator (E)')
     init_ops(Z,'fermionic parity (Z)')
  ];

  i='annihilation operator (F)';
  if ~perBC && NC<2 && ~kferm
       HAM.ops=init_ops(F,i);
  else HAM.ops=init_ops(F,i,'-ferm');
  end

  if perBC
     wblog('NB!','using periodic BC (interleaved)');

     HH=[0:L-1; ones(1,L); 2:L+1; ones(1,L); repmat(t,1,L) ]';
     HH(1)=1; HH(end,3)=L;
  else
     HH=[1:L-1; ones(1,L-1); 2:L; ones(1,L-1); repmat(t,1,L-1) ]';
  end

  if rflag
     wblog('NB!','randomizing hopping matrix elements');
     HH(:,end)=t*randn(size(HH,1),1);
  end

  if isequal(use_mpo,{'not_specified'})
     HAM=setup_mpo(HAM,HH);
  else
     HAM.info.HH=HH;
     if isequal(use_mpo,{1}), use_mpo={}; end
     HAM=setup_mpo_full(HAM,use_mpo{:});
  end

  HAM.store='DMRG_tb';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

