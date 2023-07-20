function [HAM]=setup_TLHubbard(varargin)
% function [HAM]=setup_TLHubbard(varargin)
%
%    Setup of triangular lattice Hubbard model (TLU)
%    in YC(W)x(L) geometry.
%
% Options
%
%    'L',...   length of system
%    'W',..    width of system (4)
%    '--XC'    use XC geometry instead of YC (default)
%
%    'mu',..   chemical potential (relative to half-filling; 0)
%    'U',..    Hubbard onsite interaction (relative to half-filling; 10)
%    't3',...  fermionic hopping amplitudes [vertical, horizontal, diagonal] (-1)
%
% Wb,Feb21,21

% wsys='TLHubbard'; L=6; W=6; U=9; tflag=2; tst_Hamilton1D

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     L  = getopt('L',  []);
     W  = getopt('W',   4);
     NC = getopt('NC',  1);
     pBC= getopt('-pBC');
     useXC=getopt('--XC');

     U  = getopt('U', 10.);
     mu = getopt('mu', 0.);
     t3 = getopt('t3',-1.);

     tflag= getopt('-t');

  if isempty(L)
       L=getopt('get_last',12);
  else getopt('check_error'); end

  if useXC && mod(W,2)
     wbdie('invalid usage (XC requires even width)'); end

  sym={'FermionS','Acharge,SU2spin'};

  param=add2struct('-',L,W,NC,t3,U,mu,sym);
  param.system=sprintf('%s%gx%g',iff(useXC,'XC','YC'),W,L);

  [F,Z,Sop,IS]=getLocalSpace(sym{:},'NC',NC,'-v');

  if ~isfield(IS,'istr')
  IS.istr=sprintf('TLHubbard %s',param.system); end
  IS.Sop=Sop;

  HAM=struct(Hamilton1D);

  HAM.info.istr=IS.istr;
  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.oez=[
     init_ops(IS.E,'local identity operator (E)')
     init_ops(Z,'fermionic parity (Z)')
  ];

  Nf=contract(F,'13*',F,'13');
  if mu
       Hloc=-mu*Nf; s={'-mu*n'};
  else Hloc=QSpace; s={}; end

  if U
     Q=skipzeros((Nf-IS.E),'-f');
     Hloc=Hloc + (U/2)*Q*Q;
     s{end+1}='(U/2)*(n-1)^2';
  end

  HAM.ops=init_ops(F,['annihilation operators (F)'],'-ferm');
  if norm(Hloc)>1e-9
     HAM.ops(end+1,1)=init_ops(Hloc,['Hloc = ' strjoin(s,' + ')],'~hconj');
  elseif Hloc
     wblog('WRN','ignoring Hloc @ %.3g',norm(Hloc)); 
     Hloc=QSpace;
  end

% DMRG site order for YC cylinder ---------------------------------- %
% NB! using zig-zag (and not snake pattern)
% snake pattern gets H-terms with range 2W
% whereas zig-zag has H-terms with range W
%
%   :  /   :  /   :  /      /  :
%   W --- 2W  -- 3W  -- ... - L*W
%   :  /   :  /   :  /      /  :
%   : ---- : ---- : --- ... -  :
%   :  /   :  /   :  /      /  :
%   2 --- W+2 -- 2W+2 - ... -  *
%   |  /   |  /   |  /      /  |
%   1 --- W+1 -- 2W+1 - ... - (L-1)*W+1
%
% one diagonal contribution only => equivalent to YC triangular;
% NB! alternating the diagonal order in vertical direction
% switches to the XC triangular setting! // Wb,Jul01,21
% ------------------------------------------------------------------ %
% NB! using the one above (/) has a short-ranged diagonal bond
% around the cylindrical boundary (CB), i.e. (n*W , n*W+1)
% wheras the diagonal terms in the bulk have range W+1.
% Instead, if one had taken the other diagonal contribution (\),
% then the bulk diagonal terms would be shorter (W-1), yet 
% at the price over a very long-rang CB term of distance 2W-1
% => given the present zig-zag order, prefer the triangular
%     seting above using the diagonal "/"
% ------------------------------------------------------------------ %

  q={ repmat(1:L,W,1), repmat(1:W,L,1)' };
  if useXC
     wblog(' * ','using XC cylinder');
     q{1}(1:2:W,:)=q{1}(1:2:W,:)-.5;
  end
  HAM.info.XY=[ q{1}(:), q{2}(:) ];

  SS=cell(1,3*L); N=L*W; l=1;

  H0=sparse(1:N,1:N,repmat(mu/2,1,N)); % (mu/2) since adding H0' below

  if Hloc
     S=repmat(struct('iop',[0 2],'hcpl',1.),1,N);
     for i=1:N, S(i).iop(:,1)=i; end
     SS{l}=S; l=l+1;
  end

  if numel(t3)==1, t3=repmat(t3,1,3); end

  if pBC, wblog('NB!','using full periodic BC (pBC)'); 
       oi={'-pX'}; 
  else oi={}; end

  for i=1:L, I=repmat(i,W,1); J=(1:W)';
     for w=1:3
        switch w
          case 1, I2=[I   I]; J2=[J,J+1]; t=t3(1);
          case 2, I2=[I-1,I]; J2=[J J  ]; t=t3(2);
          case 3, I2=[I-1,I];             t=t3(3);
             if ~useXC,       J2=[J,J+1];
             else J=(2:2:W)'; J2=[J J+1; J J-1];
             end

          otherwise wbdie('invalid switch'); 
        end

        kk=sub2ind_cylinder(I2,J2,L,W,oi{:}); nk=size(kk,1);
        S=repmat(struct('iop',[0 1; 0 1],'hcpl',t),1,nk);
        for j=1:nk, S(j).iop(:,1)=kk(j,:); end
        SS{l}=S; l=l+1;

        H0=H0+sparse(kk(:,1),kk(:,2),t,N,N);
        if i==1 && ~pBC, break; end
     end
  end

  HAM.info.HS=[SS{:}];

  if 1
       HAM=setup_mpo_full(HAM,'nsw',3);
  else HAM.mpo=QSpace(1,L*W);
  end

  HAM.info.mpo.H0f=(H0+H0');

  HAM.store='DMRG_TLU';

  if tflag, keyboard, end

end

