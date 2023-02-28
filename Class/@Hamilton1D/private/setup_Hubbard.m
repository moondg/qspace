function [HAM]=setup_Hubbard(varargin)
% function [HAM]=setup_Hubbard(varargin)
%
%    'L',...   length of system
%    'Ly',..   width of system (1)
%
%    'U',..    Hubbard onsite interaction [H=U*n(n-1)/2] (0.5)
%    'mu',..   chemical potential (=> epsilon_i = -mu; 0)
%              (alternative option: --ph)
%
%    't',...   hopping amplitude (-1)
%    'tx',..   hopping amplitude in x-direction (tx)
%    'ty',..   hopping amplitude in y-direction (ty)
%
%    'NC',..   number of channels (spinless, unless '--spin'; default: 1)
%    '--spin'  use spinfull fermions (by default: spinless fermions)
%    '--perBC' periodic boundary condition
%              (if Ly<=2, in x-, otherwise y-direction)
%    '--ph'    choose particle/hole symmetric sector
%              (only valid if 'mu' is not specified)
%
% Wb,Feb06,18

  if nargin<1
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  getopt('init',varargin);
     Lx=getopt('L', []);
     Ly=getopt('Ly', 1);
     NC=getopt('NC',  1);

     U =getopt('U', 1/2);
     mu=getopt('mu',[]);
     if isempty(mu)
        if getopt('--ph'), mu=U/2;
        else mu=0; end
     end

     t =getopt('t', -1);
     tx=getopt('tx',[]); if isempty(tx), tx=t; end
     ty=getopt('ty',[]); if isempty(ty), ty=t; end

     perBC= getopt('--perBC');
     if     getopt('--spin' ), Sflag=+1;
     elseif getopt('--spinB'), Sflag=-1; else Sflag=0; end
     u1flag=getopt('--U1charge');

     tflag= getopt('-t');

  if isempty(Lx)
     Lx=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(Lx), wberr('length L not specified'); end

  sym={'Fermion',''}; if ~u1flag && mu, u1flag=1; end

  if NC==1,      sym{2}='Acharge';
  elseif u1flag, sym{2}='Acharge(:)';
  else           sym{2}='Acharge,SUNchannel'; end

  if ~Sflag && U && NC<2, Sflag=1; end

  if Sflag
     sym{1}(end+1)='S';
     if Sflag>0
          sym{2}=['SU2spin,' sym{2}];
     else sym{2}=['Aspin,'   sym{2}]; end
  end

  perBCx=0; perBCy=0;

  if perBC
     if Ly<=2
        wblog(' * ','using interleaved setup for perBC-x (Ly=%g)',Ly);
        perBCx=1;
     else
        wblog(' * ','using perBC-y (Ly=%g)',Ly);
        perBCy=1;
     end
  end

  param=add2struct('-',Lx,Ly,t,tx,ty,U,mu,NC,sym,perBC);
  param.perBCxy=[perBCx, perBCy];

  [F,Z,IS]=getLocalSpace(sym{:},'NC',NC,'-v');
  if ~isfield(IS,'istr')
     if Sflag, s='spinfull'; else s='spinless'; end
     IS.istr=sprintf('%s fermions (%s @ %s, NC=%g)',s,sym{:},NC);
  end

  HAM=struct(Hamilton1D);

  s={'',''};
     if U, s{1}='Hubbard'; else s{1}='tight-binding'; end
     if     Ly==1,  s{2}=sprintf('L=%g chain',Lx);
     elseif Ly==2,  s{2}=sprintf('L=%g ladder',Lx);
     elseif perBCy, s{2}=sprintf('L=%gx%g cylinder',Lx,Ly);
     else           s{2}=sprintf('L=%gx%g stripe',Lx,Ly); end
  HAM.info.istr=[s{1} ' ' s{2} ' of ' IS.istr];

  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.oez=[
     init_ops(IS.E,'local identity operator (E)')
     init_ops(Z,'fermionic parity (Z)')
  ];

  X=sum(F);
  if numel(X.Q)>2
       nd=contract(X,'13*',X,'13');
  else nd=contract(X,'1*', X,'1' );
  end

  HAM.ops={init_ops(nd,'local occupation (n=f''f)','~hconj')};

  if numel(F)>1
       s=sprintf('sum^%g F',numel(F));
  else s='F'; end

  HAM.ops{end+1}=init_ops(X,['annihilation operator (' s ')'],'-ferm');

  if U
     Q=0.5*skipzeros(nd*(nd-IS.E),'-f');
     if ~Q, wbdie('got empty charge interaction having U=%g',U); end
     HAM.ops{end+1}=init_ops(Q,'n(n-1)/2','~hconj');
  end

  HAM.ops=cat(1,HAM.ops{:});

% DMRG site order -------------------------------------------------- %
%
%   1 --- Ly+1 -- 2Ly+1 -- ... --  *
%   |       |       |              |
%   |       |       |              |
%   2 ----Ly+2 -- 2Ly+2 -- ... --  *
%   |       |       |              |
%   :       :       :              :
%   : ----- : ----- : ---- ... --  :
%   :       :       :              :
%   |       |       |              |
%   Ly --- 2Ly --- 3Ly --- ... --Lx*Ly 
%
% ------------------------------------------------------------------ %

  XY={ repmat(0:Lx-1,   Ly,1)
       repmat(Ly-1:-1:0,Lx,1)
  };

  if perBCx
     XY{2}(2:2:end,:)=XY{2}(2:2:end,:)+Ly;
  elseif perBCy
     XY{1}=XY{1}+repmat(0.5*sin((pi/Ly)*linspace(0,Ly,Ly))',1,Lx);
  end

  HAM.info.XY=[ reshape(XY{1}, [],1), reshape(XY{2}',[],1) ];

  HH=zeros(2*(Lx-1)*Ly,5); l=1;

  ee=expand_param(-mu, Ly,Lx);
  uu=expand_param(U,   Ly,Lx);
  tx=expand_param(tx,  Ly,Lx);
  ty=expand_param(ty,  Ly,Lx);

  if perBCx && Lx<4, wberr('invalid Lx=%g !? (having perBC-x)',Lx); end

  for j=1:Lx
  for i=1:Ly, k0=i+Ly*(j-1);

     HH(l,:)=[ [k0, 1], [k0, 1], ee(i,j)]; l=l+1;
     HH(l,:)=[ [k0, 3], [k0, 3], uu(i,j)]; l=l+1;

     if i<Ly,       HH(l,:)=[ [k0,      2], [k0+1, 2], ty(i,j)]; l=l+1;
     elseif perBCy, HH(l,:)=[ [k0+1-Ly, 2], [k0,   2], ty(i,j)]; l=l+1;
     end

     if ~perBCx
        if j<Lx
           HH(l,:)=[ [k0, 2], [k0+Ly, 2], tx(i,j)]; l=l+1;
        end
     else
        if j>1,  kl=k0-Ly; else kl=k0; end
        if j<Lx, kr=k0+Ly; else kr=k0; end
        HH(l,:)=[ [kl, 2], [kr, 2], tx(i,j)]; l=l+1;
     end
  end
  end

  HH=HH(find(HH(:,end)),:);

  HAM=setup_mpo(HAM,HH);

  HAM.store='DMRG_Hubbard';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %
function q=expand_param(q,Ly,Lx);

  if numel(q)==1, q=repmat(q,Ly,Lx);
  elseif ~isequal(size(q),[Ly,Lx]), wberr(...
    'invalid usage (size mismatch in %s: %s having (Ly,Lx)=(%g,%g)',...
     inputname(1), vsprintf(size(q),'x'), Ly, Lx);
  end

end

% -------------------------------------------------------------------- %

