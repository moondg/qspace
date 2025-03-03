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
%    'tx',..   hopping amplitude in x-direction (t)
%    'ty',..   hopping amplitude in y-direction (t)
%
%    'NC',..   number of channels (spinless, unless '--spin'; default: 1)
%    '--spin'  use spinfull fermions (by default: spinless fermions)
%    '--perBC' periodic boundary condition
%              (if Ly<=2, in x-, otherwise y-direction)
%    '--dcyl'  wrap cylinder along diagonal
%    '--ph'    choose particle/hole symmetric sector
%              (only valid if 'mu' is not specified)
%
% Wb,Feb06,18

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  wcyl='square';

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

     t =getopt('t', []); if isempty(t), t=-max([1,sqrt(abs(U))]); end
     tx=getopt('tx',[]); if isempty(tx), tx=t; end
     ty=getopt('ty',[]); if isempty(ty), ty=t; end

     perBC= getopt('--perBC');
     if     getopt('--spin' ), Sflag=+1;
     elseif getopt('--spinB'), Sflag=-1; else Sflag=0; end
     u1flag=getopt('--U1charge');

     if getopt('--dcyl')
        if Ly>1, wcyl='diag';
        else wblog('WRN','got --dcyl with Ly=%d (ignore)',Ly); end
     end
     if Ly<=1, wcyl='chain'; end

     tflag= getopt('-t');
     vflag=~getopt('-q');

  if isempty(Lx)
     Lx=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(Lx), wbdie('length L not specified'); end

  sym={'Fermion',''}; if ~u1flag && any(mu(:)), u1flag=1; end

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

  perBCx=0;
  perBCy=0; if Ly<1, wbdie('invalid Ly=%d',Ly); end

  if perBC
     if Ly<2
        wblog(' * ','using interleaved setup for perBC-x (Ly=%g)',Ly);
        perBCx=1;
     else
        wblog(' * ','using perBC-y (Ly=%g)',Ly);
        perBCy=1;
     end
  end

  param=add2struct('-',Lx,Ly,t,tx,ty,U,mu,NC,sym,wcyl,perBC);
  param.PERBC=[perBCx, perBCy];

  if Sflag
       [F,Z,Sop,IS]=getLocalSpace(sym{:},'NC',NC,'-v');
  else [F,Z,    IS]=getLocalSpace(sym{:},'NC',NC,'-v');
  end
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
     init_ops(Z,   'fermionic parity (Z)')
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

  HH=zeros(2*(Lx-1)*Ly,5); l=1;

  tx=expand_param( tx, Ly,Lx);
  ty=expand_param( ty, Ly,Lx);

  if numel(mu)==1 && numel(U)==1 && U && mu
     n=numel(HAM.ops); if n~=3, wbdie('got %d HAM.ops !?',n); end
     Hloc=(-mu)*HAM.ops(1).op + U*HAM.ops(3).op;
     HAM.ops(1)=init_ops(Hloc,sprintf('Hloc(epsd=%.4g, U=%.4g)',-mu,U),'~hconj');
     HAM.ops(end)=[];

     for j=1:Lx
     for i=1:Ly, k=i+Ly*(j-1);
        HH(l,:)=[ [k, 1], [k, 1], 1.0]; l=l+1;
     end
     end
  else
     ee=expand_param(-mu, Ly,Lx);
     uu=expand_param( U,  Ly,Lx);

     for j=1:Lx
     for i=1:Ly, k=i+Ly*(j-1);
        if ee(i,j), HH(l,:)=[ [k, 1], [k, 1], ee(i,j)]; l=l+1; end
        if uu(i,j), HH(l,:)=[ [k, 3], [k, 3], uu(i,j)]; l=l+1; end
     end
     end
  end

if isequal(wcyl,'square') || isequal(wcyl,'chain')
% DMRG site order (square wrapping) -------------------------------- %
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

  if perBCx && Lx<4, wbdie('invalid Lx=%g !? (having perBC-x)',Lx); end

  for j=1:Lx
  for i=1:Ly, k=i+Ly*(j-1);

     if i<Ly,       HH(l,:)=[ [k,      2], [k+1, 2], ty(i,j)]; l=l+1;
     elseif perBCy, HH(l,:)=[ [k+1-Ly, 2], [k,   2], ty(i,j)]; l=l+1;
     end

     if ~perBCx
        if j<Lx
           HH(l,:)=[ [k, 2], [k+Ly, 2], tx(i,j)]; l=l+1;
        end
     else
        if j>1,  kl=k-Ly; else kl=k; end
        if j<Lx, kr=k+Ly; else kr=k; end
        HH(l,:)=[ [kl, 2], [kr, 2], tx(i,j)]; l=l+1;
     end
  end
  end

elseif isequal(wcyl,'diag')
% DMRG site order (diagonal wrapping) ------------------------------ %
%
%        i=1   2    3    4    5    6    7    ..      Lx
%    j=
%    1'       Ly+1    3Ly+1      5Ly+1       ..    (Lx-1)*Ly+1
%             /  \      /  \      /  \      /  \      /   
%           /      \  /      \  /      \  /      \  /    
%    1    1       2Ly+1     4Ly+1       *         ..   
%           \      /  \      /  \      /  \      /  \    
%             \  /      \  /      \  /      \  /      \  
%    2'      Ly+2     3Ly+2      5Ly+2       ..       ..
%             /  \      /  \      /  \      /  \      /   
%           /      \  /      \  /      \  /      \  /    
%    2    2       2Ly+2    4Ly+2       ..         .. 
%           \      /  \      /  \      /  \      /  \    
%    ..     ..   ..   ..   ..   ..   ..   ..   ..    ..
%             \  /      \  /      \  /      \  /      \  
%    Ly'      2Ly       4Ly       6Ly        *       Lx*Ly
%             /  \      /  \      /  \      /  \      /   
%           /      \  /      \  /      \  /      \  /    
%    Ly   Ly       3Ly       5Ly        ..        ..     
%           \      /  \      /  \      /  \      /  \    
%             \  /      \  /      \  /      \  /      \  
%             (*)        (*)       (*)      (*)       (*)
%
% ------------------------------------------------------------------ %

  if perBCx, s=['perBCx not yet implemented for ' wcyl '-cylinder'];
     wbdie('invalid usage (%s)',s); end
  if mod(Ly,2)
     wbdie('invalid Ly=%g (must be even for %s-cylinder)',Ly,wcyl);
  end

  HAM.info.xops=[HAM.ops(2).op, Sop];

  XY={ repmat(0:Lx-1,   Ly,1)
       repmat((Ly-1:-1:0)',1,Lx) };

  XY{2}(:,2:2:end)=XY{2}(:,2:2:end)+0.5;

  HAM.info.XY=[ reshape(XY{1}, [],1), reshape(XY{2},[],1) ];

  if perBCx && Lx<4, wbdie('invalid Lx=%g !? (having perBC-x)',Lx); end
  if ~perBCy, wbdie('strip BC not implemented with %s-cylinder',wcyl); end

  for j=1:Lx-1
     J=repmat(j,Ly,1); J2=[J J+1]; 
     I=(1:Ly)';

     if mod(j,2), di=[0 1]; else di=[-1 0]; end
     for w=1:2, I2=[I, I+di(w)];
        kk=sub2ind_cylinder(J2,I2,Lx,Ly);
        if w==1, t=tx(:,j); else t=ty(:,j); end
        for i=1:Ly
            HH(l,:)=[ [kk(i,1), 2], [kk(i,2), 2], t(i)]; l=l+1;
        end
     end
  end
else
  wcyl, wbdie('invalid switch'); 
end

  HH=HH(find(HH(:,end)),:);

  HAM=setup_mpo(HAM,HH,'-q');

  HAM.store='DMRG_Hubbard';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %
% allow inhomogenous setting of parameters to be specified as input
% otherwise expand

function q=expand_param(q,Ly,Lx);

  if numel(q)==1, q=repmat(q,Ly,Lx);
  elseif ~isequal(size(q),[Ly,Lx]), wbdie(...
    'invalid usage (size mismatch in %s: %s having (Ly,Lx)=(%g,%g)',...
     inputname(1), vsprintf(size(q),'x'), Ly, Lx);
  end

end

% -------------------------------------------------------------------- %

