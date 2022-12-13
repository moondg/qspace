function [HAM]=setup_KondoNecklace(varargin)
% function [HAM]=setup_KondoNecklace([L, opts])
% Wb,May24,18

% adapted from setup_HeisenbergLadder()

  if nargin<1
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  getopt('init',varargin);

     mJ2   = getopt('mJ2',[0 1 0]);
     L     = getopt('L', []);
     NC    = getopt('NC',1);
     tpar  = getopt('tpar',1);
     perBC = getopt('-perBC');
     tflag = getopt('-t');
     sym   = getopt('sym','SU2spin,Acharge');

  if isempty(L) L=getopt('get_last',[]); else
  getopt('check_error'); end

  mJ2(:,end+1:3)=0;
  [l,m]=size(mJ2); if m~=3, wberr(...
    'invalid couplings mJ2 (expecting [my, J_K, J_H]) !?'); end

  if isempty(L)
     wblog(' * ','got implicit length specification (L=%g)',l);
     L=l;
  elseif l==1, mJ2=repmat(mJ2,L,1);
  end

  if L<=1, wberr('invalid usage (system length not specified !?)'); end
  if L~=size(mJ2,1), wberr(...
    'invalid usage (L=%g; mJ2: [%s]) !?',L,vec2str(size(mJ2),'-f'));
  end

  o={'-v'};
  HAM=struct(Hamilton1D);

  param=add2struct('-',L,perBC);

  jstr={'tpar','mu','JK','JH'}; s=zeros(1,4);
  for l=1:4
     if l==1, q=-tpar; else q=unique(chopd(mJ2(:,l-1))); end
     s(l)=numel(q);
     if s(l)==1
        if l>1, q=mJ2(1,l-1); end
        param=setfield(param,jstr{l},q);
        jstr{l}=sprintf('%s=%.4g',jstr{l},q);
     else
        param=setfield(param,jstr{l},mJ2(:,l));
        jstr{l}=sprintf('%s=[%.3g %.3g]',jstr{l},q(1),q(end));
     end
  end

  jstr(end+1,:)={', ',' (',', ',')'};
  param.jstr=[jstr{:}];

  if norm(diff(mJ2,[],1))<1E-12
       param.mJ2=mJ2(1,:);
  else param.mJ2=mJ2; end

  HAM.info.mpo.jstr=param.jstr;
  HAM.info.mpo.mJ2=param.mJ2;

  [HK,Fb,Zb,Nb,S2,E1,IS]=get_ops_KondoLattice(NC,sym);

  HAM=struct(Hamilton1D);
  HAM.info.istr=sprintf('Kondo necklace having %s',IS.istr);
  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.oez=[
     init_ops(E1,'local identity operator (E2)')
     init_ops(Zb,'fermionic parity (Z)')
  ];

  HAM.ops=init_ops(Fb,'fermionic hopping','-ferm');
  p2=repmat(-tpar,size(mJ2,1),1);

  Hflag=any(mJ2(:,3));
  if Hflag
     HAM.ops(end+1,1)=init_ops(S2(1),'Heisenberg','~hconj');
     p2(:,end+1)=mJ2(:,3);
  end

  p1={};

  mJ2(:,1)=-mJ2(:,1);

  q=uniquerows(chopd(mJ2(:,1:2)));
  uloc=(size(q,1)==1);
  if uloc
     if any(mJ2(1,1:2)), q=mJ2(1,1)*Nb+mJ2(1,2)*HK;
        HAM.ops(end+1,1)=init_ops(q,'Hlocal','~hconj');
        p1=ones(size(mJ2,1),1);
     end
  else
     if any(mJ2(:,1))
        HAM.ops(end+1,1)=init_ops(Nb,'chemical potential','~hconj');
        p1{end+1}=mJ2(:,1);
     end
     if any(mJ2(:,2))
        HAM.ops(end+1,1)=init_ops(HK,'(S.S)_Kondo','~hconj');
        p1{end+1}=mJ2(:,2);
     end
     p1=[p1{:}];
  end

  wblog('==>','%s',param.jstr);

  HAM.info.xops=[HK; Nb; Fb; S2(:)];
  HAM.info.xops_istr='HK,Nb,Fb,Sa,Sb where b=bath, a=K, 2=b';

  n1=size(p1,2); n2=size(p2,2);

  if perBC
     wblog('NB!','using periodic BC (interleaved setup)'); 
     L2=ceil(L/2);

     XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
     HAM.info.XY=XY; l=1;

     for k=1:L
        for j=1:n1
           HH(l,:)=[ [k,n2+j], [k,n2+j], p1(k,j)]; l=l+1;
        end
        if k+1<L, k2=k+2; else k2=k+1; end
        for j=1:n2
           if k ==1, HH(l,:)=[ [1,j], [ 2,j], p2(k,j)]; l=l+1; end
           if k2<=L, HH(l,:)=[ [k,j], [k2,j], p2(k,j)]; l=l+1; end
        end
     end
  else
     XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
     HAM.info.XY=XY; l=1;

     for k=1:L
        for j=1:n1
           HH(l,:)=[ [k,n2+j], [k,n2+j], p1(k,j)]; l=l+1;
        end
        if k<L
        for j=1:n2
           HH(l,:)=[ [k, j  ], [k+1, j], p2(k,j)]; l=l+1;
        end
        end
     end
  end

  HAM=setup_mpo(HAM,HH);

  HAM.store='DMRG-KNL';

  if tflag, plot(HAM), end

end

% -------------------------------------------------------------------- %
function [HK,Fb,Zb,Nb,S2,E1,IS]=get_ops_KondoLattice(NC,sym)

% NB! build super-site to simplify structure
%
%      1a   2a   3a   4a      // a are the Kondo spins
%      |    |    |    | 
%      |    |    |    | 
%      1b---2b---3b---4b--    // b like b(ath)
%
% site 1    2    3    4  ...

  [F,Z,Sb,IS]=getLocalSpace('FermionS',sym,'NC',NC,'-v');

  if numel(F )>1, F =cat(F, 3); end
  if numel(Sb)>1, Sb=cat(Sb,3); end

  Sa=Sb;
  if isempty(sym)
     if numel(Sa.data)~=1 || size(Sa(1).data{1},1)~=4
        wberr('expecting NC=1 (having nosym)'); end
     Sa.data{1}=Sa.data{1}(2:3,2:3,:);
  end

  if ~isfield(IS,'istr')
     IS.istr=sprintf('spinfull fermionic site (NC=%g: %s)',NC,sym);
  end

  A1=QSpace(permuteQS(getIdentityQS(Sa,1,IS.E,1),'132'));
  A1.info.itags={'a','ab*','b'};
  E1=untag(QSpace(getIdentityQS(A1,2)));

  Q=contractQS(F,'13*',F,'13');
  Nb=untag(QSpace(contractQS(A1,'!2*',{A1,Q,'-op:b'})));

  Zb=untag(QSpace(contractQS(A1,'!2*',{A1,Z,'-op:b'})));

  S2=QSpace(2,1);
  S2(1)=contractQS(A1,'!2*',{A1,Sa,'-op:a'});
  S2(2)=contractQS(A1,'!2*',{A1,Sb,'-op:b'});

  for i=1:numel(S2)
     S2(i).info.otype='operator';
     S2(i)=untag(S2(i));
  end

  Fb=F;
  for i=1:numel(F)
     Fb(i)=contractQS(A1,'!2*',{A1,F(i),'-op:b'});
     Fb(i).info.otype='operator';
     Fb(i)=untag(Fb(i));
  end

  HK=untag(QSpace(contractQS(A1,'!2*',{Sa,'-op:a', {A1,Sb,'*','-op:b'}})));

end

% -------------------------------------------------------------------- %

