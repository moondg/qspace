function [HAM]=setup_HeisenbergTriLadder(varargin)
% function [HAM]=setup_HeisenbergTriLadder([L, opts])
%
%    setup triangular ladder, also allowing for chiral 3-site
%    correlations for up/down triangles.
%
%    The triangular ladder is equivalent to a 1D chain
%    with NNN coupling J2, where with
%    J = [ Jleg, Jtri, chiral ] \equiv [ J2 J1, chiral ]
%    J2 comes first in the parameter listing!
%
%    This setup also permits wider system, where for Ly>1
%    the triangular ladders are coupled in Kagome like fashion
%    [e.g. cf. Kalmeyer-Laughlin spin liquid]
%
% Wb,Mar23,16

% adapted from setup_HeisenbergLadder()

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     sym  = getopt('sym','SU2');
     qloc = getopt('qloc', []);
     J    = getopt('J', [1 1 sqrt(2)]);
     Dz   = getopt('Dz',[]);
     L    = getopt('L', []);
     Ly   = getopt('Ly', 1);
     pBCy = getopt('-pBCy');
     perBC= getopt('-perBC');
     tflag= getopt('-t');

     Aflag= getopt('-A');
     DzI  = getopt('DzI',[]);

  if isempty(L) L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), L=size(J,1);
     if L>1, wblog(' * ','got implicit length specification (L=%g)',L);
     else wbdie('invalid usage (system length not specified)'); end
  elseif size(J,1)==1, J=repmat(J,L,1); end

  if ~isempty(Dz)
     if numel(Dz)>2, Dz, wbdie('invalid usage');
     elseif all(Dz==1), Dz=[];
     elseif Dz(1)==Dz(end), Dz=Dz(1);
     end
  end
  if ~isempty(DzI)
     if numel(DzI)~=1, DzI, wbdie('invalid usage');
     elseif DzI==1 || (Ly<2 && ~pBCy), DzI=0;end
  end

  if ~perBC, J(end,1)=J(end-1,1); end

  if size(J,1)~=L && size(J,2)~=3 && size(J,2)~=4, wbdie(...
     'invalid usage (L=%g; J: [%s]) !?',L,vec2str(size(J),'-f'));
  end

  o={'-v'};
  if regexp(sym,'^SU\d+')
     N=str2num(sym(3:end)); if N>4, o={'-V'}; end
  end

  if ~Aflag
     if     ~isempty(Dz ), Aflag=2;
     elseif ~isempty(DzI), Aflag=4; end
  end

  param=add2struct('-',L,sym,J,Dz,DzI,perBC,qloc);
  q=uniquerows(J); if size(q,1)==1, param.J=q; end

  HAM=struct(Hamilton1D);

  q=uniquerows(J); J3=[];
  if size(q,1)==1
     jstr=sprintf('J=[%g, %g]',q);
     HAM.info.mpo.J=q;
  else
     jstr=sprintf('J=[%.4g %.4g]..[%g %g]',min(J,[],1),max(J,[],1));
     HAM.info.mpo.J=J;
  end
  HAM.info.mpo.Dz=Dz;

  if size(q,1)==1 && size(J,2)>=3
     J3=J(1,1:3);
  end

  [S1,S2,Sleg,Hp1,Sp2,Sch,Slad,E2,IS]=...
  get_ops_TriLadder(sym,qloc,J3,Dz,DzI,Aflag);

  if Aflag, param.sym=getsym(S1(1)); end

  dloc=getDimQS(E2);
  s=sprintf('(%s) using supersite @ d=%g',IS.istr,dloc(end,1));
  if ~isempty(Dz), q=sprintf(' %g',Dz);
     if numel(Dz)==1, q=q(2:end); else q=['[' q(2:end) ']']; end
     s=[s, 10 'NB! Dz=' q ' incorporated in S.S!' 10];
  end

  HAM=struct(Hamilton1D);
  HAM.info.istr=['Heisenberg tri-ladder ' s];
  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.info.lops=[S2];

  if ~isempty(Dz)
      for i=1:numel(Dz), sD{i}=sprintf(' @ Dz=%g',Dz(i)); end
  else sD={''}; end

  HAM.oez=init_ops(E2,'local identity operator (E2)');

  HAM.ops=init_ops(Hp1,['rung Heisenberg term (S.S)_p1' sD{end}],'~hconj');

  if isempty(J3), J3=J;
     if numel(Sleg)~=1, wbdie('invalid usage'); end
     HAM.ops=[ HAM.ops
        init_ops(Sleg,['spin operator for (S.S)_legs' sD{1}],'~hconj')
        init_ops(Sp2, ['spin operator for (S.S)_p2' sD{end}],'~hconj')
        init_ops(Sch, ['spin operator for S.(SxS)'],         '~hconj')
     ];
  else J3=J(1,:);
     if numel(Dz)>1 
        sD=sprintf(' %g',Dz); sD={[' @ Dz=[' sD(2:end) ']']};
     end
     HAM.ops=[ HAM.ops
        init_ops(Slad,['combined NN Slad ops' sD{1}],'~hconj','-q')
     ];
  end

  if Ly>1 || pBCy
     HAM.ops=[ HAM.ops
        init_ops(S2,'rung spin ops for (S.S) with impurities', '~hconj')
     ];

     HAM.ops(1,2)=...
        init_ops(S1,'impurity spin ops for (S.S) with ladder','~hconj');
     HAM.oez(1,2)=init_ops(getIdentityQS(S1),'local identity operator (E2)');
  end

  if perBC 
     if Ly==1
        [HH,HAM.info.XY]=setup_HH_perBC(J3,L); stype=[];
     else wbdie('invalid usage (perBC only for Ly=1)'); end
  else
     [HH,HAM.info.XY,stype]=setup_HH_openBC(J3,L,Ly,pBCy);
  end

  HAM=setup_mpo(HAM,HH,stype);

  HAM.store='DMRG_HBTriL';

  if tflag, plot(HAM), end

end

% -------------------------------------------------------------------- %
% NB! build super-site to simplify structure of chiral term
%  +---+    +---+    +---+    +---+  
%  | 1 | -- | 3 | -- | 5 | -- | 7 | --  ...
%  | | |  / | | |  / | | |  / | | |  /
%  | | | /  | | | /  | | | /  | | | / 
%  | 2 | -- | 4 | -- | 6 | -- | 8 | --  ...
%  +---+    +---+    +---+    +---+ 
%  site_1   site_2   site_3   site_4    ...
%
% NB! changed triangular coupling from (1-4) to (2-3) etc. in order
% to have *consecutive* triangular (=NN) path: 1-2-3-4 ... // TRIANGULAR
% as this permits a natural interpretation of a 1D chain with NNN coupling
% Wb,Sep28,19
%
% See also / adapted from MEX/tst_HLTriTsvelik.m
% Wb,Mar24,16
% -------------------------------------------------------------------- %

function [S1,S2,Sleg,Hp1,Sp2,Sch,Slad,E2,I1]=get_ops_TriLadder(sym,qloc,J3,Dz,DzI,Aflag)

% returned operators
%
%    S1 spin operator of single site (on a rung; single QSpace, with Dz encoded!)
%    S2 spin operator in super-site basis for upper and lower rung (2 QSpaces)
%
%    Sleg    split of leg = (S.S)_upper_leg + (S.S)_lower_leg
%    Hp1     (S.S)_local    <= Jperp, part 1
%    Sp2     (S.S)_local    <= Jperp, part 1
%    Slad    split of total Hcpl, incl. chiral term
%
% NB! the 3-particle chiral term breaks particle-hole symmetry
% => S.(SxS) leads to an *antisymmetric* `Hamiltonian' in Hilbert space
%    [since S.(SxS) always contains a product Sx*Sy*Sz with Sy complex!]
% => based on real CGTs, this requires an imaginary coupling strength!
% e.g. see Archive/tst_HLTriTsvelik_160323.m  // CHIRAL_COMPLEX

  if ~isempty(DzI) && ~Aflag, Aflag=1; end
  if ~isempty(Dz)  && ~Aflag, Aflag=2; end

  if ~Aflag
     [S3,I2]=getLocalSpace('Spin',sym,2,'-v');
  else
     [S3,I2]=getLocalSpace('Spin',2/2,'-A','-v');
     S3=sum(S3); S3.info.otype='operator';
  end

  Z1=QSpace(getIdentityQS(S3,2,'-0'));
  if Aflag
     for i=find(Z1.Q{1}==0)', Z1.data{i}=-Z1.data{i}; end
  end
  Y3=QSpace(contractQS(contractQS(S3,2,Z1,1),2,Z1,1));

  if ~Aflag
       [S1,I1]=getLocalSpace('Spin',sym,qloc,'-v');
  else [S1,I1]=getLocalSpace('Spin',qloc/2,'-A','-v'); sym='U1';
  end
  if ~isfield(I1,'istr')
     I1.istr=sprintf('%s Sloc=%s',sym,num2rat(qloc/2));
  end

  Sl=[]; Sr=[];
  if Aflag
     if norm(S1(1).Q{3}), S1, wbdie('unexpected location of Sz-op'); end
     if ~isempty(Dz)
        if numel(Dz)>1, s='_leg  [=Dz2] '; else s=''; end
        if Dz(1)~=1
           wblog('-->','applying Dz%s=%g',s,Dz(1));
           Sl=Dz(1)*S1(1) + sum(S1(2:end));
           Sl.info.otype='operator';
        end
        if numel(Dz)==1, Sr=Sl;
        elseif Dz(end)~=1
           wblog('-->','applying Dz_rung [=Dz1] =%g',Dz(end));
           Sr=Dz(end)*S1(1) + sum(S1(2:end));
           Sr.info.otype='operator';
        end
     end
     S1=sum(S1); S1.info.otype='operator';

  elseif numel(S1)~=1, wbdie('unexpected S1'); 
  end

  A2=QSpace(permuteQS(getIdentityQS(S1,S1),[1 3 2]));
  E2=QSpace(getIdentityQS(A2,2));

  HH=QSpace(1,4);             S2=getS2(A2,S1);
  if isempty(Sl), SL=S2; else SL=getS2(A2,Sl); end
  if isempty(Sr), SR=S2; else SR=getS2(A2,Sr); end

  Hp1=QSpace(contractQS(SR(1),'13*',S2(2),'13'));

  HH(1)=QSpace(contractQS(SL(1),3,S2(1),'3*'))+contractQS(SL(2),3,S2(2),'3*');
  HH(2)=QSpace(contractQS(SR(2),3,S2(1),'3*'));

  for i=1:2 % (11',2'2) => (12,1'2')
     HH(i)=permuteQS(skipzeros(HH(i)),[1 4 2 3]);
  end

  C1=contractQS(contractQS(contractQS(S2(1),2,S2(2),1),'24',Y3,'12'),3,S2(1),'3');
  C2=contractQS(S2(2),'3',contractQS(contractQS(S2(2),2,S2(1),1),'24',Y3,'12'),3);
  HH(3)=permuteQS(1i*skipzeros(QSpace(C1)+C2),[1 3 2 4]);

  if numel(J3)==3
     HH(4) = J3(1)*HH(1) + J3(2)*HH(2) + J3(3)*HH(3);
  elseif isempty(J3)
     HH(4) = sum(HH(1:3));
  else J3, wbdie('invalid J3'); end

  [Q2,X2]=splitH4_SdotS(HH,'-v');

     if isempty(Dz) || all(Dz>=0)
        Sleg=Q2(1);
        e=normQS(Q2(1)-X2(1))/normQS(X2(1));
        if e>1E-12, wbdie('unexpected (S.S)_legs (e=%.3g !?)',e); end
     else
        Sleg=[ Q2(1), X2(1) ];
     end

     Sp2  = [ Q2(2), X2(2) ];
     Sch  = [ Q2(3), X2(3) ];
     Slad = [ Q2(4), X2(4) ];

  if ~isempty(DzI)
     Z=QSpace(getIdentityQS(S1,3));
     if ~isequal(Z.Q{1},[-2;0;2]) || ~isequal(Z.data,{1;1;1})
        wbdie('unexpected spin operator');
     end
     Z.data{2}=DzI;

     S1=QSpace(contractQS(S1,3,Z,1));x
     S1.info.otype='operator';
  end

end

% -------------------------------------------------------------------- %

function S2=getS2(A2,S1)

   S2=QSpace(1,2);
   S2(1)=contractQS(A2,'13*',contractQS(A2,1,S1,2),'32'); % LRs -> RsL'
   S2(2)=contractQS(A2,'13*',contractQS(A2,3,S1,2),'13');

end

% -------------------------------------------------------------------- %
%
%    [1]=[1+m]=[1+n]=[..]    double line = triangluar ladder
%      \ /       \ /
%      [2]      [3+n]        n=m+ml;
%      / \       / \
%    [3]=[3+m]=[4+n]=[..]
%      \ /       \ /
%      [4]       [6]
%      / \       / \
%    [5]=[5+m]=[7+n]=[..]
%      \ /       \ /
%      [6]      [6+n]       <= only if pBCy
%      / \       / \
%    ... ...   ... ...
%
% Wb,Mar25,16
% -------------------------------------------------------------------- %

function [HH,XY,stype]=setup_HH_openBC(JJ,L,Ly,pBCy)

  i=1; l=1;

  if size(JJ,2)<4 && (Ly>1 || pBCy), wbdie( ... 
    'invalid usage (J requires 4 elements having Ly=%g)',Ly); end
  if L<2, L, wbdie('invalid usage (L too short)'); end

  if size(JJ,1)==1, uniJs=1; Jx=JJ; ic=[4 5];
  else uniJs=0; ic=[7 8];
     if size(JJ,1)~=L, wbdie('size inconsistency L vs. J'); end
  end

  if Ly>1 && mod(L,2)
     wblog('WRN','reducing ladder length L to even (%g->%g)',L,L-1);
     L=L-1;
  end

  m2=[2*Ly+(pBCy-1), Ly];

  XY=zeros(sum(m2)*ceil(L/2),2); stype=zeros(1,size(XY,1));
  HH=zeros(6*L*Ly,5);

  for ix=1:L, odd=mod(ix,2);
     if odd
          mr=m2(1); ml=m2(2);
     else mr=m2(2); ml=m2(1);
     end; m_=mr;

     for iy=1:Ly
        if iy<Ly || pBCy
             if odd, gotimp=2; else gotimp=1; end
        else gotimp=0; end

        if ~uniJs, Jx=JJ(ix,:); end

        XY(i,  :)=[ix-1,  2*(Ly-iy)   ]; stype(i  )=1; if gotimp>1
        XY(i+1,:)=[ix-.5, 2*(Ly-iy)-1 ]; stype(i+1)=2; end

        HH(l,:)=[ [i,1], [i,1], Jx(2)]; l=l+1; % triangular rungs: rung (J')

        if ix<L
           if uniJs
              HH(l,:)=[ [i, 2], [i+mr,3], 1.   ]; l=l+1;
           else
              HH(l,:)=[ [i, 2], [i+mr,2], Jx(1)]; l=l+1;
              HH(l,:)=[ [i, 3], [i+mr,4], Jx(2)]; l=l+1; % triangular rungs: diag (J')
              HH(l,:)=[ [i, 5], [i+mr,6], Jx(3)]; l=l+1;
           end
        end

        if gotimp
           if odd, im=i+1; else im=i+1-ml; end
           if odd, i2=i+2; else i2=i+1;    end; if iy==Ly, i2=i2-m_; end

           HH(l,:)=[ [i,ic(2)], [im, 1], Jx(4)]; l=l+1;
           HH(l,:)=[ [im,1], [i2,ic(1)], Jx(4)]; l=l+1;
        end

        if gotimp>1
             i=i+2; mr=mr-1; ml=ml+1;
        else i=i+1; mr=mr+1; ml=ml-1; end
     end
  end

  XY(i:end,:)=[]; stype(i:end)=[];
  HH(l:end,:)=[];

end

% -------------------------------------------------------------------- %
% interleaved perBC setup // Wb,Apr01,16
% -------------------------------------------------------------------- %

function [HH,XY]=setup_HH_perBC(JJ,L)

  [l,m]=size(JJ);
  if m<3 || m>4, wbdie('invalid usage (J requires 3 columns having Ly=1)');
  elseif m==4, wblog('WRN','ignoring J4 coupling (having Ly=1)'); end

  if nargin<2
       uniJs=0; L=size(JJ,1);
  else uniJs=1;
     if l==1, JJ=repmat(JJ,L,1);
     elseif l~=L, wbdie('invalid usage (L=%g/%g !?)',L,l); end
  end

  if L<2, L, wbdie('invalid usage (L too short)'); end

  wblog('NB!',['setting up triangular ladder (L=%g, Ly=1)\n' ...
    'using interleaved perBC'],L);

  XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
  HH=zeros(3*L,5); l=1;

  for k=1:L
     HH(l,:)=[ [k,1], [k,1], JJ(k,2)]; l=l+1; % triangular rungs: rung (J')

     if k>1
          k1=k-1; if k<L, k2=k+1; else k2=k; end
     else k1=k; k2=k+1;
     end

     if mod(k,2)==0
        q=k1; k1=k2; k2=q;
     end

     if uniJs
        HH(l,:)=[ [k1,2], [k2,3], 1.     ]; l=l+1;
     else
        HH(l,:)=[ [k1,2], [k2,2], JJ(k,1)]; l=l+1;
        HH(l,:)=[ [k1,3], [k2,4], JJ(k,2)]; l=l+1; % triangular rungs: diag (J')
        HH(l,:)=[ [k1,5], [k2,6], JJ(k,3)]; l=l+1;
     end
  end

  HH(l:end,:)=[];

end

% -------------------------------------------------------------------- %

