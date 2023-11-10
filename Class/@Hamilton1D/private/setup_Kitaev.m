function [HAM]=setup_Kitaev(varargin)
% function [HAM]=setup_Kitaev([opts])
%
%    Setup of 2D Kitaev model on 2D honeycomp lattice.
%    as zig-zag (ZZ) tube (like vertical bricklayer)
%
% Options
%
%    'L',...   length of system
%    'W',..    width of system (4)
%    'NC',..   number of channels (1: default, or 2 with J2)
%    'J2',..   include next-nearest neighbor Heisenberg interaction
%
%    '--YL'    Yao-Lee spin liquid Hamiltonian
%              include Heisenberg spin operator with NN Kitaev terms
%    '--AC'    use arm-chair (AC) geometry instead of ZZ.
%
% Wb,Oct23,21

% wsys='Kitaev'; L=6; W=6; J2=0.5; tflag=2; tst_Hamilton1D

  global param

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     L  = getopt('L',  []);
     W  = getopt('W',  4);
     pBC= getopt('--pBC');

     zigzag=1;
     if getopt('--AC'), zigzag=0; end

     YLHam=getopt('--YL');
     use_mpo = getopt('--mpo',{'nsw',3});

     NC = getopt('NC',1);
     J2 = getopt('J2',0);
     K3 = getopt('K3',-1.);

     tflag= getopt('-t');

  if isempty(L)
       L=getopt('get_last',12);
  else getopt('check_error'); end

  if zigzag
       wstr='zigzag';
  else wstr='armchair'; end

  if mod(W,2), wbdie(...
    'invalid usage (%s requires even width; got W=%g)',wstr,W); end
  N=L*W; 

  if zigzag, s='ZZ'; else s='AC'; end
  param.system=sprintf('%s%gx%g',s,W,L);

  gotJ2=(~isempty(J2) && (ischar(J2) || norm(J2)));
  s='Kitaev';
     if NC>1 && (YLHam && gotJ2 || ~YLHam || NC~=2), s=[s num2str(NC)]; end
     if YLHam, s=[s 'YL']; end
  param.istr=[s,' ',param.system];

  wblog('<i>','setting up %s',param.istr); 

  sym={'Spin',1/2,'--Z2'};
  [T3,IS]=getLocalSpace(sym{:},'-v');

  if norm(T3(1).Q{3})
     wbdie('unexpected spin operator order'); end; T3_=T3;

  if NC==1, E=IS.E;
     if YLHam, wbdie('invalid usage (--YL requires NC=2)'); end
  elseif NC==2
     T3=addSymmetry(T3,'SU2');

     [Sop,IS(2)]=getLocalSpace(sym{1:end-1},'-v');
     Sop=addSymmetry(Sop,'Z2','pos',1);

     A=getIdentity(sum(T3),Sop,[1 3 2]);
     E=getIdentity(A,2);

     for i=1:numel(T3)
        T3(i)=contract(A,'13*',contractQS(T3(i),2,A,1),'14',[1 3 2]);
     end
     Sop=contract(A,'13*',contractQS(A,3,Sop,2),'13');

     if YLHam
        A=getIdentity(sum(T3),3,Sop,3);
        for i=1:numel(T3)
           T3(i)=contract(contract(T3(i),2,Sop,1),'24',A,'12');
        end
     end

  else wbdie('invalid usage (NC=%g)',NC); 
  end

  sym=E.info.qtype;

  add2struct(param,sym,L,W,YLHam,pBC,zigzag,NC,K3,J2);

  HAM=struct(Hamilton1D);

  HAM.info.istr=param.istr;
  HAM.info.IS=IS;

  HAM.info.lops=T3;
  if NC>1, HAM.info.lops(end+1)=Sop; end

  HAM.oez=init_ops(E,'local identity operator (E)');

  T3 = [ (T3(2)-T3(3))/sqrt(2)
         (T3(2)+T3(3))/sqrt(2)
          T3(1)
  ];

  if numel(K3)==1, K3=repmat(K3,1,3); end

  Q=QSpace; for i=1:numel(T3), Q=Q+contractQS(T3(i),'13*',T3(i),'13'); end
     q=3/4; if YLHam, q=q*q; end
     e=normQS(Q-q*E);
  if e>1E-12, wbdie('unexpected spin operator (e=%.3g)',e); end

  if YLHam,
       s={'YaoLee ','.S'};
  else s={'Kitaev ',''  }; end

  HAM.ops     =init_ops(T3(1),[s{1} 'Tx' s{2}],'~hconj');
  HAM.ops(2,1)=init_ops(T3(2),[s{1} 'Ty' s{2}],'~hconj');
  HAM.ops(3,1)=init_ops(T3(3),[s{1} 'Tz' s{2}],'~hconj');

  if NC>1
  HAM.ops(4,1)=init_ops(Sop,'spin Sop','~hconj'); end

  Hloc=QSpace;

  q=1/6;
  x4=[ 1+q, 2-q
       1-q, 2+q ]; 

  if zigzag
  % DMRG site order for zigzag (ZZ) cylinder ------------------------- %
  % ZZ honeycomb cylinder is like brick layer in vertical direction
  % around the cylinder, with x,y,z bonds in T3 around every site
  % Setup for W=4; on the particular case of W=4:
  %
  %   - the vertical next-nearest neighbor J2 bonds appear
  %     twice for the same pair of sites, with one also
  %     wrapping around the cylinder via perBC-y
  %
  %   - the longest NNN bond stems from wrapping around the perBC-y
  %     boundary, e.g. (W+1) - 3W => range = 3W - (W+1) + 1 = 2W
  %     (having range 2 for shortest NN interaction)
  %
  %   - th case W=4 actually turns out gapped
  %     as it does not include the K points => for critical,
  %     width must be multiple of 3, and for consisteny
  %     with perBC-y, actually a multiple of 6 (!)
  %  => smallest width of interest for ZZ is 6 (!)
  %
  %  (1)----(W+1)   (2W+1)-   ...((L-1)*W+1)   ZZ::
  %   :       :       :            :
  %   z       y       |            |
  %   |       |       |            |
  %   W      2W  -x- 3W       ... L*W
  %   |       |       |            |
  %   y       z       |            |
  %   |       |       |            |
  %   3 -x-- W+3     2W+3 -x- ... (L-1)*W+1
  %   |       |       |            |
  %   z       y       |            |
  %   |       |       |            |
  %   2      W+2 -x- 2W+2     ... (L-1)*W+1
  %   |       |       |            |
  %   y       z       |            |
  %   |       |       |            |
  %   1 -x-- W+1     2W+1 -x- ... (L-1)*W+1
  %
  % ------------------------------------------------------------------ %

    q=0:2:L;
    q=[ sort([x4(1,1)+q, x4(1,2)+q])
        sort([x4(2,1)+q, x4(2,2)+q]) ];
    q={ repmat(q(:,1:L),W/2,1), repmat(1:W,L,1)' };

  else
  % DMRG site order for armchair (AC) cylinder ----------------------- %
  % AC honeycomb cylinder is like brick layer in horizontal direction
  % along the cylinder, with x,y,z bonds in T3 around every site;
  %
  %  - advantage of AC over ZZ cylinder:
  %    -> *always* critical (e.g. metallic for armchair carbon nanotubes)
  %       i.e. K points are always included
  %    => can use W=4 here, even W=2 (ladder!)
 %
  %  - longest range bond is horizontal NNN interaction has r=2W+1
  %    i.e., is comparable to ZZ; yet there are significantly more
  %    of these long-range bonds, because they are intrinsic 
  %    do not relate to perBC-y, so the MPO will have signifianctly
  %    larger bond dimension (like 2W+2 for zigzag => 4W for armchair)
  %
  %  - the circumference of the AC cylinder is C_ac = 3W/2 [=(W/2)*1 + (W/2)*2]
  %    which for a given W is smaller as compared to the ZZ cylinder
  %    which has C_zz = sqrt(3)*W = 1.732*W. In particular,
  %    AC4 (i.e. W=4) has C=6, as compared to ZZ6 with C=6*sqrt(3) = 10.39
  %    The MPO bond dimension of AC4 (Dmpo=16) is comparable, however,
  %    to ZZ6 (Dmpo=14)
  %
  % Site setup (shown here for W=4): // AC::
  %
  %  (1) x- (W+1) y (2W+1) x (3W+1) y (4W+1) -x-    ...((L-1)*W+1)
  %           :                :                 :
  %           z                z                 z
  %           |                |                 |
  %   W -y-- 2W --x-- 3W --y-- 4W --x-- 5W --y--    ... L*W
  %   |                |                |                |
  %   z                z                z                z
  %   |                |                |                |
  %   3 -x-- W+3 -y- 2W+3 -x- 3W+3 -y- 4W+3 -x-     ... (L-1)*W+1
  %           |                |                 |
  %           z                z                 z
  %           |                |                 |
  %   2 -y-- W+2 -x- 2W+2 -y- 3W+2 -x- 4W+2 -y-     ... (L-1)*W+1
  %   |                |                |                |
  %   z                z                z                z
  %   |                |                |                |
  %   1 -x-- W+1 -y- 2W+1 -x- 3W+1 -y- 4W+1 -x-     ... (L-1)*W+1
  %
  % ------------------------------------------------------------------ %

    q=0:2:W;
    q=[ sort([x4(1,1)+q, x4(1,2)+q])
        sort([x4(2,1)+q, x4(2,2)+q]) ];
    q={ repmat(1:L,W,1), repmat(q(:,1:W),ceil(L/2),1)' };
    if mod(L,2), q{2}=q{2}(1:N); end
  end

  HAM.info.XY=[ q{1}(:), q{2}(:) ];

  SS=cell(1,3*L); l=1;

  if pBC, wblog('NB!','using full periodic BC (pBC)'); 
       oi={'-pX'}; 
  else oi={}; end

  oJ2={};
  if gotJ2
     if ischar(J2), % e.g. '--lin:0:2'
        q=textscan(J2,'%s','whitespace',':, '); q=q{1}.';
        if isequal(q{1},'--lin'), oJ2={q{1}, [str2num(strjoin(q(2:end)))]};
        else J2, wbdie('invalid usage'); end
     else
        n=numel(J2);
        if n==2, oJ2={'--lin',J2};
        elseif n==N, oJ2={'--map',J2};
        elseif n==1, oJ2={'--const',J2};
        else J2, wbdie('invalid usage (L=%g)',L); end
     end
     if norm(oJ2{2})==0, oJ2={}; end
  end

  if ~isempty(oJ2)
     if NC<2, wbdie(...
       'invalid usage (got J2 data but having NC=%g)',NC); end
     get_J2_coupling('--init',HAM.info.XY);
  end
  param.oJ2=oJ2;
  HAM.info.param=param;

  for i=1:L
     I=repmat(i,W/2,1); J=(1:2:W)';
     if mod(i,2), w=1; else w=2; end

     for p=1:3
        if zigzag
           if i==1 && ~pBC && p>2, break; end
           switch p
             case 1, I2=[I, I];  J2=[J,J+1]; j=w+1;
             case 2,             J2=J2+1;    j=4-w;
             case 3, I2=[I-1,I]; J2=[J, J]+(2-w); j=1;
             otherwise wbdie('invalid switch'); 
           end
        else
           if i==1 && ~pBC && p>1, break; end
           switch p
              case 1, I2=[I, I];   J2=[J-1, J]+w; j=3;
              case 2, I2=[I-1, I]; J2=[J, J]; j=w;
              case 3,              J2=J2+1; j=3-w;
              otherwise wbdie('invalid switch'); 
           end
        end

        kk=sub2ind_cylinder(I2,J2,L,W,oi{:});  nk=size(kk,1);
        S=repmat(struct('iop',[0 j; 0 j],'hcpl',K3(j)),1,nk);
        for j=1:nk, S(j).iop(:,1)=kk(j,:); end
        SS{l}=S; l=l+1;
     end

     if isempty(oJ2), continue; end
     I=repmat(i,W,1); J=(1:W)';

     for p=1:3
        if zigzag
           if i==1 && ~pBC && p>1, break; end
           switch p
             case 1, I2=[I  , I]; J2=[J,J+2];
             case 2, I2=[I-1, I]; J2=[J,J-1];
             case 3, I2=[I-1, I]; J2=[J,J+1];
             otherwise wbdie('invalid switch'); 
           end
        else
           if i<=2 && ~pBC
              if i==1 || i==2 && p>2, break; end
           end
           switch p
             case 1, I2=[I-1, I]; J2=[J,J-1];
             case 2, I2=[I-1, I]; J2=[J,J+1];
             case 3, I2=[I-2, I]; J2=[J,J  ];
             otherwise wbdie('invalid switch'); 
           end
        end

        kk=sub2ind_cylinder(I2,J2,L,W,oi{:});  nk=size(kk,1);
        J2cpl=get_J2_coupling(kk,oJ2{:});

        S=repmat(struct('iop',[0 4; 0 4],'hcpl',[]),1,nk);
        for j=1:nk
           S(j).iop(:,1)=kk(j,:);
           S(j).hcpl=J2cpl(j);
        end
        SS{l}=S; l=l+1;
     end
  end

  if L*W<=4 || isempty(use_mpo)
     HH=convert_HS([SS{:}]);
     HAM=setup_mpo(HAM,HH);
  else
     HAM.info.HS=[SS{:}];
     HAM=setup_mpo_full(HAM,use_mpo{:}); % 'nsw',3
  end

  s='DMRG_Kitaev'; if NC>1, s=[s num2str(NC)]; end
  HAM.store=s;

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

function J2=get_J2_coupling(varargin)

  persistent XY x1 x2

  if isequal(varargin{1},'--init')
     XY=varargin{2}; x1=min(XY(:,1)); x2=max(XY(:,1));
     return
  end

  n=numel(varargin);
  if n~=3, wbdie('invalid usage (%g/3 args)',n); 
  elseif isempty(XY), wbdie('invalid usage (not yet initialized)'); end

  kk=varargin{1}; wJ2=varargin{2}; J2=varargin{3}; n=numel(J2);
  switch wJ2
     case '--lin'
         if n~=2, wbdie('invalid usage (%s: %g/2)',wJ2,n); end
         x=mean(reshape(XY(kk(:),1)-x1,size(kk)),2)/(x2-x1);
         J2=J2(1)+x*diff(J2);
     case '--map'
         if n~=2, wbdie('invalid usage (%s: %g/2)',wJ2,n); end
         J2=mean(reshape(J2(kk(:)),size(kk)),2); e=numel(J2)<=2;
     case '--const'
         if n~=1, wbdie('invalid usage (%s: %g/2)',wJ2,n); end
         J2=repmat(J2,size(kk,1),1);
     otherwise wJ2, wbdie('invalid switch'); 
  end

end

% -------------------------------------------------------------------- %

