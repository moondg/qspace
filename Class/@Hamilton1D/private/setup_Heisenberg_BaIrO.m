function [HAM]=setup_Heisenberg_BaIrO(varargin)
% function [HAM]=setup_Heisenberg_BaIrO([L, opts])
%
%    setup 2-leg brick-ladder with frustration term
%    see Ba4Ir3O10_structure_model.pdf
%
%    NB! the setup can be mapped to 3-site periodic 2-leg system
%    where a rung is made up of 3 sites.
%
%    Using isotropic Heisenberg couplings:
%      J(1)   coupling along legs (=J_1)
%      J(2)   direct rung coupling of legs (=J_2)
%      J(3)   coupling to center site (=J_3)
%      J(4)   incase case coupling to center site is asymmetric (=J_3')
%      J(5)   coupling of center sites (introduces 3rd leg; default 0)
%   
% Wb,Jun15,20

% [12/07/22]
%  - reversed order of (J2,J3) to align with paper notation
%  - removed option dJ2 and alpha=(J3-J2)/J1
%    merged into new column 4 in J(:,[3 4])
%    see new dJ3a in tst_Hamilton1D.m
% adapted from setup_HeisenbergTriLadder.m

  if nargin<1
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  getopt('init',varargin);
     sym  = getopt('sym','SU2');
     qloc = getopt('qloc', 1);
     J    = getopt('J', [1, 0.2, 0.2]);
     L    = getopt('L',[]);
     Dz   = getopt('Dz', 1); % motivated by Weiguo's Ising analysis
     perBC= getopt('--perBC');

     if getopt('-X'); xflag=1;
     elseif getopt('--Xzz'); xflag=2;
     else xflag=0; end

     tflag= getopt('-t');

     Aflag= getopt('-A');

  if isempty(L), L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), wberr('invalid usage (L not specified)'); end
  if L<4, wberr('invalid usage (L=%g)',L); end

  [nJ,m]=size(J); j=5; J_=J;
  if m<3 || m>j, wberr('invalid J data (%g x %g)',nJ,m); end

  if nJ==L-1 && ~perBC, l=L; else l=nJ; end
  if l>L || mod(L,l), wblog('WRN got length(J) = %d / %d',nJ,L); end

  q=min(abs(J),1); q(end+1:4)=0; q(5)=norm(q([3:4])); s='';
  if q(1)<1E-8, s=sprintf('decoupled rung(s) @ %.3g', q(1)); end
  if q(3)<1E-8, s=sprintf('invalid small J3=%.3g',q(3)); end
  if q(5)<1E-8, s=sprintf('decoupled center spin(s) @ %.3g',q(5)); end
  if ~isempty(s), wberr('invalid usage (got %s)',s); end

  if m<j, J(:,end+1:j)=0;
  if m<4, J(:,4)=J(:,3); end, end % J3=J3'

  Jc=norm(J(:,j));
  if Jc, Jc=Jc/sqrt(nJ); dJc=std(J(:,j));
     if dJc
          s=sprintf('Jc~%.3g @ %.3g',Jc,dJc);
     else s=sprintf('Jc=%.4g',Jc); end
     wblog('NB!','using 3-leg ladder with %s',s);
  end

  if ~Aflag && Dz~=1, Aflag=2; end

  istr=mfilename;
  param=add2struct('-',istr,L,J,Dz,sym,qloc,Aflag,xflag,perBC);

  [Sleg,Hloc,Eloc,Sloc,I1,param] = get_ops_BrickLadder(sym,param);
  param.J=J_;

  if Aflag, param.sym=getsym(Hloc); end

  dloc=getDimQS(Eloc); dloc=dloc(:,1);
  sx=I1.istr; if ~isempty(sx), sx=[sx ' with ']; end
  sx=['using ' sx 'supersite @ '];
  if numel(dloc)>1
       sx=[sx 'd^\ast=' sprintf('%g (%g)',dloc)];
  else sx=[sx 'd='      sprintf('%g',     dloc)];
  end

  jstr=cell(1,m); ndiff=zeros(size(jstr));
  for j=1:m, x=J(:,j);
     if j<=m
          q=[x(1), std(x)];
     else q=[0 0]; end

     if ~q(2), jstr{j}=sprintf('%.4g',q(1));
     else
        ndiff(j)=1;
        q(3)=mean(x); dx=x-q(3);
        if std(abs(dx))<1E-12
           if dx(1)>0, s='±'; else s='∓'; end
           jstr{j}=sprintf(['%.3g' s '%.3g'],q(3),abs(dx(1))/2);
        else
           q(3)=norm(x)/sqrt(nJ);
           if ~norm(diff(sign(x))), s=sign(q(1)); if ~s, s=1; end
                jstr{j}=sprintf( '%.3g @ %.3g',s*q(3),q(2));
           else jstr{j}=sprintf('~%.3g @ %.3g',q([3 2])); end
        end
     end
  end

  if any(ndiff), s=', '; else s=','; end
  jstr(2,:)={s};
  if m==4 && norm(diff(J(:,[3 4]),[],2))<1E-12, jstr(:,end)=[];
  elseif m>=4, jstr(2,[2:2:end])={'; '}; end

  HAM=struct(Hamilton1D);
  HAM.info.istr=['Heisenberg BaIrO ladder [' jstr{1:end-1} '] ' sx];
  HAM.info.param=param;
  HAM.info.IS=I1;

  HAM.info.lops=Sloc;

  if Aflag
     q=split_ops(Sloc);
     HAM.info.xops=fixScalarOp(q(:,2));
  else
     HAM.info.xops=Sloc;
  end

  n=numel(Sloc);
  for i=1:n, for j=i+1:n
     HAM.info.xops(end+1)=contract(Sloc(i),'13*',Sloc(j),'13');
  end, end

  HAM.oez=init_ops(Eloc,'local identity operator (E3)');

  HAM.ops=init_ops(Sleg,'(S.S)_legs','~hconj');
  for i=1:nJ
     if nJ==1, s='Hloc'; else s=sprintf('Hloc(%d)',i); end
     HAM.ops(i+1,1)=init_ops(Hloc(i),s,'~hconj');
  end

  stype=[]; l=1; HH=zeros(2*L-1,5);

  if ~perBC
     for i=1:L, i_=mod(i-1,nJ)+1;
        HH(l,:)=[ [i, 1+i_],  [i,   1+i_ ], J(i_,2) ]; l=l+1;  if i<L
        HH(l,:)=[ [i, 1   ],  [i+1, 1    ], J(i_,1) ]; l=l+1;  end
     end
  else
     wblog('NB!','using periodic BC (interleaved setup)'); 

     XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
     HAM.info.XY=XY;

     for i=1:L, i_=mod(i-1,nJ)+1;
        HH(l,:)=[ [i,   1+i_], [i,   1+i_ ], J(i_,2) ]; l=l+1;
        HH(l,:)=[ [i-1, 2   ], [i+1, 2    ], J(i_,1) ]; l=l+1;
     end

     HH(  2,1)=1;
     HH(end,3)=L;
  end

  HAM=setup_mpo(HAM,HH,stype);
  HAM.store='DMRG_BaIrO';

  if tflag, plot(HAM), end

end

% -------------------------------------------------------------------- %
% NB! build super-site to permit NN real-time evolution
% NB! the 6-site unit cell of the 2-leg ladder ...
%                                   
%     ┌───────────┐───────────┐───────
%     │ (1)---(6)-│-(7)---(C)-│-(D)-..
%     │  |        │  |        │  |
%     │ (2)       │ (8)       │ (E)
%     │  |        │  |        │  |
%     │ (3)---(4)-│-(9)---(A)-│-(F)-..
%     │        |  │        |  │
%     │       (5) │       (B) │
%     │        |  │        |  │
%     │ [1]...[6].│.[7]...[C].│.[D]...
%     └───────────┘───────────┘───────
%      SITE_1      SITE_2      SITE_3  
%
% can be equivalently redrawn as 2-leg ladder with 3-site periodicity (!) ...
%
%     ┌─────┐─────┐─────┐─────┐──────
%   a │ (1)-│-(6)-│-(7)-│-(C)-│-(D)..
%     │  |  │  |  │  |  │  |  │  |
%   b │ (2) │ (5) │ (8) │ (B) │ (E)        tags: b_center
%     │  |  │  |  │  |  │  |  │  |
%   c │ (3)-│-(4)-│-(9)-│-(A)-│-(F)..
%     └─────┘─────┘─────┘─────┘──────
%       s1    s2    s3    s4    s5 ..
%
% Wb,Jun15,20
% -------------------------------------------------------------------- %

function [Sleg,Hloc,Eloc,Sloc,I1,p]=get_ops_BrickLadder(sym,p)

  if size(p.J,2)~=5, wberr('invalid usage (J %dx%d)',size(p.J)); end
  istr={}; j3=[]; jc=0;

  if p.Aflag, p.sym='U1';
     [S1,I1]=getLocalSpace('Spin',p.qloc/2,'-A','-v');
     if p.Dz~=1, i=1;
        if p.Dz<0, wberr('invalid Dz=%g',p.Dz); end
        if numel(S1)~=3 || numel(S1(i).Q)~=3 || norm(S1(i).Q{3})
           wberr('unexpected S1 operators');
        end
        S1(i)=sqrt(p.Dz)*S1(i);
        istr{end+1}=sprintf('Dz=%g',p.Dz);
     end
     S1=sum(S1);
  else
     [S1,I1]=getLocalSpace('Spin',p.sym,p.qloc,'-v');
  end
  E1=I1.E;

  if size(p.J,2)>=5
     jc=p.J(:,5)./p.J(:,1); q=[norm(jc), std(jc)];
     if all(q), wberr(...
       'invalid usage (got varying Jc couplings @ %.3g)',q(2));
     end
     jc=jc(1);
     if jc, istr{end+1}=['jc=',num2rat(jc)]; end
  end

  [nJ,m]=size(p.J);

  if any(abs(p.J(:,2))<1E-3)
     wbdie('invalid usage (got small J2=%g)',min(abs(p.J(:,2)))); end

  j3=p.J(:,2:4)./repmat(p.J(:,2),1,3);

  q=mean(p.J(:,3:4));
  if nJ==2 && ~norm(diff(p.J(:,1:2),[],1)) && ~norm(diff(q,[],2))
     q=[ q(1), diff(p.J(1,3:4),[],1)/2 ];
     istr{end+1}=sprintf('jr=[%.3g, %.3g, %.3g±%.3g, %.3g∓%.3g]', ...
     [ p.J(1,1:2), q, q(1),-q(2) ]);
  elseif nJ<p.L-1
     if norm(diff(p.J(:,3:4),[],2)), j=1:3; else j=1:2; end
     istr{end+1}=vec2str(j3(:,j),'jr=[%.3g,]');
  else
     istr{end+1}='non-uniform';
  end
  I1.istr=strjoin(istr,', ');         % `jr' = normalized Js for rung

  A2=getIdentityQS(E1,E1);
  A3=getIdentityQS(A2,3,E1,1);
  A3=contractQS(A2,3,A3,1);

  X=contractQS(A3,'321*',A3,'123');
  [pp,IX]=eigQS(X); 
  Ux=IX.AK; Ux.data{1}=fliplr(Ux.data{1});
  A3=setitags(QSpace(contractQS(A3,4,Ux,1)),{'a','b','c','abc'});

  I1.Ux=Ux;
  I1.A3=A3; X=untag(A3); 
  I1.X3=contract(X,'321*',X,'123');

  SS=[ contract(A3,'!4*',{S1,'-op:a','*',{S1,'-op:c',A3}})
       contract(A3,'!4*',{S1,'-op:a','*',{S1,'-op:b',A3}})
       contract(A3,'!4*',{S1,'-op:b','*',{S1,'-op:c',A3}}) ]; % Sb.Sc J3'

  Hloc=QSpace(1,nJ);
  for i=1:nJ
     Hloc(i) = j3(i,1)*SS(1) + j3(i,2)*SS(2) + j3(i,3)*SS(3);
  end

  Sa=contract(A3,'!4*',{A3,S1,'-op:a'});
  Sb=contract(A3,'!4*',{A3,S1,'-op:b'});
  Sc=contract(A3,'!4*',{A3,S1,'-op:c'}); Sloc=[Sa Sb Sc];

  A6=setitags(QSpace(getIdentityQS(A3,4,A3,4)),{'s1','s2','s12'});
  H6 = contract(A6,'!3*',{Sa,'-op:s1','*',{Sa,'-op:s2',A6}}) ...
     + contract(A6,'!3*',{Sc,'-op:s1','*',{Sc,'-op:s2',A6}});
  if jc, H6 = H6 ...
  + jc*contract(A6,'!3*',{Sb,'-op:s1','*',{Sb,'-op:s2',A6}});
  end

  H4=QSpace(contractQS(A6,'!12',{H6,A6,'*'}));
  [Sleg,Sl_]=splitH4_SdotS(H4,'-v');

  e=norm(Sleg-Sl_);
  if e>1E-12, error('Wb:ERR','unexpected setting (e=%g)',e); end

  Eloc=getIdentity(A3,4);
  untag(Sleg,Hloc,Eloc,Sloc);

  if p.xflag, q=iff(all(p.J(:,2:3)>=2),'==>','WRN');
     wblog(q,'projecting out S=3/2 multiplet space');

     project_out_q3(Sleg,Hloc,Eloc,Sloc);

     if p.xflag>1
     end
  end
end

% -------------------------------------------------------------------- %
function project_out_q3(varargin)

   for l=1:nargin
      X=varargin{l}; n=inputname(l); 
      if isempty(n), error('Wb:ERR','invalid usage'); end
      for j=1:numel(X)
         q=[X(j).Q{1}, X(j).Q{2}]; i=find(sum(q<=1,2)==2);
         X(j)=getsub(X(j),i);
      end
      assignin('caller',n,X);
   end

end

% -------------------------------------------------------------------- %

