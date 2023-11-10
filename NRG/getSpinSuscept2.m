function [chi,Iout]=getSpinSuscept2(Simp,varargin)
% function [chi,I]=getSpinSuscept2(Simp [,NRG,opts])
%
%    calculate mixed spin susceptibility derived from Matsubara
%    expression for <Simp||Stot> = <Simp*Stot>. Here Simp is the spin
%    operator at the impurity represented in the local basis of A0.
%
% Wb,Jun19,13

% adapted from getSpinSuscept2_iter.m
% see also rnrg_tst_chi2.m

  getopt('init',varargin);
     vflag=getopt('-v'); if ~vflag
     qflag=getopt('-q'); else qflag=0; end
     TT=getopt('T',[]);
     iS=getopt('iS',[]);
  args=getopt('get_remaining'); nargs=numel(args);

  NRG=[]; in_mem=0;
  if nargs
     if ischar(args{1})
        NRG=args{1};
     elseif iscell(args{1}) && numel(args{1})==2 && ...
        isfield(args{1}{1},'AK') && isfield(args{1}{2},'EScale')
        NRG=args{1}; Inrg=NRG{2};
     end
     if ~isempty(NRG), args(1)=[]; nargs=nargs-1; end
  end
  if isempty(NRG), NRG=[getenv('LMA') '/NRG/NRG']; end

  if nargs
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if ischar(NRG)
     f=sprintf('%s_info.mat',NRG);
     if ~exist(f,'file')
        wbdie('invalid NRG data (%s)',NRG);
     end
     Inrg=load(f);
  end

  if isempty(TT), mt=4;
     TT=(250*Inrg.EScale(end))*Inrg.Lambda.^((0:mt-1)/mt);
  end

  nT=numel(TT); rw=zeros(Inrg.N,nT);

  i=find(TT<=0); n=numel(i);
  if n
     if n>1, wbdie('invalid T (nan may only be set once)'); end
     rw(:,i)=repmat(Inrg.rhoNorm(:,2),1,n);
     TT(i)=Inrg.rhoT;
  end

  if n<nT
     if vflag, o={'-v'}; else o={}; end
     if ~qflag, wblog(' * ','get FDM weights for %g temperatures',nT); end

     [rw,TT]=getWeightsNRG(NRG,'T',TT,o{:});
     if n
        e=norm(rw(:,i)-repmat(Inrg.rhoNorm(:,2),1,n));
        if norm(TT(i)-Inrg.rhoT)~=0 || e>1E-12
           wbdie('FDM weigth inconsistency (@ %.3g)',e); 
        end
     end
  end

  i=find(rw(end,:)>0.01);
  if numel(i)
     j=find(rw(end,:)==max(rw(end,:)));
     if numel(j)>1, j=j(find(TT(j)==min(TT(j)),1)); end
     wblog('WRN','got too small a temperature (T=%.3g @ %.3g) !??',...
     TT(j),rw(end,j));
  end

  N=size(rw,1); chi=zeros(nT,N);

  if ~qflag, wblog(' * ',...
     'calculate spin susceptibility <Simp||Stot> ...',nT); end

  wblog('WRN','USING STOT_SIGNS');

  for k=1:N
     if vflag, fprintf(1,...
       '\r   k=%2g/%g (T=%.3g .. %.3g; %g values) ... \r',...
        k,N,min(TT),max(TT),nT);
     end
     if ischar(NRG)
          q=load(sprintf('%s_%02g.mat',NRG,k-1),'AK','AD','HD','ES');
     else q=NRG{1}(k);
     end

     if k>1
        if ~isempty(q.AD.data)
           SDD=contractQS(q.AD,'!2*',{q.AD,SKK});
        end
        if k<N
           SKK=contractQS(q.AK,'!2*',{q.AK,SKK});
        end
        if 0
           if k==2
                STOT=SKK;
           else STOT=contractQS(q.AK,'!2*',{q.AK,STOT});
           end

           t=q.AK.info.itags{end};
           STOT=QSpace(STOT)+contractQS(q.AK,'!2*',{q.AK,Simp,['-op:' t]});
        end
     else
        if ~isempty(q.AD.data), wbdie(...
           'got NRG truncation already with first Wilson site !?');
        end
        d=getDimQS(q.AK);
        if d(end,1)>1
             SKK=contractQS(q.AK,'!2*',{q.AK,Simp,'-op:L00'}); kimp='L';
        else SKK=contractQS(q.AK,'!2*',{q.AK,Simp,'-op:s00'}); kimp='s';
        end

        if isempty(SKK.data)
           wbdie('invalid Simp contraction (empty)');
        end
     end

     if isempty(q.AD.data), continue; end

     for it=1:nT
        if vflag
        end

        if ~isempty(SDD.data)
           [ctk,iS] = getSpinSuscept2_iter(SDD,q.HD, TT(it)/q.ES, iS);
           chi(it,k) = (rw(k,it)/q.ES) * ctk;
        end
     end
     if vflag, fprintf(1,'\r%50s\r',''); end
  end

  if nargout>1
     Iout=add2struct('-',chi,rw,TT,kimp,iS);
  end

  chi=reshape(sum(chi,2),size(TT));

end

% -------------------------------------------------------------------- %

function [chi,iS]=getSpinSuscept2_iter(Simp,H,T,iS)
% function [chi,iS]=getSpinSuscept2_iter(Simp,H,T,iS)
%
%    calculate (contribution to) mixed spin susceptibility
%    <Simp||Stot> = <Simp*Stot>. Here Simp is the matrix elements
%    of S_imp represented in the effective basis of given NRG
%    iteration, H contains the energy of each state, which
%    is required for the construction of the thermal density
%    matrix at temperature T.
%
%    For the construction of Stot, the position
%    of the SU(2) spin symmetry in info.qtype is required.
%    if iS is not specified, this routine tries to
%    determine iS from the structure of Simp.
%
% Wb,Jun19,13

  if nargin==1 && isa(Simp,'QSpace')
     chi=get_Stot(Simp); return
  elseif nargin<3 || nargin>4
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  Simp=struct(Simp);

  if nargin<4 || isempty(iS)
     if numel(Simp.Q)<3, wbdie([...
       'cannot determine iS from abelian Simp (Sz only)\n' ...
       '=> must specify iS upon input']);
     end

     q=std((Simp.Q{1}-Simp.Q{2}).^2,[],1); i=find(q); n=numel(i);
     if n>1
        s=sprintf('got %g symmetries with varying q-labels (%g)',n);
        if numel(find(q>0.99*max(q)))>1, wbdie(s);
        else wblog('WRN',s); end
     end
     if max(q)>1E-12
          iS=find(q==max(q));
     else iS=(find(Simp.Q{3}(1,:)==2));
     end
     if numel(iS)~=1
        iS, wbdie('failed to determine iS');
     end
  elseif numel(iS)~=1
     iS, wbdie('invalid iS');
  end

  if numel(Simp.Q)>2
     if norm(diff(Simp.Q{3},[],1))>1E-8 || norm(Simp.Q{3}(:,iS)-2)>1E-8
        wbdie('invalid Simp (expecting iROP here; %g)!',iS);
     end
  end

% -------------------------------------------------------------------- %
% build Stot-operator *implicitly* from Simp
%
% (1) Simp already has all the correct CGC spaces!
% (2) note that Stot is *diagonal* in given effective state space!
%     that is, Q{1}==Q{2} must be equal throughout! (this is
%     justified by the fact S_+ only affects the Sz label of
%     *given* multiplet!)
%
%  => can reduce Simp to block-diagonal
%     (only these enter in the final contraction).
% -------------------------------------------------------------------- %

  i=find(sum(abs(Simp.Q{1}-Simp.Q{2}),2)<1E-8);
  Simp=getsub(QSpace(Simp),i);

  dS=getDimQS(Simp); dS(:,end+1:3)=1;
  dop=dS(end);

  aflag=iff(diff(dS(:,3)),0,1);

  if ~isa(H,'QSpace'), H=QSpace(H); end
  if isdiag(H,'-d')~=2
     wbdie('invalid Hamiltonian (expected in diagonal format!)');
  end

  [i1,i2,I]=matchIndex(Simp.Q{1},H.Q{1});
  if ~isempty(I.ix1)
     wbdie('Simp appears incomplete/incompatible w.r.t. given H!');
  end

  if T<=0, wbdie('invalid temperature (%g)',T); end
  beta=1/T;

  R=getrhoQS(H,beta);
  X=contractQS(R,2,Simp,1); if isempty(X.data), chi=0; return; end

  ns=size(X.info.cgr,2);

  for i=1:numel(X.data)
     S=X.Q{1}(i,iS)/2;

   % NB! need matrix elements of Sz_tot or Stot only
   % since Simp is already present; the CGC contraction
   % of Simp.Stot eventually yields 1. by normalization
   % convention on CGCs within QSpace // Wb,Aug15,17
   % for red. matrix element of Stot see e.g. Hanl14, App. C.1]
     if aflag
        a=S;
     else
        a=sqrt(S*(S+1));
        if S<=1, a=-a; end
     end

     m=1;
     w=1;

     for j=1:ns, c=X.info.cgr(i,j);
        q=reshape(mpfr2dec(c.cgw),[],1);
        if norm(q'*q-1)>1E-12
           wblog('WRN','got norm(cgw)=%4g !?',q'*q); end
        w=w*sum(q);
        m=m*double(c.size(1));
     end
     if norm(w-1)>1E-12, wblog('WRN','got w=%4g !?',w); end

     X.data{i}=(a*w*sqrt(m))*trace(X.data{i});
  end

  chi=(beta/dop)*sum([X.data{:}]);

end

% -------------------------------------------------------------------- %
function [Stot,iS]=get_Stot(Simp,iS)

  if nargin<2 || isempty(iS), iS=[];
     if numel(Simp.Q)==3
        q=uniquerows(Simp.Q{3}); i=find(q(1,:));
        if size(q,1)==1 && numel(i)==1, iS=i; end
     end
  end

  Stot=Simp(end,1);
  i=find(sum(abs(Stot.Q{1}-Stot.Q{2}),2)<1E-8);
  Stot=getsub(Stot,i);

  ns=size(Stot.info.cgr,2);
  Stot=struct(Stot);

  for i=1:numel(Stot.data)
     S=Stot.Q{1}(i,iS)/2; a=sqrt(S*(S+1));
     s=min(size(Stot.data{i}));
     m=1;

     for j=1:ns
        c=Stot.info.cgr(i,j);
        m=m*double(c.size(1));
     end
     Stot.data{i}=diag(repmat(a*sqrt(m),1,s));
  end
  Stot=skipzeros(QSpace(Stot));

end

% -------------------------------------------------------------------- %

