function [mps,Iout,IL,E3]=MPS_add(varargin)
% function [mps,Iout]=MPS_add(mps1[,fac1],mps2[,fac2],... {,opts})
%
%    Add set of MPSs based on the respectors factors fac_k
%
% Options
%
%    'nsw', ...   number of sweeps used for MPS fitting (3)
%    'rtol',,..   truncation threshold on svd^2 values (1E-24);
%    'Nkeep',..   maximum MPS dimension to keep (default -1,
%                 i.e., truncate based on rtol)
% Wb,Feb11,25

% adapted from MPO_add // Wb,Feb11,25

% ensure that simple operations such as subtracting a constant etc.
% are always accurate irrespective of the input MPS
% --> do not set Nkeep (i.e., keep Nkeep=-1)
  Nkeep=-1; rtol=1E-24; nsw=3; % default values

  t0=tic();

  opts={}; gotE3=0;
  nargs=numel(varargin);
  if nargs
     if iscell(varargin{end}) || isstruct(varargin{end})
        opts=varargin{end};
        nargs=nargs-1;
     else
        for i=1:nargs
           if ischar(varargin{i})
              opts=varargin(i:end);
              nargs=i-1; break
           end
        end
     end
  end
  
  if nargs<2, wbdie('invalid usage (%d args)',nargs); end
  iM=[]; l=0; fac=zeros(0,2);
  for i=1:nargs
     if isQSpace(varargin{i}), l=l+1; iM(l)=i;
        si=size(varargin{l});
        if l==1, s1=si;
        elseif ~isequal(s1,si), wbdie('inconsistent MPS sizes'); end
     elseif isnumber(varargin{i})
        if ~l, wbdie('invalid usage (first entry must be MPS)'); end
        fac(l,1)=varargin{i};
        fac(l,2)=fac(l,2)+1;
     else wbdie('invalid usage (arg #%d)',i); 
     end
  end

  if any(fac(:,2)>1), wbdie(...
    'invalid usage (MPS factor specified multiple times?)'); end
  fac(end+1:l,:)=0;
  i=find(fac(:,2)==0); fac(i,1)=1;
  fac=fac(:,1).';

  if s1(1)~=1, s1, wbdie('unexpected MPS QSpace array size'); end
  N=prod(s1); nterms=l; % number of states

  MPS=cat(1,varargin{iM}); 

  vflag=1; kflag=0;

  if ~isempty(opts)
     if isstruct(opts) % convert to regular cell
        if isfield(opts,'sweep'), I=opts.sweep; else I=opts; end
        opts={}; % safeguard: nsw -> nsw_mps
        for f={'nsw_mps','Nkeep','rtol'} % fields
           if isfield(I,f{1}), x=getfield(I,f{1});
              if ~isempty(x), opts=[opts, {f{1},x}]; end
           end
        end
     end

     getopt('init',opts);
        if     getopt('-v'), vflag=2;
        elseif getopt('-q'), vflag=0; end
        kflag=getopt('-k');

        nsw  = getopt('nsw',  nsw);
        Nkeep= getopt('Nkeep',Nkeep);
        rtol = getopt('rtol', rtol);

     getopt('check_error');
  end

  stol=sqrt(rtol); if rtol<0, wbdie('invalid rtol=%g',rtol); end
  i=find(abs(fac(:,1))<min(rtol,1E-20));
  if ~isempty(i), q=[ numel(i), max(abs(fac(i,1))) ];
     wblog('WRN','skipping %d terms (|fac|<%.3g)',q); 
     fac(i)=[]; 
     MPS(i,:)=[];
  end

  E1=QSpace(nterms,N); % local identity
  for l=1:nterms
     for k=1:N
        E1(l,k)=getIdentityQS(MPS(l,k),3);
     end
     if l>1
      % safeguard to ensure complete local state space
        E1(1,k)=E1(1,k)+E1(l,k);
        if l==nterms, E1(1,k)=getIdentity(E1(1,k)); end
     end
  end
  E1=E1(1,:);

  osw={'stol',stol,'Nkeep',Nkeep}; % stol is taken absolut

% initialize output with largest MPS
  s=zeros(nterms,1);
     for l=1:nterms, s(l)=sizeof(MPS(l,:)); end
     i=find(s==max(s),1);
  mps=MPS(i,:); iM=i;

% initialize: normalize mps and obtain XR
  XL=QSpace(nterms,N); XR=XL;
  for k=N:-1:2
   % make sure mps is R->L orthonormalized
     [mps(k),X,Io]=orthoQS(mps(k),1,'<<',osw{:});
     mps(k-1)=contract(mps(k-1),X,[1 3 2]);

   % calculate overlaps XR
     for l=1:nterms
        if k<N, X_=XR(l,k+1); else X_='bdry'; end
        XR(l,k)=update_overlap(1,l,mps(k),MPS(l,k),X_);
     end
  end

  tt=get_time(t0);
  istr=''; xpo=[]; xol=[];

  for isw=0:nsw % isw=0 only computes initial overlap
     if nargout>1, mps_last=mps; end

     if isw % <----- outer bracket!

     for dk=[1 -1] % sweep direction
        if dk>0, ksw=1:N-1; else ksw=N:-1:2; end
        for k=ksw
           if dk>0
                k1=k; k2=k+1;
           else k1=k-1; k2=k; end

           A1=mps(k1); r1=rank(A1);
           A2=mps(k2); r2=rank(A2);

           [tb,c,m]=getitags(mps(k1),r1-1); % c=conj, m=mark
         % whether to apply mark on bond [eventually using mod(m,2)]
         % i.e., if ~c && dk || c && dk<=0
           if xor(c~=0,dk>0), m=m+1; end

         % tl/tr: L/R itags for bond tensor; tb: final itag on bond
           tl=[tb 'l'];
           tr=[tb 'r']; s='''';

           if mod(m,2)
                tl=[tl s]; tb=[tb s];
           else tr=[tr s]; end

         % switch 2-site -> bond picture
           if r1==3 % typical case (within MPS)
              B1=getIdentity(A1,1,E1(k1),1,tl,[1 3 2]);   p1=[2 1 3];
           elseif k1==1 && r1==2
            % [B1,X1]=orthoQS(A1,1,'<<','itag',tl,osw{:}); p1=[];
            % NB! m01 may not be complete => consider full operator Id!
              B1=permute(E1(k1),[2 1]); setitags(B1,1,tl); p1=[];
           else
              wbdie('invalid usage (unexpected rank r1=%d at k=%d)',r1,k1);
           end

           if r2==3 % same for the right tensor A2->B2
              B2=getIdentity(A2,2,E1(k2),1,tr,[3 1 2]);
           elseif k2==N && r2==2
              B2=permute(E1(k2),[2 1]); setitags(B2,1,tr);
           else
              wbdie('invalid usage (unexpected rank r2=%d at k=%d)',r2,k2);
           end

           X1=contract(B1,'*',A1); 
           X2=contract(B2,'*',A2);

         % old X12 (actually not required since replaced right below!)
           X12_=contract(X1,X2);
           r=rank(X12_); if r~=2
             wbdie('invalid contractions (got rank-%d QSpace X)',r); end

         % update overlaps
           for l=1:nterms
              if k1>1, X1=XL(l,k1-1); else X1='bdry'; end
              if k2<N, X2=XR(l,k2+1); else X2='bdry'; end
              XL(l,k1)=update_overlap(2,l,B1,MPS(l,k1),X1);
              XR(l,k2)=update_overlap(1,l,B2,MPS(l,k2),X2);
           end

         % new X12 -> mps(k,k+1)
           for l=1:nterms
              Q=fac(l)*contract(XL(l,k1),XR(l,k2));
              if l>1 X12=X12+Q; else X12=Q; end
           end

           if dk>0
              [U,X,Il(k1,1)]=orthoQS(X12,2,'<<','itag',tb,osw{:});
              mps(k1)=contract(U,B1,p1);
              mps(k2)=contract(X,B2);
           else
              [U,X,Il(k1,2)]=orthoQS(X12,1,'<<','itag',tb,osw{:});
              mps(k1)=contract(X,B1,p1);
              mps(k2)=contract(U,B2);
           end

         % update overlaps (L or R depending on dk)
           if dk>0
              for l=1:nterms
				 if k1>1, X1=XL(l,k1-1); else X1='bdry'; end
				 XL(l,k1)=update_overlap(2,l,mps(k1),MPS(l,k1),X1);
              end
           else
              for l=1:nterms
                 if k2<N, X2=XR(l,k2+1); else X2='bdry'; end
                 XR(l,k2)=update_overlap(1,l,mps(k2),MPS(l,k2),X2);
              end
           end
        end % ksw (k range for half-sweep)
     end % dk (sweep direction)

        if nargout>1
         % this contains |mps - mps_last|^2 when combined with xol(isw-1,:)
           xol(isw,1:2)=[ MPS_overlap(mps_last,mps), MPS_overlap(mps) ];
        end

     else k1=1; k2=2;
     end % if isw

   % being at k=1, i.e., (k1,k2)=(1,2)
   % simply also compute overlap xpo (for info purposes only)
     k=k1;
     for l=1:nterms
        XR(l,k)=update_overlap(1,l,mps(k),MPS(l,k),XR(l,k+1));
        xpo(isw+1,l)=getscalar(XR(l,k));
     end

     if ~isw, continue; end

     s2=getdatafield(Il,'svd2tr','-0'); % -0 => default value: 0
     ds2max=max(s2(:));
     nk=getdatafield(Il,'Nkeep','-0'); Nkept=max(nk(:));

     if nargout>2
        q=struct('svd2tr',s2,'Nkeep',nk,'ds2max',ds2max,'Nkept',Nkept);

        [k,j]=find(nk==max(nk(:)),1); % first iteration with largest Nkept
        Iq=Il(k,j); Iq.k=[k j, N]; q.Ik=Iq;

        [k_,j_]=find(s2==ds2max,1); % iteration with largest svd2tr
        if ~isequal([k j],[k_ j_])
           Iq=Il(k_,j_); Iq.k=[k_ j_, N]; q.Ik(2)=Iq;
        end

        q.xpo=xpo(isw+1,:);
        IL(isw)=q;
     end

     if ds2max<rtol
        istr=sprintf('converged at isw=%d/%d @ svd_tr <= %.2g / %.2g',...
          isw,nsw,sqrt(ds2max),stol);
        if vflag>1, wblog('-->',istr); end
        break
     end
     tt(end+1)=get_time(t0);
  end % n_sweeps

% copy itags from from first input MPS
  for k=1:N
      mps(k)=setitags(mps(k),MPS(1,k));
  end

  if nargout>1
     se=zeros(size(Il));            % normalize! --v
     for i=1:numel(Il), se(i)=SEntropy(Il(i).svd,'-n'); end

     if isempty(istr), istr=sprintf(...
        'target MPSs at svd2tr <= %.3g / %g (nsw=%d)',ds2max,rtol,isw);
     end

     Iout=add2struct('-',istr,'converged=0',fac,... % NORM_OPS
         N,nterms,iM,Nkeep,Nkept,rtol,'mps_stol=0','Dtot?',se,ds2max,isw,xol);
     Iout.converged = (ds2max<rtol);
     Iout.mps_stol=osw{2};
     if isw<nsw, Iout.isw(2)=nsw; end
   % if isvar('Dtot'), Iout.rtol=Iout.rtol/Dtot; end

     if nargout<3
          Iout.xpo=xpo; 
     else Iout.initial_xpo=xpo(1,:); % initial overlap only
     end
     % else contained in IL anyways
  end

  if kflag, wbstop; end
  if ~vflag, return; end

end

% -------------------------------------------------------------------- %
% see also HAM/initRho.m

function S=get_time(t0)

   S=dbstack();
   S=struct('line',S(2).line,'time',toc(t0));

end

% -------------------------------------------------------------------- %
% upate overlap // Wb,Feb11,25

function X=update_overlap(odir,iM,mps,MPS,X_)

   if odir, ic=sprintf('!%d*',odir); else ic='*'; end

 % differentiate bond itags of input MPSs
   if isnumeric(iM), iM=char('a')+(iM-1); end
   trep=['--itag:s/^([A-R]+\d+)/$1' iM '/'];
   % include K,D,H, etc., exclude like S##

   gotX=(nargin>4 && ~isequal(X_,'bdry'));
   if gotX
      MPS={MPS,trep,X_};
   end

   X=contract(mps,ic,MPS);

   if ~gotX % need to fix itag at boundary
      i=find(trep=='/'); % got 's/*/*/'
      p=trep( i(1)+1 : i(2)-1 );  % pattern
      r=trep( i(2)+1 : i(3)-1 );  % replace
      X=itagrep(X,2,p,r);
   end
end

% -------------------------------------------------------------------- %

