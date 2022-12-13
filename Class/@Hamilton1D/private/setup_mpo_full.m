function [HAM]=setup_mpo_full(HAM,varargin)
% function [HAM]=setup_mpo_full(HAM [,HS [,stype, opts]])
%
%    Automated routine to generate full MPO (HAM.mpo)
%    for the Hamiltonian described by HS together with
%    HAM.ops, HAM.oez, and HAM.info.stype (site type).
%
%    By default, if hermitian conjugate terms are required
%    the full Hamiltonian is Hfull = H_ + H_^\dagger
%    where H_ is encoded in the returned MPO HAM.info.mpo (!)
%
% Wb,Nov17,20

   getopt('INIT',varargin);
      nsw  =getopt('nsw',[]);
      kflag=getopt('-k');
      vflag=getopt('-v');
      tflag=getopt('-t');
      Fflag=getopt('-F');
      addHc=getopt('--Hc');
   args=getopt('get_remaining'); nargs=length(args);

   gotH=0;

   if nargs>2, wberr('invalid usage'); end
   if nargs, HS=args{1};
   else
      if isfield(HAM.info,'HS')
           HS=HAM.info.HS; gotH=1;
      else HS=HAM.info.HH; end
   end
   if nargs>1, stype=args{2};
   else
      if isfield(HAM.info,'stype'), stype=HAM.info.stype;
      else stype=[]; end
   end

   if isnumeric(HS) && size(HS,2)==5, HS=convert_HH(HS); gotH=0; end
   if ~isstruct(HS) || ~isfield(HS,'iop') || ~isfield(HS,'hcpl')
      wberr('invalid HS (struct array expected)');
   end

   nH=numel(HS); mH=zeros(nH,1);

   for i=1:nH, q=HS(i).iop;
      if ~isnumeric(q) || ndims(q)>2 || isempty(q) || size(q,2)~=2
          disp(HS(i)); wblog('ERR','invalid HS(%g)',i); 
      end
      mH(i)=size(HS(i).iop,1);
   end

   q=zeros(1,nH); for i=1:nH, q(i)=norm(HS(i).hcpl); end
   i=find(q<1E-15);
   if ~isempty(i)
      if numel(i)==nH
         wbdie('invalid usage (got all-zero Hamiltonian terms)');
      end
      HS(i)=[]; q=[numel(i),norm(q)];
      wblog(' * ','eliminating %g zero-terms in Hamiltonian (@ %.3g)',q);
      nH=numel(HS);
   end

   q=cell(nH);
   for i=1:nH, [q{i},is]=sort(HS(i).iop(:,1));
      if numel(unique(q{i}))~=numel(q{i})
         disp(HS(i)), wberr('got non-unique site index in HS(%g)',i);
      end
      HS(i).iop=HS(i).iop(is,:);
   end

   k=unique(cat(1,q{:})); L=k(end); ntype=0;
   if ~isequal(k,(1:L)')
      wberr('invalid usage [HS is not addressing all sites !?]');
   end

   ndig='2';

   oez=getdatafield(HAM.oez,'op');
   ops=getdatafield(HAM.ops,'op'); nops=size(ops,1);

   ntype=size(oez,2); q=[ntype,size(ops,2)];
   if diff(q), wberr('invalid usage (ntype=%g/%g !?)',q); end

   q=cell(nH); for i=1:nH, q{i}=HS(i).iop(:,2); end
   q=cat(1,q{:}); q=[min(q), max(q), nops];
   if q(1)<1 || q(2)>q(3), wberr('HS.iop out of bounds (%g,%g/%g)',q); end

   if ~isempty(stype)
      if unique(stype)==1, stype=[];
      else
         q=[numel(stype), L];
         if diff(q), wberr('stype size mismatch (l=%g/%g)',q); end
         ntype=max(stype);
         if ntype<=1, wberr('got stype with ntype=%g',ntype); end
      end
   end

   q=[ntype, size(ops,2), size(oez,2)];
   if diff(q(2:3)), wberr(...
     'HAM.ops and HAM.oez size mistmatch (%g/%g)',q(2:3)); end
   if q(1)>q(2), wberr('stype out of bounds (%g/%g)',q(1:2)); end

   if isempty(oez), wberr('invalid usage (HAM.oez not yet setup)');
   else
     q=size(oez,1);
     if q(1)>2, wberr('invalid HAM.oez (unexpected size %d/2)',q(1)); end
   end

   if ~isempty(stype), s=sprintf(', %g',stype);
      wblog('NB!','got stype = %s',s(3:end)); end

   if ~gotH, HAM.info.HS=HS;
      if ntype, HAM.info.stype=stype; end
   end

   if isempty(stype), stype=ones(1,L); end
   for i=1:nH, HS(i).iop(:,3)=stype(HS(i).iop(:,1)); end

   iop=cat(1,HS.iop);

   qops=ones(size(ops));
   for i=1:numel(ops)
      if isempty(ops(i)) || normQS(ops(i))<1E-6 || ...
         isempty(ops(i).Q), qops(i)=0;
      end
   end

   qoez=ones(size(oez));
   for i=1:numel(oez)
      if isempty(oez(i)) || normQS(oez(i))<1 || ...
         isempty(oez(i).Q), qoez(i)=0; end
   end

   if any(qoez(:)==0)
      wberr('invalid usage (HAM.oez not yet setup !?)'); end

   for j=1:size(oez,2)
      if ~isIdentityQS(oez(1,j))
      wberr('invalid oez (expecting E as first operator)'); end
   end

   isFerm=getdatafield(HAM.ops,'fermionic');
   isferm=max(isFerm,[],1); i=find(isferm);
   if ~isempty(i)
      if size(oez,1)~=2 || any(~qoez(2,i))
         wberr('invalid oez (z-ops required based on HAM.ops)');
      end
   end

   if size(qoez,1)>1 && all(qoez(2,:)==0)
      wblog('mpo','unsetting z-ops since not required'); 
      oez=oez(1,:);
   end
   noez=size(oez,1);

   isHconj=getdatafield(HAM.ops,'hconj');
   if any(isHconj(find(qops)))
        addHc=bitor(addHc,2);
   else addHc=0; end

   HR=zeros(nH,5);
   for i=1:nH
      q=HS(i).iop(:,1); q=[ min(q), max(q), numel(q) ];
      HR(i,:)=[ q(1:2), q(2)-q(1), q(3), HS(i).hcpl ];
      HS(i).nops=q(3);
   end
   iloc=find(HR(:,3)==0)';
   mops=max(HR(:,4)); 
   hcpl=[HS.hcpl];

   [q,I,D]=uniquerows(HR(:,[4 3]));
   if isempty(nsw)
      nsw=max(3,numel(I));
   end
   q(:,2)=q(:,2)+1;

   wblog(' * ','Hamilton setup HS (L=%g, max.range=%g)',L,max(q(:,2)));

   isb=isbatch(); s='      #terms  #ops  range     hcpl';
   if isb, fprintf(1,'\n%s\n',s);
   else printfc('\n\e[38;5;8m %s \e[0m\n',s); end

   for i=1:size(q,1), ss={'',''};
      if     q(i,2)==1, ss{1}='onsite';
      elseif q(i,2)==2, ss{1}='NN';
      elseif q(i,2)==3, ss{1}='NNN'; end

      h=unique(hcpl(I{i}));
      if numel(h)<6, s=vec2str(h,'%.4g','sep',', ');
         if length(s)<50, ss{2}=s; end
      end
      if isempty(ss{2})
         s=vec2str(h([1,end]),'%.4g','sep',' ..');
         ss{2}=[s, sprintf(' (%g)',numel(h)) ];
      end

      if i>1 && diff(q(i-1:i)) && q(i)>2, fprintf(1,'\n'); end
      fprintf(1,' %12d  %4d  %2d %-6s  %3s\n',D(i),q(i),q(i,2),ss{1},ss{2});
   end
   fprintf(1,'\n');

   if mops>2
      wblog('NB!','got %g-site operator terms',mops);
   elseif mops~=2, wberr('got mops=%g',mops); end

   q=uniquerows(iop(:,2:end));
   i=q(:,1) + (q(:,2)-1)*size(ops,1);
   if any(~qops(i)), wberr('got referenced zero ops'); end

   q=cat(1,HS(iloc).iop); q=uniquerows(q(:,2:end));
   for i=1:size(q,1)
      X=HAM.ops(q(i,1),q(i,2)); r=numel(X.op.Q);
      if r~=2, wberr('got non-scalar local Hamiltonian term'); end
      if X.fermionic, wberr('got fermionic local Hamiltonian term'); end
      if ~X.hconj, e=normQS(X.op-X.op'); if e>1E-8
         wberr('got non-hermitian local Hamiltonian term'); end
      end
   end

   [x,I,D]=uniquerows(HR(:,4));
   for i=1:size(x,1)
      if x(i)<2, continue; end
      q=cat(3,HS(I{i}).iop); s=size(q);
      q=uniquerows(reshape(q(:,2:end,:),[],D(i))');
      q=reshape(q',[s(1), s(2)-1, size(q,1)]);

      j=q(:,1)+nops*(q(:,2)-1); opl=[HAM.ops(j).op]; l=0;
      n=numel(opl); if n>4, wberr('unexpected #ops = %g',n); end
      for k=1:n
         r=numel(opl(k).Q); 
         if r<3, continue; elseif ~l, l=k; end

         t=sprintf('s%g',k); t={t,t}; t(3:r)={'op'};
         opl(k)=setitags(opl(k),t);

         if l>0 && l<k, X=contract(opl(l),'*',opl(k)); l=-2;
         elseif l<0, X=contract(X,opl(k)); l=l-1; end
      end
      q=norm(X); if q<1E-8, wberr(...
         'invalid %g-site interaction term (nrm=%g)',n,q); end
      q=numel(X.Q); if q~=2*(-l), wberr(...
         'invalid %g-site interaction term (contracts to rank %g)',n,q);
      end
   end 

   IH.HR=HR;

   if mops<=2
      H=zeros(1,L); if ~isempty(iloc)
         k=cat(1,HS(iloc).iop); k=k(:,1);
         H(k)=[HS(iloc).hcpl]/2; % see (H+H') below!
      end
      H=diag(H); e=0;

      l2=1:nH; l2(iloc)=[];
      for l=l2, q=HS(l).iop; h=HS(l).hcpl;
         if numel(q)~=6, disp(HS(l)), wberr('unepected iop'); end

         k1=q(1,1); k2=q(2,1); i=q(1,2); j=q(1,3);
         if ~HAM.ops(i,j).hconj, h=h/2; end % see (H+H') below!
         if H(k1,k2) || H(k2,k1), e=e+1; end
         H(k1,k2)=h;
      end
      H=H+H';

      IH.H=H; if e
      IH.H_WRN='(got multiple interactions for the same pair of sites)'; end
   else
      IH.H=[];
      IH.H_WRN='(got multi-site terms)';
   end

 % build *normalized* local operator basis which can be used to fuse
 % local mpo bra/ket indices into a `fat local super-state space'
 % => ||Id|| = 1. (rather than d^L)
 % => for consistency, must also scale each local operator
 %    in HAM.ops or HAM.oez by 1/sqrt(d)
 % => eventually, need to re-multiply each mpo(k) with 1/dfac

   n=max(1,ntype); dloc=zeros(1,n);
   for i=1:n, q=getDimQS(HAM.oez(1,i).op); dloc(i)=q(end); end
   dfac=1/sqrt(max(dloc));

   for i=1:numel(oez), oez(i)=dfac*oez(i); end

   for i=1:numel(ops), q=dfac*ops(i);
      if rank(q)<3 && ~isscalarop(q), q=makeIrop(q); end
      d=getqdir(q);
      if numel(d)>=3 && d(3)>0 || numel(d)==4 && d(4)<0
         wberr('invalid operator qdir');
      end
      ops(i)=q;
   end

   EE=oez(1,:); E0=EE; 
   for j=1:size(EE,2)
      u=getIdentityQS(EE(j),'-0');
      E0(j)=contractQS(getIdentityQS(u,1,u,2),2,u,'2*',[1 3 2]);

      if 1
         x2=QSpace; X2=QSpace;

         for i=1:size(oez,1), q=appendSingletons(oez(i,j),'  -');
            Q=contract(q,'12*',E0(j),'12');
            x2=x2+contract(Q,'1*',Q,1);
         end

         l=i+1;
         for i=1:size(ops,1), q=appendSingletons(ops(i,j),'  -+');
            if isFerm(i,j), Q=contract(oez(2,j),2,q,1);
            else Q=q; end

            Q=contract(Q,'12*',E0(j),'12');
            x2=x2+contract(Q,'12*',Q,'12'); l=l+1;  

            if isHconj(i,j)
               Q=contract(q,'21',E0(j),'12');
               x2=x2+contract(Q,'12*',Q,'12');
            end
         end

         [ee,Ie]=eigQS(x2,'Rtrunc',1E-9); X2(j)=x2;
         for i=1:numel(Ie.AK.data)
            q=Ie.AK.data{i}; [n,m]=size(q); if n<2, continue; end
            if n==m
               Ie.AK.data{i}=eye(n);
            elseif m>1
               for j=1:m-1
                  u=zeros(n,1); l=find(sum(q.^2,2)>1E-8,1); u(l)=1;
                  u=q*(q'*u); x=norm(u); if x<1E-12, break; end
                  U(:,j)=u/x;
                  [q,s,vd]=svd(q-U*(U'*q),'econ'); s=diag(s);
                  q=q(:,find(s>1E-12));
               end
               if x>1E-12 && size(q,2)==1, Ie.AK.data{i}=[U,q];
               else wblog('WRN',...
                 'got x=%.3g with %g vectors in q',x,size(q,2));
               end
            end
         end
         EE(j)=rmitags(contract(E0(j),3,Ie.AK,1));
      else
         EE(j)=rmitags(E0(j));
      end
   end
   EE_=EE; EE=EE(1,stype);

   for i=1:numel(ops)
    % NB! in interactions terms the dagger is always applied
    % to left-most (the `1st') operator; since irops have an
    % outgoing irop index, by convention, this means that 
    % for consistency, operator chains in the Hamiltonian
    % must point right-to-left // OPS_POINT_R2L
    % => keep origional irop index with the left side
    % => index order convention for ops // INDEX_ORDER_MPO
    %    bottom,top,left,right (btlr) where the original operator
    %    index points to the left(!) at index position 3
      if numel(ops(i).Q)<4
         ops(i)=appendSingletons(ops(i),'  -+');
      end
   end

   oez=oez(:,stype);
   for i=1:size(oez,1)
      for k=1:L
         if k>1, qdir='-+'; else qdir='++'; end
         oez(i,k)=appendSingletons(oez(i,k),['  ' qdir]);
         oez(i,k)=contract(oez(i,k),'12',EE(k),'12*');
         oez(i,k)=setitags(oez(i,k),['-A:H,H,S@' ndig],k);
      end
   end

   mm=zeros(nH,L);
   SS=QSpace(nH,L);

   zL=getIdentity(oez(1,  1),1,'-0');
   zR=getIdentity(oez(1,end),2,'-0');

   mL=-6;
   mR=-9;
   mz=0.5;

   for l=1:nH
      I=HS(l).iop; n=size(I,1); k1=I(1,1); k2=I(end,1); 

    % mm(l,:) = [ [-6 ... -6=mL], 2 0 1 0 2, [-9 ... -9=mR] ]
    % => for longer range terms
    %    intermediate sites with no operators acting have m=0
    % => sites to the left (right) of a Hamiltonian term have -10 (-20)
    % => intermediate fermionic `Jordan-Wigner' string has m=mz=0.77777
    %    which needs to be positive since m>0 indicates
    %    that SS(l,k) must be applied // Z_STRING
    % => mm>=0 indicates Hamiltonian term

      mm(l,1:k1-1)=mL; mm(l,k1:k2)=0; mm(l,k2+1:end)=mR;
      mm(l,I(:,1))=I(:,2); nz=0;

      for i=1:n, j=I(i,:); k=j(1);
         S=ops(j(2),j(3)); isf=isFerm(j(2),j(3));
         if isf && nz
            S=contract(HAM.oez(2,j(3)).op,2,S,1);
         end
         if n>1 && i==1
              S=contract(S,'21*',EE(k),'12*',[2 1 3]);
         else S=contract(S,'12', EE(k),'12*');
         end

         e=normQS(S); if e<1E-12
            wberr('got S projected out by EE(%g) (l=%g @ e=%.3g)',k,l,e);
         end

         if k==1, S=contract(zL,2,S,1); end

         SS(l,k)=setitags(S,['-A:H,H,S@' ndig],k);

         if isf, nz=nz+1; end

         if i<=1, continue; end

         kk=I(i-1,1)+1:I(i,1)-1; if isempty(kk)
            if nz==2 && isf, nz=0; end
            continue
         end

         Em=getIdentity(SS(l,I(i-1,1)),2);
         Em=appendSingletons(untag(Em),'  +-');

         if ~nz, iz=1;
         else    iz=2;
            if nz==2, nz=0; elseif nz~=1, wberr('invalid nz=%g',nz); end
            mm(l,kk)=mz;
         end

         for k=kk
            Q=contract(Em,'34',untag(oez(iz,k)),'12',[2 1 3]);
            SS(l,k)=setitags(Q,['-A:H,H,S@' ndig],k);
         end
      end

      if nz, wberr('invalid setup (got open z-string !?)'); end

      h=HS(l).hcpl;
      if n==1
         if addHc, h=(h/2); end
         SS(l,k)=SS(l,k)*h;
      else
         h=sqrt(abs(h))*[1, sign(h)];
         SS(l,k1)=SS(l,k1)*h(1);
         SS(l,k2)=SS(l,k2)*h(2);
      end
   end

   i=find(max(mm,[],1)<1);
   if ~isempty(i), q=sprintf(', %s',i);
      wberr('got non-referenced sites %s',q(3:end));
   end

   dk1=L;
   for i=find(mm(:,1)>=1)'
         dk1=min(dk1,find(mm(i,2:end)>=1,1));
   end
   if isempty(dk1), dk1=0; end

   dkL=0;
   for i=find(mm(:,L)>=1)'
        dkL=max(dkL,find(mm(i,1:end-1)>=1,1,'last'));
   end; dkL=L-dkL;
   if isempty(dkL), dkL=0; end

   if dk1>3, wblog('WRN',...
     'first site has range >= %g interactions',dk1+1); end
   if dkL>3, wblog('WRN',...
     'last site has range >= %g interactions',dkL+1); end

   Ham2=zeros(nH); fflag=(nH<75);

   if fflag
        fprintf(1,'   computing full |H|^2 ...  \r');
   else fprintf(1,'   computing diagonal contributions to |H|^2 ... \n');
   end
   tlast=now; dt=5/(24*3600);

   for ia= 1:nH, K1=HS(ia).iop(:,1)'; if fflag, iB=ia:nH; else iB=ia; end
   for ib=iB,    K2=HS(ib).iop(:,1)'; k1=min([K1,K2]); k2=max([K1,K2]);
       t=now; if t>tlast+dt, tlast=t;
       fprintf(1,'X  computing full |H|^2 (%.1f%%) ...  \r',100*ia/nH); end

       for k=k1:k2
          S2=[SS(ia,k), SS(ib,k)];
          for i=1:2, if isempty(S2(i)), S2(i)=oez(1,k); end; end
          Q=S2(2); if k>k1, Q={Xl,Q}; end
          Xl=contract(S2(1),'!2*',Q); if norm(Xl)<1E-8, break; end
       end
       Ham2(ia,ib)=getscalar(Xl);
   end
   end
   if fflag, fprintf(1,'%60s\r',''); end

   ham2_=diag(Ham2); 
   ham2=2*real(sum(Ham2(:)))-sum(ham2_);

   q=sort(sqrt(abs(ham2_))); k=0;
   while k<3, k=k+1;
      if     k==1, i=find(q<1E-12); s='very small';
      elseif k==2, i=find(q<1E-5);  s='small';
      else         i=find(q>1E5);   s='large'; end

      n=numel(i); if n, s={s,''};
         if n>1, if n>2, s{2}=' ..'; else s{2}=','; end
              s{2}=sprintf('%.3g%s %.3g',q(i(1)), s{2}, q(i(end)));
         else s{2}=sprintf('%.4g',q(i)); end

         wblog('WRN','got %g/%g %s HAM terms (%s)',n,nH,s{:});
         if k==1, k=k+1; end
      end
   end

   SA=QSpace(1,L);  % cummulative sum for either L or R block (`sum')
   SM=QSpace(nH,L);

   EM=QSpace(1,L);

   init_ops=1; Dmax=[]; NKEEP=[];
   if isb || vflag, nl={'',''}; else nl={'   \r\\',' ...'}; end

   mpo=oez(1,:);

   for isw=0:nsw, if init_ops, sdir_=-1; else sdir_=[+1,-1]; end
      for sdir=sdir_, if sdir>0, kk=1:L-1; else kk=L:-1:2; end
         for k=kk, q1=max(0,k-1); q2=k+1; if q2>L, q2=0; end
            if sdir>0
                 k1=k; k2=k+1; kp=q1; kn=q2; odir='>>'; ic='!2*';
            else k1=k-1; k2=k; kp=q2; kn=q1; odir='<<'; ic='!1*';
            end

            if vflag || k1>1 || sdir>0 || ~isw
                 q=getDimQS(mpo(k1)); q=q(1,2);
            else q=max([Isvd.Nkeep]); end

            k1_=max(0,k1-1);
            k2_=k2+1; if k2_>L, k2_=0; end

            ssw=sweep2str(isw,k,L,sdir);
            if sdir>0 && k2_ || sdir<0 && k1_
                 if init_ops, nl{3}='(init)'; else nl{3}=nl{2}; end
            else nl{3}='   '; end

            wblog('MPO',[ssw ' %3g %16.12g  %s' nl{1}], ...
               q,normQS(mpo(k))^2, iff(init_ops,'(init)',nl{3}));

          % ---------------------------------------------------------- %
            if ~init_ops
               E2=getIdentity(EE(1,k1),3,EE(1,k2),3);
                  t3={ getitags(oez(1,k1),3), getitags(oez(1,k2),3) };
                  t3{3}=sprintf('S%04g',(k1+k2)/2);
               E2=setitags(E2,t3); X=QSpace;

               m2=mm(:,[k1 k2]); nL=0; nR=0;

               m4=repmat(-99,nH,2); % NN outside => `4'-site range
               if k1_, m4(:,1)=mm(:,k1_); end
               if k2_, m4(:,2)=mm(:,k2_); end

               for l=1:nH, q2=m2(l,:);
                  if all(q2<0)
                     if     all(q2==mL), nL=nL+1; if nL>1, continue; end
                     elseif all(q2==mR), nR=nR+1; if nR>1, continue; end
                     else wberr('invalid q2 data [%g,%g] (l=%g)',q2,l); end
                  end

                  if   all(q2>=0), Q={  SS(l,k1),  SS(l,k2) };
                  elseif q2(1)>=0, Q={  SS(l,k1), oez(1,k2) };
                  elseif q2(2)>=0, Q={ oez(1,k1),  SS(l,k2) };
                  else             Q={ oez(1,k1), oez(1,k2) }; end
                  S2=contractQS(Q,E2,'*');

                  if all(q2<0), a=[]; p=[];
                     if q2(1)==mR
                        Q=S2; if k2_, Q={Q,EM(k2_),'*'}; p=[1 3 2]; end
                        a=contractQS(SA(k1_),'*',Q,p);
                     elseif q2(1)==mL
                        Q=S2; if k1_, Q={EM(k1_),'*',Q}; end
                        a=contractQS(Q,SA(k2_),'*',[1 3 2]);
                     else wberr('!?'); end

                     X=X+a; a=normQS(a);
                     if a<1E-6, s='  '; if isw<2, s=[s '\r\\']; end
                        wblog(' * ',['%s fixing |X(l=%2g)|=%g' s],ssw,l,a);

                        Q=S2; if k1_, Q={EM(k1_),'*',Q}; end
                              if k2_, Q={Q,EM(k2_),'*',[1 3 2]}; end
                        X=X+1E-3*contract(Q{:});
                     end
                  else
                     Q=S2; q4=m4(l,:)>=0; p=[];
                     if k1_
                        if q4(1), Xl=SM(l,k1_); else Xl=EM(k1_); end
                        Q={Xl,'*',Q};
                     end
                     if k2_
                        if q4(2), Xr=SM(l,k2_); else Xr=EM(k2_); end
                        Q={Q,Xr,'*'}; p=[1 3 2];
                     end
                     a=contractQS(Q{:},p);
                     X=X+a;
                  end
               end
               if isempty(X), wberr('got empty MPO !?'); end

               X4=contract(X,E2);

               if sdir>0
                  t=getitags(oez(1,k),2); t=['K' t(2:end)];

                  A2=getIdentity(X4,2,X4,4,t);
                  A1=contract(A2,'*',X4);
                  [u,xv_,Isvd(k1)]=orthoQS(A1,1,'<<');
                  mpo(k1)=permuteQS(u,[2 1 3]);
                  mpo(k2)=contractQS(xv_,2,A2,3);
               else
                  t=getitags(oez(1,k),1); t=['K' t(2:end)];

                  A1=getIdentity(X4,1,X4,3,t);
                  A2=contract(A1,'*',X4);
                  [u,xv_,Isvd(k1)]=orthoQS(A2,1,'<<');
                  mpo(k2)=u;
                  mpo(k1)=contractQS(xv_,2,A1,3,[2 1 3]);
               end
            elseif isw
               [mpo(k1),mpo(k2),r2,Io]=ortho2site(mpo(k1),mpo(k2),odir);
               Isvd(k1)=rmfield(Io,{'DX','Dk','dir'});
            end

          % ---------------------------------------------------------- %
            if kp
               Q={ SA(kp), mpo(k) }; if sdir<0, Q=Q([2 1]); end
               SA(k)=contract(oez(1,k),ic,Q);

               r=numel(SA(k).Q);
               if r>4, wberr('got growing tensor rank (r=%g)',r); end
            else
               SA(k)=QSpace;
            end

            for l=find(mm(:,k)>=0)', Q=mpo(k);
               if kp
                  if mm(l,kp)>=0, Q={SM(l,kp),Q}; else Q={EM(kp),Q}; end
                  if sdir<0, Q=Q([2 1]); end
               end

               SM(l,k)=contractQS(SS(l,k),ic,Q);

               if kn && mm(l,kn)<0
                  SA(k)=SA(k)+SM(l,k);
               end
            end

            Q=mpo(k); if kp, Q={ EM(kp), Q }; if sdir<0, Q=Q([2 1]); end; end
            EM(k)=contractQS(oez(1,k),ic,Q);

            x=norm(EM(k));
            if x<1E-8 && (sdir<0 && k>1+dk1 || sdir>0 && k<L-dkL)
             % NB! may have <MPO|E>=0 globally if Hamiltonian terms
             % have no overlap with identity matrix, e.g. tight-binding
             % at half-filling in which case turning on epsd!=0, indeed,
             % prevents this case to happen!
             % => assuming NN term between first or last pair of sites
             % the MPO must at have overlap with E at least for N-2 sites!
             % Wb,Nov28,20
               if isbatch, wberr('%s |EM|=%.4g',ssw,x);
               else wblog('WRN','%s |EM|=%.4g',ssw,x); wbstop
               end
            end
         end

         fprintf(1,'\n'); % '\r680s\r','');
         if isw, NKEEP(:,end+1)=[Isvd.Nkeep]; end
         init_ops=0;

      end

      hmpo=zeros(1,nH); EL=EM; kc=max(1,k-1);
      for k=1:L
          Q=mpo(k); if k>1, Q={EL(k-1),Q}; end
          EL(k)=contractQS(oez(1,k),'!2*',Q);
      end

      for l=1:nH, kk=HS(l).iop(:,1); k1=kk(1); k2=kk(end);
         for k=k1:k2
            Q=mpo(k); if k>k1, Q={Xl,Q}; elseif k>1, Q={EL(k-1),Q}; end
            if k==k2 && k<L, Q={Q,EM(k+1)}; end
            Xl=contract(SS(l,k),'!2*',Q);
         end
         hmpo(l)=getscalar(Xl);
      end

      if isw<1, continue; end
      Dmax(isw)=max([Isvd.Nkeep]);

    % add all longer range terms (NNN or further) within sweep isw=1
    % NB! if longer range terms e.g. contain Q=0 scalar operator
    % (e.g. as mixed operators together with S.S terms) one may have
    % n=0 while not really having the long-range terms fully included!
    % hence always inlcude *all* longer range terms with isw=1 (!)
    % even if terms may duplicated, the mpo is compressed later anyway!
      ll=find(HR(:,3)>1)';
      i=find(abs(hmpo)<1E-12); n=numel(i); q=norm(hmpo(i));

      if n
         if    n==1, s=sprintf('term %g',i);
         elseif n<4, s=['terms' sprintf(' %g',i)];
         else s=sprintf('%g terms [%g .. %g]',numel(i),min(i),max(i)); end

         s={' * ', s}; if isw>1 || isempty(ll), s{1}='WRN'; end
         wblog(s{1},'%s not yet in mpo @ %.2g',s{2},q);

      elseif isw>1, wblog(' * ',...
         'all terms contained in mpo (max hmpo @ %.3g)',min(abs(hmpo)));
         continue
      end
      if isempty(ll) || isw>=max(3,nsw), continue; end

      q=HR(ll,3); q=[min(q),max(q)];
      if diff(q), q=sprintf('%g..%g',q); else q=num2str(q(1)); end
      wblog('==>','%s adding %g/%g terms (range dk=%s)',ssw,numel(ll),nH,q);

    % NB! need to explicitly add long-range terms, as they will be
    % missed (projected out) with 2-site updates. Longer range terms
    % may also act non-trivially on arbitrarily many sites in between
    % such that simple fixes like temporary NN-mpo, etc.
    % are therefore also useless! // HAM_XTERMS

      [K1,is]=sort(getdatafield(HS(ll),'iop',1));
      ll=ll(is); kc=1; nx=numel(ll);

      for j=1:numel(ll)
         l=ll(j); k=HS(l).iop(:,1); k1=k(1); k2=k(end); clear Q_
         fprintf(1,'\r %8g / %g  [ %3g .. %3g ] @ %2g ... \r',...
            j,nx,k1,k2,k2-k1);

         if kc<K1(j) 
            for k=kc:K1(j)-1
               [mpo(k),mpo(k+1)]=ortho2site(mpo(k),mpo(k+1),'>>'); 

               Q=mpo(k); if k>1, Q={ EM(k-1), Q }; end
               EM(k)=contractQS(oez(1,k),'!2*',Q);
            end
            kc=K1(j);
         end

         for k=k1:k2, Q=SS(l,k); if mm(l,k)>0, Q_=Q; end
          % NB! using projection |M><M|E> of E (all identities) onto M (mpo)
          % where Q=<M|E> then connects to the remainder of |M> via oplus
          % => E and M must be orthonormalized
          % => by construction, the multiplet bond states in M must contain E
          %    hence EM must be normalized to 1. for k<kc and k>kc
          %    with the current side (kc) excempt.
          % tags: ORDER_START_STOP_NORM
            if k==k1
               if k>1, Q=contract(EM(k-1),'1*',Q,1); end
               mpo(k)=oplus(mpo(k),Q,'first');
            elseif k==k2
               if k<L, Q=contract(EM(k+1),'1*',Q,2,[2 1 3]); end
               mpo(k)=oplus(mpo(k),Q,'last');
            elseif mm(l,k)>0
               mpo(k)=oplus(mpo(k),Q);
            else
               Q=appendSingletons(getIdentity(Q_,2),'  +-');
               Q(2)=oez(1,k); Q=untag(Q);
               Q=contract(Q(1),'34',Q(2),'12',[2 1 3]);
               Q=setitags(Q,['-A:H,H,S@' ndig],k);
               mpo(k)=oplus(mpo(k),Q);
            end
         end

         for k=k2:-1:k1, Q=mpo(k);
            if k>k1, ic='!1*'; if k<L, Q={ Q, EM(k+1) }; end
            else     ic='!2*'; if k>1, Q={ EM(k-1), Q }; end; end
            EM(k)=contractQS(oez(1,k),ic,Q);
            if numel(EM(k).Q)>2, wberr('ERR',''); end
         end

      end

      for k=kc:L-1, [mpo(k),mpo(k+1)]=ortho2site(mpo(k),mpo(k+1),'>>'); end

      init_ops=1;

   end

   mpo_0=mpo;

   for k=1:L % tr(H*H')
      Q=mpo(k); if k>1, Q={x2,Q}; end
      x2=contract(mpo(k),'!2*',Q);
   end
   mpo2=getscalar(x2);

   for k=1:L
      u=contract(EE(k),'12',EE(k),'21');
      Q=contract(mpo(k),3,u,2);
      if k>1, Q=contract(h2,'2', Q,1);
      else    Q=contract(zL,'2*',Q,1); end
      h2=contract(mpo(k),'13',Q,'13');
   end
   h2=getscalar(h2);

   q=[mpo2, h2]; e=abs(diff(q))/norm(q);
   if e<1E-9
        wblog('ok.','got hermitian HAM @ e=%.3g',e); 
   else wblog('NB!','got non-hermitian HAM (@ %.6g)',e); end

   hmpo_=hmpo; hmpo=sum(hmpo);

   q=[ mpo2, real(hmpo) ]; if fflag, q(3)=ham2; end
   e=norm(diff(q))/q(1);

   s=sprintf('|HAM|^2 = %.8g (D=%g) @ %.3g',mpo2,max(Dmax),e); 

   if e<1E-8
      s={'-->',['got complete ' s]}; if fflag, s{1}='==>'; end
      wblog(s{:});
   else
      wblog('WRN','got difference %s',s);      if fflag
      wblog(' * ','  |HAM|^2 = %12.6g',q(3));  end
      wblog(' * ','  |mpo|^2 = %12.6g',q(1)); 
      wblog(' * ','<HAM|mpo> = %12.6g',q(2)); 
      kflag=kflag+1;
   end

   isferm=any(isferm);
   Nkeep=NKEEP(:,end)'; hmpo=hmpo_; nrm2=q;

   svd2tr=[Isvd.svd2tr];

   for k=numel(Isvd):-1:1, q=Isvd(k); s=q.svd;
      if size(s,2)>1, s2=s(:,2)'*s(:,1).^2; else s2=sum(s(:,1).^2); end
      svd_range(k,:)=[q.svd([1,q.Nkeep],1); s2];
   end

   if isfield(HAM.info,'mpo')
        Iout=HAM.info.mpo; if isempty(Iout), Iout=struct; end
   else Iout=struct; end

   Iout=structmerge(Iout,IH);
   Iout=add2struct(Iout,L,dloc,mm,'use_hconj:addHc',isferm,'EE:EE_',...
     Nkeep,svd_range,svd2tr, Ham2,hmpo,mpo2,nrm2,e);

 % Identify global energy offset in HAM propto Id
 % since this generates ambiguity with which site is going to
 % be associated with!
 % (1) obtain e0 for entire mpo (= trace(mpo) !)
 % (2) remove it from mpo while making it canonical, before finally,
 % (3) again adding it to the MPO uniformly over all sites

   if k~=1, wberr('got k=%d !?',k); end
   EM(k)=contract(oez(1,k),'!1*',{mpo(k),EM(k+1)}); EM_0=EM;

   xR=getscalar(EM(k));

   e0=xR;
   e1=e0/L;

   gotId=(abs(e1)>1E-12);
   if gotId
    % with Id for local identity and dfac=1/sqrt(dloc) above
    % and both, the mpo scaled by dfac as well as Id
    % tr[ ( Id/sqrt(dloc) )^2 ] = 1 => got term without any extra factor
    % e.g. wsys='tb3' (spinless fermion): epsd*n_i for all sites,
    % then with n_i = 1/2*(Id + sigma_z))
    % => total contribution to Id = L*(epsd/2)
    % => furthermore, if ~addHc, then half of this => (L/2)*(epsd/2) - ok!
    % e0 can be complex if ~addHc => use num2str()
    % => e1 = uniform contribution per site

      q=1; if addHc, q=2; end % NB! if addHc: mpo+mpo' => 2*e1
      wblog('-->','got Id term in HAM.mpo (e0=%s)',num2str(q*e1,'%.4g')); 
   end

 % Canonicalize mpo by inserting 1=U*U' along bonds // MPO_CANONICAL
 % NB! in order to evaluate block matrix elements entering <A|mpo|A>
 % it is useful to be able to split off, and thus to identify
 % the identity operator which corresponds to completed terms
 % need to bring (0,0,0) sector in mpo to canonical form:
 %     1st state: `start' state from the left // ORDER_START_STOP_NORM
 %     2nd state: `stop' state towards the right
 % Then
 %  => the cumulative matrix elements of mpo(k).data{l}(1,1)
 %     with l the (0,0*,0) sector describes the weight for
 %     the identity operator in the left block
 %  => the cumulative matrix elements of mpo(k).data{l}(2,2) (= 1!)
 %     with l the (0,0*,0) sector describes the weight for
 %     the identity operator in the right block; by the RL
 %     orthonormalization of the mpo, this coefficient is 1.
 % also check sign to make sign in EM(k) positive

   for k=L:-1:2
    % build weighted local Hermitian super-operator
    % to respect local operator order, as well as incoming ops;
    % => expect a more transparent structure of the mpo;
    % the unitary U below is irrelevant, otherwise;
    % what is important only, is the structure of the q=0 sector
    % which will be explicitly set (overwritten) anyways
      A=getIdentity(mpo(k),3); l=0;
      [~,is]=sortrows(abs(A.Q{1}));
      for i=is', s=size(A.data{i},1);
         A.data{i}=diag(l+1:l+s); l=l+s;
      end
      X=contract(mpo(k),'!1*',{mpo(k),A});

      A=getIdentity(mpo(k),2); l=0.5;
      [~,is]=sortrows(abs(A.Q{1})); % start with `smallest' q-labels
      for i=is', s=size(A.data{i},1);
         A.data{i}=diag(l+1:l+s); l=l+s;
      end
      X=X+contract(mpo(k),'!1*',{mpo(k),A});

      [ee,Ie]=eigQS(X); ee=ee(:,1);
      if any(ee<1), ee, wbdie('got min(ee)=%.3g [k=%g]',min(ee),k); end

      % NB! EM represents the identity applied all sites to the right
      % which is a well-defined single(!) super-state for the R-block!
      % Now Id in the R-block can be linked in at most two ways
      % to the left (L-block), such that EM(k) is row-vector:
      % (1) to HAM terms (H-terms) in the L-block that are completed
      %     which must come with Id for the R-block, or
      % (2) to the identity operator on the entire L-block itself
      %     e.g. in case of gotId where HAM has finite trace,
      %     i.e. has a trivial energy offset term propto Id globally.
      %
      %     Since the mpo can be tuned continuously w.r.t bare
      %     parameters in the Hamiltonian (1) must act independently
      %     from (2); therefore these must correpsond to independent
      %     `states' (except for trivial trace-full part of H-terms
      %     which actually represent the Id operator).
      %
      % Note, however, that the Q=0 bond state space must have
      % two states even if the Hamiltonian is traceless (~gotId)
      % because the distinction between completed terms and terms
      % to come (not yet started) must exist in any case [similar
      % to the previous argument, the mpo is continuous in bare
      % paraemeters of HAM, and so it must distinguish between
      % terms completed and terms still to come].
      % For example: wsys='tb3'; use_mpo=1; epsd=0; has two states
      % within the bulk of EM objects to start with, with significant
      % SVD weight on both Q=0 states so this is not by accident,
      % except for k<=2 (first bond 1-2), and k=L (last bond).
      %
      % Because of (1), EM must be non-zero as soon as H-terms
      % have been completed in the L-block (e.g. k > ~ 2 .. L)
      %
      % => (1) can be absent for small k only; save to assume that
      %    this cannot be the case k>L/2 for sufficiently large systems.
      % => (2) can be absent if there is no global Id term
      %    i.e. tr(H)=0; this must hold globally!

      q=norm(EM(k));
      if q<1E-12
       % NB! EM(k) can only be zero (empty) if both (1) and (2)
       % above are absent: this can only occur for small k if there are
       % no H-terms yet completed for k'<k (the L-block); e.g. this
       % may occur for k=2, which describes the bond 1-2).
       % Completed terms require Id for the R-block [see (1) above]!
       % [todo] Add infinitesimal local term to all sites so that
       % identity operator is not fully absent? // Wb,Feb15,21
         q={' * ', q }; if k>2, q{1}='WRN'; end
         if ~EM(k)
            wblog(q{1},'<%g| got empty EM',k); break
         end
      end

      if numel(EM(k).data)~=1 || size(EM(k).data{1},1)~=1
         wberr('unexpected EM(%g) dimensions',k);
      end

      v=EM(k).data{1}.'; v2=v*v'/(v'*v); p=eye(size(v2))-v2;
        z=zeros(1,size(EM(k).Q{1},2));
        i=matchIndex(X.Q{1},z);
      X.data{i} = p*X.data{i}*p + 0.5*v2;

      [ee,Ie]=eigQS(X);
      U=Ie.AK; i=matchIndex(U.Q{1},z); u=U.data{i};

      IE(k)=add2struct('-',X,U,ee,u,v);

      if v'*u(:,1)<-1E-10, u=-u; end

      if numel(u)>1
         u=u(:,[2 1,3:end]);
      elseif numel(u)~=1, u, wberr('unexpected setting'); 
       % NB! numel(u)==1 may only occur close to the boundary where
       % no H-terms are yet completed close to the L (or R) boundary
       % and where there is no global energy off set (~gotId)
       % => should be interpreted w.r.t. L-block where either the
       %    start or the stop state is present, but no longer both
       % => artificially inflate bond dimension // ORDER_START_STOP_NORM
      elseif k>L/2
         u=[0, u]; % expect u to represent `stop' state close to R-end
      else
         u=[u, 0]; % expect u to link to `start' state close to L-end
      end

      U.data{i}=u;

      mpo(k)=contractQS(U,'1',mpo(k),1);
      mpo(k-1)=contractQS(mpo(k-1),U,'*',[1 3 2]);
   end

   for k=L:-1:1
      Q=mpo(k); if k<L, Q={Q,EM(k+1)}; end
      EM(k)=contractQS(oez(1,k),'!1*',Q);
   end

 % NB! normalize EL: so far mpo is RL orthonormalized
 % => `stop' state is well-defined as second state within (0,0,0)
 % normalizing EL here, the mpo will no longer be strictly
 % RL orthonormal, at the benefit that Id_Left (like Id_Right)
 % is well-defined and normalized // ORDER_START_STOP_NORM
   EL=QSpace(1,L); xp=zeros(1,L);
   for k=1:L
      Q=mpo(k); if k>1, Q={EL(k-1),Q}; end
      EL(k)=contractQS(oez(1,k),'!2*',Q);

      if EL(k)
         q=EL(k).data{1}; a=abs(q); i=find(a==max(a),1);
         l=length(q); p=1:l;
         if i>2
            p([1 i])=[i 1];
            e=abs(q(2))>1E-8; if e || tflag
               s={' * ',sprintf('|%g> swaping i=(1,%g) / %g',k,i,l)};
               if e, s{1}='WRN'; end
               wblog(s{:}); displ(q')
            end
         elseif i==2
            if k>1 && k+dkL<L, displ(q');
               wblog('WRN','|%g> got max xL at i=%g !?',k,i); end
            i=1; 
         end
         q=q(i);
      else q=0; end

      if abs(q)>1E-6 && k<L
         a=struct(mpo(k));
         for i=matchIndex(a.Q{2},z)
            x=a.data{i}(:,p,:); x(:,1,:)=(1/q)*x(:,1,:);
            mpo(k).data{i}=x;
         end

         a=struct(mpo(k+1));
         for i=matchIndex(a.Q{1},z);
            x=a.data{i}(p,:,:); x(1,:,:)=(  q)*x(1,:,:);
            mpo(k+1).data{i}=x;
         end

      elseif k+2<L
         wblog('WRN','|%g> got EL->(1) = %.3g',k,q);
      end

      Q=mpo(k); if k>1, Q={EL(k-1),Q}; end
      EL(k)=contractQS(oez(1,k),'!2*',Q);
   end

   irep=0; if gotId, irep=1; end
   for irep=0:irep

      if irep, xL_=xL; EM_=EM; xR_=xR; EL_=EL;
      for k=1:L
         Q=mpo(k); if k>1, Q={EL(k-1),Q}; end
         EL(k)=contractQS(oez(1,k),'!2*',Q);
      end, end
      q=struct(EL); q=[q.data]; xL=cat2(1,q{:});
      xL(end,:)=flip(xL(end,:));

      for k=L:-1:1
         Q=mpo(k); if k<L, Q={Q,EM(k+1)}; end
         EM(k)=contractQS(oez(1,k),'!1*',Q);
      end
      q=struct(EM); q=[q.data]; xR=cat2(1,q{:}); 

      q=[xR(1); diag(xL(1:end-1,:)*xR(2:end,:).'); xL(end); xR(1); xL(end)];
      e=[ norm(q-e0), abs(e0) ]; e=e(1)/max(1,e(2));
      if any(e>1E-12), s=sprintf(', %.3g',e);
         wberr('inconsistent Id term (e=%s)',s(3:end));
      end

      if ~irep
         if gotId, e1=repmat(e0/L,1,L); else e1=[]; end

         i=find(abs(xL(:,1))>1E-12,1,'last'); e=norm(xL(1:i,1)-1);
         if e>1E-10
            s=sprintf('failed to canonicalize start state (EL @ e=%.3g)',e);
            if Fflag, wblog('ERR',s); else wberr(s); end
         elseif i<L, k=i+1;
            wblog(' * ','got start states for k<L-%g @ e=%.1g',L-k,e);
            if k<L && gotId, e1(k)=sum(e1(k:L)); e1(k+1:L)=0; end
         end

         i=find(abs(xR(:,2))>1E-12,1); e=norm(xR(i:end,2)-1);
         if e>1E-10
            s=sprintf('failed to canonicalize stop state (ER @ e=%.3g)',e);
            if Fflag, wblog('ERR',s); else wberr(s); end
         elseif i>1, k=i-1;
            wblog(' * ','got stop  states for k>%-3g @ e=%.1g',k,e);
            if k>1 && gotId, e1(k)=sum(e1(1:k)); e1(1:k-1)=0; end
         end

         continue
      end

      check_mpo_overlap(mpo_0,mpo, mpo2); mpo_1=mpo;

      for k=1:L, if abs(e1(k))
         if k>1, a=EL(k-1).data; else a={1}; end
         if k<L, l=k+1;
            Q=mpo(l); if l<L, Q={Q,EM(l+1)}; end
            EM(l)=contractQS(oez(1,l),'!1*',Q);
            b=EM(l).data;
         else b={1}; end

         P=getIdProj(mpo(k),a,b,oez(1,k));
         if P
            xp(k)=getscalar(contract(P,'*',mpo(k)));
            mpo(k)=mpo(k)+(e1(k)-xp(k))*P;
         end
      end; end
   end

   check_mpo_overlap(mpo_0,mpo, mpo2,'-v');

   xR=fliplr(xR);

   HAM.info.mpo=add2struct(Iout,xL,xR,e0,e1);

   dfac=1/dfac;
   for k=1:L, tt=getitags(mpo(k));
    % NB! think about the contraction below as mpo(k)*[dfac*EE(k)]
    % rather than  [dfac*mpo(k)]*EE(k); here EE => dfac*EE restores
    % the original normalization of the operators; e.g. for matrix
    % elements which are say 1. in mpo(k) may still be 1. after
    % contracting (dfac*EE)
      Q=contract(dfac*mpo(k),3,EE(k),3);

      tt(1:2)=regexprep(tt(1:2),'^[A-Z]','m'); % 'mpo' 'h' 'H'
      tt(3:4)={ regexprep(tt{3},'^S*(\d+).*','s$1') };
      mpo(k)=setitags(Q,tt);
   end
   if 0 || kflag, wbstop; end

   mpo(1)=squeeze(mpo(1),1);
   mpo(L)=squeeze(mpo(L),2);

   HAM.mpo=mpo;

end

% -------------------------------------------------------------------- %

function s=sweep2str(isw,k,L,sdir)

   if L<=99
        s=sprintf('%g.%02g',isw,k);
   else s=sprintf('%g.%03g',isw,k); end

   if sdir>0, s=['|' s '>']; else s=['<' s '|']; end

end

% -------------------------------------------------------------------- %
function P=getIdProj(mpo_k,El,Er,oez_k)

   n=[ numel(El), numel(Er), numel(oez_k.data) ];
   if any(n>1), n, wberr('invalid usage'); end
   if any(~n), P=[]; return; end

   q=[El, Er, oez_k.data];
   for k=1:2
      n=numel(q{k});
      if n~=2
         if n~=1, wblog('WRN','got %s of size %s',s{2},sizestr_(q{k})); end
         continue
      end

      if k==1
           i=1; s={'start','EL'};
      else i=2; s={'stop', 'ER'}; end

      if abs(q{k}(i)-1)<1E-12, j=1:n; j(i)=[]; q{k}(i)=1; q{k}(j)=0;
      else wblog('WRN',...
        'invalid %s vector with %s (%s)',s{:},num2str(q{k}(i),'%.4g'));
      end
   end

   P=setitags(oez_k,getitags(mpo_k));

   s=zeros(1,3);
   for i=1:3
      q{i}=reshape(q{i},[],1); s(i)=numel(q{i});
      nrm=norm(q{i}); if nrm<1E-12, nrm=-nrm; break; end
      q{i}=q{i}/nrm;
   end

   if nrm>0
      P.data{1}=reshape(kron(q{3},kron(q{2},q{1})),s);
   end
end

% -------------------------------------------------------------------- %
function [x,D]=mpo_overlap_0(SS,mpo,oez)

   L=length(mpo); nH=size(SS,1); x=[];

   if nH>10, fprintf(1,...
   '  computing mpo overlap with %g SS terms ...  \r',nH); end

   XL=QSpace(1,L);
   for k=1:L
       Q=mpo(k); if k>1, Q={XL(k-1),Q}; end
       XL(k)=contractQS(oez(1,k),'!2*',Q);
   end
   XR=QSpace(1,L);
   for k=L:-1:1
       Q=mpo(k); if k<L, Q={Q,XR(k+1)}; end
       XR(k)=contractQS(Q,oez(1,k),'!1*');
   end

   SS=struct(SS); mm=zeros(size(SS)); K2=zeros(nH,2);
   for i=1:numel(SS)
      if ~isempty(SS(i).Q), mm(i)=1; end
   end
   for l=1:nH
      k=find(mm(l,:));
      if isempty(~k) || ~isequal(k,k(1):k(end))
         wberr('unexpected S.S data'); end
      K2(l,:)=k([1,end]);
   end

   for l=1:nH, k1=K2(l,1); k2=K2(l,2);
      for k=k1:k2, Q=mpo(k);
         if k>k1, Q={Xl,Q}; elseif k>1, Q={XL(k-1),Q}; end
         if k==k2 && k2<L,  Q={Q,XR(k+1)}; end
         Xl=contractQS(SS(l,k),'!2*',Q);
      end
      x(l)=getscalar(QSpace(Xl));
   end

   if nargout>1, D=[];
      for k=L:-1:1
         q=getDimQS(mpo(k));
         D(k,:)=q(:,2);
      end
   end
end

% -------------------------------------------------------------------- %
function [x2,e]=check_mpo_overlap(A,B,x2ref,varargin)

   getopt('init',varargin);
      vflag=getopt('-v');
   getopt('check_error');

   if ~isequal(size(A),size(B))
      wberr('invalid usage (size mismatch)'); end
   L=numel(A);

   for k=1:L
      Q=B(k); if k>1, Q={Xl,Q}; end
      Xl=contractQS(A(k),'!2*',Q);
   end

   x2=getscalar(QSpace(Xl));

   if nargin<3 || isempty(x2ref), return; end

   if numel(x2ref)~=1 || ~isfinite(x2ref) || ~abs(x2ref)
      x2ref, wberr('invalid x2ref');
   end

   e=abs(x2-x2ref)/max(1,norm(x2ref));
   if nargout<2
      if e>1E-12, wberr('mpo changed @ %.3g',e);
      elseif vflag, wblog('MPO','consistent @ e=%.3g',e);
      end
      if ~nargout, clear x2; end
   end
end

% -------------------------------------------------------------------- %

