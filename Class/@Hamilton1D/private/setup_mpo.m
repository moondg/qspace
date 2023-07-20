function HAM=setup_mpo(HAM,HH,stype,heps)
% function HAM=setup_mpo(HAM,HH [,stype,heps])
%
%    Automated routine to generate pseudo MPO (HAM.mpo) for the
%    Hamiltonian described by HH together with HAM.ops and stype.
%
%    The pseudo MPO generated here is complete *except* for
%    the fermionic z-operator to be applied to closing coupling
%    operators, i.e. F'*(ZF); the latter Z is applied for
%    fermionic operators in updateHK(). // FERMIONIC_ZF
%    (in principle, ZF could be provided as an independent extra
%    operator in the ops set, which would make the MPO description
%    complete; however, this artificially blows up the number of
%    operators and hence the MPO in d_loc).
%
% Wb,Aug25,15

   if size(HH,2)~=5
    % 5 columns required:
    % [ site1, iop1, site2, iop2, coupling strength ]
      wbdie('invalid usage (HH requires 5 columns)');
   end

   if nargin<4, heps=1e-15; end
   i=find(abs(HH(:,end))<heps);
   if ~isempty(i), e=norm(HH(i,end));
      if e, wblog(' * ',...
        'eliminating %g near-zero terms in Hamiltonian (@ %.3g)',numel(i),e);
      else wblog(' * ','eliminating %g zero-terms in Hamiltonian',numel(i));
      end
      HH(i,:)=[];
   end

   H0=HH;

   i=unique(reshape(HH(:,[1 3]),[],1));
   L=i(end);
   if ~isequal(i,(1:L)')
      wbdie('invalid usage [HH is not addressing all sites !?]');
   end

   i=find(HH(:,1)>HH(:,3));
   HH(i,1:4)=HH(i,[3 4,1 2]);
   nH=size(HH,1);

   HU=uniquerows(HH(:,1:4));
   if size(HU,1)~=size(HH,1)
      wbdie('invalid HH (got non-unique records !?)');
   end

   check_ortho_ops(HAM.oez,HAM.ops,HH);

   nops=size(HAM.ops,1); ntype=size(HAM.ops,2);
   if size(HAM.oez,2)~=ntype, wbdie(['invalid usage '...
     '(same number of sets for ops and oez required: %d/%d)'],...
     size(HAM.oez,2),ntype);
   end

   if nargin<3 || isempty(stype)
      stype=ones(1,L);
      if ntype~=1, wblog('ERR',...
         'HAM.ops needs to be column vector (having stype=1 / %d)',ntype);
      end
   else
      if numel(stype)~=L, wbdie(...
        'length inconsistency (L=%g/%g) !?',L,numel(stype));
      end
      q=unique(stype(:));
      if any(q<1 | q>ntype) wbdie(['invalid usage '...
        '(stype out of bounds: %g..%g/%g)'],q([1 end]),ntype);
      end
      if ~isequal(q,(1:ntype)')
         wblog('WRN','unused local state space [%s; %g]',...
         vec2str(setdiff((1:ntype)',q),'-f'),ntype);
      end
   end

   qops=ones(size(HAM.ops)); n=numel(HAM.ops);
   for i=1:n
      if isempty(HAM.ops(i).op) || isempty(HAM.ops(i).op.Q)
      qops(i)=0; end
   end

   if isempty(HAM.oez), wbdie('invalid usage (HAM.oez not yet setup)');
   elseif size(HAM.oez,1)>2, wbdie(...
     'invalid usage (HAM.oez got %d operators !?',size(HAM.oez,1));
   end

   qoez=ones(size(HAM.oez)); n=numel(HAM.oez);
   for i=1:n
      if isempty(HAM.oez(i).op) || isempty(HAM.oez(i).op.Q)
      qoez(i)=0; end
   end

   if any(qoez(1,:)==0)
      wbdie('invalid usage (HAM.oez not yet setup !?)'); end
   if size(qoez,1)==1, qoez(2,:)=0; end

   for j=1:size(HAM.oez,2)
      if ~isIdentityQS(HAM.oez(1,j).op)
      wbdie('invalid oez (expecting E as first operator)'); end
   end

   q=uniquerows([ HH(:,3)-HH(:,1), HH(:,[2 4 5])]);
   n=size(q,1);
   if n==1
      wblog(' * ','Hamilton setup HH -> [dk,op1,op2,J]=[%g,%g,%g,%g]',q);
   elseif n<8
      wblog('mpo','Hamilton setup HH -> [distance, op1/op2, coupling]');
      s=sprintf('%4g %3g/%g %6g\n',q');
      wblog('   ',s(1:end-1));
   else
      wblog(' * ','HH got %g different entries',n);
   end

   il=find(HH(:,3)<=0 & HH(:,4)<=0);
   HH(il,3:4)=HH(il,1:2);

   il=find(HH(:,1)==HH(:,3)); i=find(HH(il,2)-HH(il,4));
   if ~isempty(i), HH(il(i),:), wbdie(...
      'got invalid local terms (inconsistent ops index)');
   end

   I=[ HH(il,2), reshape(stype(HH(il,1)),[],1) ];
   I=uniquerows(I);
   for i=1:size(I,1)
      X=HAM.ops(I(i,1),I(i,2));
      r=numel(X.op.Q); if r~=2
         wbdie('got non-scalar local Hamiltonian term'); end
      e=normQS(X.op-X.op'); if e>1E-8 && ~X.hconj
         wbdie('got non-hermitian local Hamiltonian term');
      end
   end

   HAM.info.HH=HH; 

   q=HH(:,[1 3 end]);
     i=find(q(:,1)==q(:,2)); if ~isempty(i), q(i,3)=q(i,3)/2; end
   H=sparse(q(:,1),q(:,2),q(:,3),L,L); H=H+H';

   IH.H0=H0;
   IH.H=H;

   [dx,~,d2]=uniquerows(abs(HH(:,1)-HH(:,3)));
   IH.ldist=[dx,d2'];
   IH.lmax=dx(end);

   h2=uniquerows([HH(:,[1 2]); HH(:,[3 4])]);
   if any(h2(:,2)<1 | h2(:,2)>nops)
      wbdie('invalid usage (operator index out of bounds)'); end

   [k,I1,D1]=uniquerows(h2(:,1));
   for l=1:numel(I1);
      i=unique(h2(I1{l},2)); j=stype(k(l));
      if any(qops(i,j)==0), wbdie( ...
        'invalid usage (addressing empty operator in HAM.ops !?)');
      end

      ff=[HAM.ops(i,j).fermionic];
      if any(ff)
         if ~qoez(2,j), wbdie(['invalid usage ' ...
            '(missing or mixed up fermionic parity operator Z)']);
         end
      else qoez(2,j)=0; end
   end

   if all(qoez(2,:)==0)
      if size(HAM.oez,1)>1
         wblog('mpo','unsetting z-ops since not required'); 
      end
      HAM.oez=HAM.oez(1,:);
   end
   noez=size(HAM.oez,1);

   for i=1:nH
      k1=HH(i,1); i1=HH(i,2); j1=stype(k1); o1=HAM.ops(i1,j1);
      k2=HH(i,3); i2=HH(i,4); j2=stype(k2); o2=HAM.ops(i2,j2);
      if k1~=k2
         [i1,i2,I]=matchIndex(o1.qop,o2.qop);
         if isempty(i1), wbdie(...
           'coupling term cannot be contracted to scalar !?'); end
         if o1.fermionic~=o2.fermionic, wbdie(...
           'coupling term got fermionic mismatch !?'); 
         end
      else
         e=[norm(o1.qop), norm(o2.qop)];
         if norm(e), wbdie('got non-scalar local operator !?'); end
         if o1.fermionic || o2.fermionic
            wbdie('local fermionic terms !?'); 
         end
      end
   end

   for l=1:2
      if l==1, j=[1 2, 3]; else j=[1, 3 4]; end
      [Q,I,D]=uniquerows(HH(:,j)); i=find(D>1);
      if ~isempty(i)
        if l==1, s='start'; else s='target'; end
        wblog('WRN','combine operators that link to same %s site',s);
        fprintf(1,'\n'); disp(HH([I{i}],:))
      end
   end

 % -------------------------------------------------------------------- %
 % Conventions:
 % * the first two states in M, i.e. M(1:2,1:2,:), always stand
 %   for the `start and stop states' |0> and |e>, respectively
 % * the operators along the 3rd dimension in M(:,:,k) are
 %   [ HAM.oez, HAM.ops ] in this order, that is:
 %
 %   -> k=1: always links to the local Id operator;
 %      this propagates the states |0> and |e> without altering them.
 %   -> k=2 (in the presence of fermionic operators only) links
 %      to the fermionic parity operator to take care of fermionic
 %      Jordan Wigner like strings
 %   -> remaining operators for k>3 (or 2) are linked to ops(k-noez)
 % -------------------------------------------------------------------- %

   M=zeros(2,2,noez+nops); M(:,:,1)=eye(2);

 % -------------------------------------------------------------------- %
 % NB! for propagated operators, M by itself looses the information
 % from which operator in HAM.ops the matrix elements were computed
 % => use extra index iop that points into HAM.ops(iop)
 % => iop has two columns (for 1st and 2nd index in M, respectively)
 %
 % iop(:,:,k={1,2}) indicates start (1) or end (2) point of coupling term
 % iop(:,j={1,2},:) indicates left (1) or right (2) bond index of M-tensor
 % => iop(i,1,:) directly matches to M(i,:,:)
 % => iop(j,2,:) directly matches to M(:,j,:)
 % -------------------------------------------------------------------- %

   mpo=setup_empty('L',L,'M',M); mpo=mpo.mpo;

   if any(stype~=1)
        for k=1:numel(mpo), mpo(k).stype=stype(k); end
   end; clear stype

   sfac=2^-10;
   if sfac*(size(HAM.oez,1)+size(HAM.ops,1))>=1
      wbdie('invalid usage (got %g+%g local operators !?)',...
      size(HAM.oez,1), size(HAM.ops,1))
   end

   ib=find(HH(:,1)==HH(:,3));
   for l=1:length(ib), i=ib(l);
      k=HH(i,1);
      io=noez+HH(i,2);
      mpo(k).M(1,2,io)=HH(i,end);

   end

   HM={HH(ib,:)};
   HH(ib,:)=[];

   ifac=2^-10; qfac=[1 ifac]'; qfac=blkdiag(qfac,qfac);
   if ifac*ntype>=1, wbdie(...
     'invalid usage (too many operator sets %g/%g !?)',ntype,1/ifac);
   end

   for dk=1:L, k=1; nfound=0;
      while 1
         if k>L-dk
            if nfound, k=1; nfound=0;
            else break; end
         end

         ik=find(HH(:,1)==k);

         [Q,I,D]=uniquerows(HH(ik,1:2));
         for i=1:numel(I), l=ik(I{i});
            D(i)=min(abs(HH(l,1)-HH(l,3)));
         end

         i=find(D==dk,1);

            if isempty(i), k=k+1; continue;
            else nfound=nfound+numel(I{i}); end

         l=ik(I{i}); [D,is]=sort(abs(HH(l,1)-HH(l,3)));
         ik=l(is); dd=[];

         nk=numel(ik); k1=HH(ik(1),1); is=mpo(k1).stype;

       % ----------------------------------------------------------- %
       % NB! the closing operator is assumed to match the opening
       % operator in terms of symmetries and that they are combined
       % with a dagger, i.e. S1'*S2. // Wb,Mar31,16
       % For fermionic operators, however, this requires F'*(ZF).
       % where the parity operator in the 2nd term cannot be taken
       % care of here! Instead, it will be added in updateHK()
       % *in addition* to the MPO specification setup here!
       % Wb,Aug26,15 // tags: FERMIONIC_ZF
       % ----------------------------------------------------------- %
         io1=HH(ik(1),2); i31=noez+io1;
         ko1=k1+io1*sfac;
         if HAM.ops(io1,is).fermionic, z=2; else z=1; end

         if diff(HH(ik,3))<0, HH(ik,:)
            wbdie('got sorting mismatch (%g..%g) !?',k1,k2); end

         for j=1:nk, i=ik(j); k2=HH(i,3);
            io2=HH(i,4); i32=noez+io2; ko2=k2+io2*sfac;
            k12=[ko1 ko2];

            for l=k1:k2, M=mpo(l).M; iop=mpo(l).iop; y=1:2;
               s=size(M); r=s(1); s=s(2);
               if l==k1
                   if l<k2
                      if j==1, M(1,s+1,i31)=1;
                      else     M(r,s+1,z)=1; end;                     iop(s+1,2,y)=k12;
                   else
                            if M(r,2,i32), wbdie('overwriting value in M !?'); end
                            M(r,2,i32)=HH(i,end);
                   end
               elseif l<k2, M(r+1,s+1,z)=1;         iop(r+1,1,y)=k12; iop(s+1,2,y)=k12;
               else         M(r+1,2,i32)=HH(i,end); iop(r+1,1,y)=k12;
               end
               mpo(l).M=M; mpo(l).iop=iop; dd=[dd; size(M)];
            end

            for l=k1:k2, M=mpo(l).M; iop=mpo(l).iop; y=1:2;
               if l>1 && size(mpo(l-1).M,2)~=size(mpo(l).M,1)
               wbdie('invalid mpo.M !?'); end
               if l<length(mpo) && size(mpo(l).M,2)~=size(mpo(l+1).M,1)
               wbdie('invalid mpo.M !?'); end
            end

            k1=k2;
         end

         HM{end+1}=HH(ik,:);

         HH(ik,:)=0;
         k=HH(ik(1),3);
      end
      HH(find(sum(HH.^2,2)==0),:)=[];
   end

   if ~isempty(HH), HH, wbdie('failed to include all H terms !?'); end

   stol=1E-12;
   [NK,SS,mpx]=check_ortho(mpo,stol);
   ss=cat(1,SS{:}); smin=min(ss); smax=max(ss);

   nk=max(NK,[],1);
   if diff(nk)
      wblog('WRN','got compressible mpo !?');
   else
      wblog(' * ','dimension of pseudo-mpo: D=%g',nk(1));
      wblog(' * ','having S=[ %.4g .. %.4g ], |mpo|^2 = %.8g',...
      smin,smax,norm(reshape(mpx(end).M,[],1))^2);
   end

   HAM.info.mpo=add2struct(IH,HM,NK,SS);

   mpr=reduce_RL(mpo);

   q=[ mpsoverlap({mpr.M},'-q'),  mpsoverlap({mpo.M},'-q'), ...
       mpsoverlap({mpx.M},'-q'),  mpsoverlap({mpr.M},{mpo.M},'-q'), ...
       norm(reshape(mpx(end).M,[],1))^2
   ];
   e=norm(diff(q))/max(abs(q));
   if e>1E-10, q
        wblog('WRN','got MPO norm inconsistency @ %.3g !?',e); 
   else wblog('mpo','MPO norm consistency @ %.3g',e); 
   end

   HAM.mpo=mpr;

end

% -------------------------------------------------------------------- %
function mpo=reduce_RL(mpo)

   L=numel(mpo); nx=zeros(L,1);

   for k=L:-1:2
      M=mpo(k).M; q=M(3:end,:,:);
      [q,I,D]=uniquerows(q(:,:));
      if any(D>1)
         X=eye(sum(D)); Ik=zeros(size(D));
         for i=1:length(D)
            l=I{i}(1); Ik(i)=l;
            X(:,l)=sum(X(:,I{i}),2); I{i}(1)=[];
         end

         Ix=[I{:}]; nx(k)=length(Ix); M_=M;
         X(:,Ix)=[]; X=blkdiag(eye(2),X);
         Ik=sort([1 2, Ik+2]); M=M(Ik,:,:);

         x=contract(X,M,2,1); e=norm(reshape(x-M_,[],1));
         if e>1E-12, wbdie('got e=%.4g !?',e); end

         mpo(k-1).M=contract(mpo(k-1).M,X,2,1,[1 3 2]);
         mpo(k).M=M;

         Ix=Ix+2; nk=numel(Ik);
         q=mpo(k  ).iop; q(1:nk,1,:) = q(Ik,1,:); q(nk+1:end,1,:)=0;
         mpo(k  ).iop=q;
         q=mpo(k-1).iop; q(1:nk,2,:) = q(Ik,2,:); q(nk+1:end,2,:)=0;
         mpo(k-1).iop=q;

         wblog(' * ','reduced mpo(%g) D=%g->%g',k,size(X));
      end
   end
end

% -------------------------------------------------------------------- %
function [NK,SS,mpo]=check_ortho(mpo,stol)

   L=numel(mpo);
   for k=L:-1:2
      NK(k-1,1)=size(mpo(k-1).M,2);
      [mpo(k-1).M,mpo(k).M,I(k)]=mpo_ortho(mpo(k-1).M,mpo(k).M,'<<',stol);
   end
   for k=2:L
      [mpo(k-1).M,mpo(k).M,I(k)]=mpo_ortho(mpo(k-1).M,mpo(k).M,'>>',stol);
      NK(k-1,2)=size(mpo(k-1).M,2);
   end

   SS={I(2:end).S};
end

% -------------------------------------------------------------------- %
function [M1,M2,I]=mpo_ortho(M1,M2,dir,stol)

   X=contract(M1,M2,2,1); s=size(X); s(end+1:4)=1;
   [U,S,V]=svd(reshape(X,prod(s(1:2)),prod(s(3:4))),'econ');

   S=diag(S); i=find(S>stol);
   U=U(:,i); V=V(:,i); S=S(i); I.S=S; S=diag(S);
   I.Nkeep=[length(i), length(I.S)];

   s=[s(1:2), length(i), s(3:4)];

   if isequal(dir,'>>')
      M1=reshape(U,   s(1:3  ));
      M2=reshape(S*V',s(3:end));
   elseif isequal(dir,'<<')
      M1=reshape(U*S, s(1:3  ));
      M2=reshape(V',  s(3:end));
   else dir
      wbdie('invalid direction');
   end

   M1=permute(M1,[1 3 2]);

end

% -------------------------------------------------------------------- %
function check_ortho_ops(oez,ops,HH)

  q=QSpace(size(oez)); for i=1:numel(q), q(i)=oez(i).op; end; oez=makeIrop(q);
  q=QSpace(size(ops)); for i=1:numel(q), q(i)=ops(i).op; end; ops=makeIrop(q);

  i1=unique(HH(:,2))'; n1=size(oez,1);
  i2=unique(HH(:,4))'; nt=size(ops,2); e=0;

  for k=1:nt
     X=[oez(:,min(k,end)); ops(:,k)]; n=numel(X); ee=zeros(n);
     for i=1:n
        for j=i+1:n
            ee(i,j)=normQS(contractQS(X(i),'12*',X(j),'12'));
        end
     end

     i=[1:n1,i1+n1]; q=norm(ee(i,i));
     i=[1:n1,i2+n1]; q(2)=norm(ee(i,i));

     e(k)=norm(ee);
     if any(q>1E-12)
        if nt>1, s=sprintf('(:,%d)',k); else s=''; end
        [i,j]=find(abs(ee)>1E-12); i=unique([i;j])-size(oez,1);
        wblog([ 'NB! [oez; ops%s] contains non-orthogonal ops ' ...
          '(%s @ e=%.3g)'], s, vec2str(i(i>0),'%d,'), e(k));
     end
  end

end

% -------------------------------------------------------------------- %

