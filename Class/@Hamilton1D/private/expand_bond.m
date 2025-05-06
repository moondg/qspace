function [X1,X2,Psi,Iout]=expand_bond(HAM,X1,X2,varargin)
% function [X1,X2,Psi]=expand_bond(HAM,X1,X2 [,J,hconj])
%
%    Increase / expand bond dimension based on H|psi>
%    hence within the 2-site setting sharing the bond in between.
%    This assumes a full MPO representation of HAM.
%    In case of a pseudo-MPO, this temporarily expands
%    to the full 2-site setting for H|Psi>.
%
%    This uses HAM.info.sweep.[Nkeep|facCBE]
%    where Nkeep is expected the final target dimension
%    after the DMRG update. Therefore the intermediate
%    bond dimension will be increased to Nkeep*facCBE
%    assuming facCBE >= 1. Assuming that the current
%    dimension is Nkeep, this extends the kept space (K)
%    by Nkeep*(facCBE-1) discarded states with dominan
%    weigth in H|Psi>.
%
%    Default: facCBE=max(1.10,1+8/Nkeep);
%    i.e., keep at least 8 additional states / multiplets
%
% Wb,Sep14,23

% adpated from get_HPsi.m and update_psi_2site.m

  f='--full-mpo'; wCBE=1;

  if got_gauge(HAM);
     wbdie('CBE not yet implemented for gauge fields');
  end

  nargs=numel(varargin); mark=zeros(1,nargs);
  for i=1:nargs, if ischar(varargin{i}), mark(i)=1; end; end

  getopt('init',varargin(find(mark)))
     nx=getopt('nx',8);
     fullMPO=getopt(f);
     if     getopt('--smin'), wCBE=2;
     elseif getopt('--smax'), wCBE=1;
     elseif getopt('--sall'), wCBE=0;
     end
  getopt('check_error');

  args=varargin(find(~mark)); nargs=numel(args);
  if fullMPO
     if nargs, disp(args);
        wbdie('invalid usage (additional args with %s)',f); end
     args={f};
  elseif nargs~=2
     wbdie('invalid usage (got %d / 2 args with pseudo MPO)',nargs);
  end

  k1_=get_kidx(X1.AK); k1=[ abs(k1_(1:2)), k1_(3) ]; n   =numel(k1_);
  k2_=get_kidx(X2.AK); k2=[ abs(k2_(1:2)), k2_(3) ]; n(2)=numel(k2_);
  if ~isequal(diff(k1),[1 0]) || ~isequal(k1+1,k2) || any(n>4)
     wbdie('unexpected input for 2-site bond expansion');
  end
  k1=k1(3);
  k2=k2(3);

  I=HAM.info.sweep;
  Nkeep=I.Nkeep;

  if isfield(I,'facCBE')
       facCBE=I.facCBE;
  else facCBE=1.10; end

  if I.facCBE<1
     wbdie('invalid facCBE=%g (expecting >=1)',I.facCBE); end
  facCBE=max(facCBE,1+nx/Nkeep);

  Ntot=round(facCBE*Nkeep);

  El=get_local_id(HAM,k1,X1.AK,3);
  Er=get_local_id(HAM,k2,X2.AK,3);

  [t,c,m]=getitags(X1.AK,2);
  t=regexprep(t,'^[a-zA-Z]+(\d+)','E$1','emptymatch');
  t={[t 'l'], [t 'r']}; if xor(m,c), i=1; else i=2; end
  t{i}=[t{i} ''''];

  E1=getIdentity(X1.AK,1,El,['-m:' t{1}],[1 3 2]);  d1=getDimQS(E1);
  E2=getIdentity(X2.AK,2,Er,['-m:' t{2}],[3 1 2]);  d2=getDimQS(E2);

  Dfac=sqrt(prod([ d1(end,2)/d1(1,2), d2(end,1)/d2(1,1) ]));

  d3=min([ d1(1,2), d2(1,1) ]);
  if d3<=Ntot
     if     numel(X1.AK.Q)>3, Eg=getIdentityQS(X1.AK,4);
     elseif numel(X2.AK.Q)>3, Eg=getIdentityQS(X2.AK,4); else Eg=QSpace; end

     if Eg, q=getDimQS(Eg);
        if q(end)>1
           if d1(end,2)<=d2(end,1)
                E3=getIdentity(E1,2,Eg);
           else E3=getIdentity(E2,1,Eg);
           end
           d3=getDimQS(E3); d3=d3(1,end);
        end
     end
  end

  complete=(Ntot>=d3);
  if complete, Ntot=d3; end

% NB! in case of pseudo-MPO, there is no MPO bond dimension yet!
% then when representing the partially contracted Hamiltonian
% as in (Ls,Rs) without intermediate indices, this is equivalent
% already to full 2-site expansion! // Wb,Sep08,23
% --> there is not much gain with CBE having pseudo MPO!
% to the extent that *this may be additional overhead !?

  iflag=0; if nargout>3
     if fullMPO && wCBE==0, iflag=2; else iflag=1; end
  end

  if iflag
     Iout=add2struct('-',k1,k2,fullMPO,complete,wCBE,facCBE,...
         'nrmPsi','nrmDD','dPsi');
     Iout.N2site=d3;
     Iout.Nkeep=[ Nkeep Ntot ];
     Iout.nk=[];
  end

  if complete || ~fullMPO || ~wCBE
     r=[ numel(X1.AK.Q), numel(X2.AK.Q) ];
     if any(r<3) || any(r>4) || any(r>3) && sum(r)~=7
        wbdie('invalid usage (got [%g,%g]-rank A-tensors)',r);
     end
     p3=[]; if r(1)>3, p3=[1 3 2]; end

     if iflag>1, X1_AK=X1.AK; X2_AK=X2.AK; end

     Psi=contract({E1,'*',X1.AK},{E2,'*',X2.AK},p3);

     X1.AK=E1; X1=updateHK(HAM,k1,'>>',X1);
     X2.AK=E2; X2=updateHK(HAM,k2,'<<',X2);

     if complete, return; end

     HPsi=get_HPsi(HAM,Psi,X1,X2,args{:});
     rPsi=rank(Psi);

     if rPsi==3
        [~,~,U ,I1]=svdQS(HPsi,1,'Nkeep',Ntot);
        [~,~,Vd,I2]=svdQS(HPsi,2,'Nkeep',Ntot); 
     else
        [U,S,Vd,I0]=svdQS(HPsi,2,'Nkeep',Ntot); 
     end

     nrm=normQS(Psi);
     if iflag, Iout.nrmPsi=nrm; end

     Psi=contract({U,'*',Psi},Vd,'*',p3);

     nrm(2)=normQS(Psi); q=-diff(nrm)/nrm(1);
     if abs(q)>1E-15
        if q<-1E-12 || q>1E-3
           wblog('WRN','adjusting norm by %.3g',q); end
        Psi=(nrm(1)/nrm(2))*Psi;
     end

     X1.AK=contract(U,X1.AK,[2 1]);
     X1.HK=contract(U,'*',{X1.HK,2,U},1);
     if ~fullMPO
        for i=1:numel(X1.OP), if X1.OP(i)
           X1.OP(i)=contract(U,'*',{X1.OP(i),2,U},1); end
        end
     end
     X1.Psi=[];

     X2.AK=contract(Vd,X2.AK);
     X2.HK=contract(Vd,'*',{X2.HK,2,Vd},1);
     if ~fullMPO
        for i=1:numel(X2.OP), if X2.OP(i)
           X2.OP(i)=contract(Vd,'*',{X2.OP(i),2,Vd},1); end
        end
     end
     X2.Psi=[];

     if iflag>1
        HPsi_=contract({X1.AK,{{U,'*',HPsi},Vd,'*'}},X2.AK);
        [~,Iout.nrmDD]=project_DD(HPsi_,X1_AK,X2_AK);
     end
     return
  end

  A1=X1.AK; A2=X2.AK;

  if ~itags2odir(X1.AK)
     [X1.AK,Psi]=ortho(X1.AK,2,'<<'); if iflag, Iout.k1=-Iout.k1; end
     A1(2)=X1.AK;
  elseif ~itags2odir(X2.AK)
     [X2.AK,Psi]=ortho(X2.AK,1,'<<'); if iflag, Iout.k2=-Iout.k2; end
     Psi=permute(Psi,[2 1]); A2(2)=X2.AK;
  else wbdie('failed to identify OC'); end

  if iflag, Iout.nrmPsi=normQS(Psi); end

  q={ X1.info.hconj, X2.info.hconj };
     if isequal(q{:}), q=q{1};
     else disp(q), wbdie('inconsistent hconj'); end
     if numel(q)~=1 || q<0, disp(q), wbdie('invalid hconj'); end
  hconj=q;

  dMPO=getDimQS(HAM.mpo(k1)); dmpo=dMPO(:,2);

  dPsi=getDimQS(Psi); rPsi=rank(Psi);
     if     rPsi==2, nPsi=1;            NPsi=1;
     elseif rPsi==3, nPsi=dPsi(:,end);  NPsi=nPsi(end);
     else wbdie('invalid Psi (rank %d)',rPsi); end
  nk=max(nx,ceil(Ntot/dmpo(end)));

  if rPsi==3, p3=[1 3 2]; else p3=[]; end

  if wCBE==1

     if rPsi==3
        [~,~,U ]=svdQS(Psi,1,'Nkeep',nk);
        [~,~,Vd]=svdQS(Psi,2,'Nkeep',nk); 
        Psi=contract({U,'*',Psi},Vd,'*',p3);
     else
        [U,S,Vd]=svdQS(Psi,2,'Nkeep',nk);
        Psi=diag(QSpace(S));
     end

     X1.AK=contract(U, X1.AK,[2 1]);
     X2.AK=contract(Vd,X2.AK);

  elseif wCBE==2

     [~,I1]=eigQS(contract(Psi,Psi,'!1*'),'Nkeep',nk);
     X1.AK=contract(I1.AK,X1.AK,[2 1]);

     [~,I2]=eigQS(contract(Psi,Psi,'!2*'),'Nkeep',nk); 
     X2.AK=contract(I2.AK,X2.AK);

     Psi=contract({I1.AK,'*',Psi},I2.AK,'*',p3);

  else wbdie('invalid switch (wCBE=%d)',wCBE);
  end

  X1=updateHK(HAM,k1,'>>',X1,'-x');
  X2=updateHK(HAM,k2,'<<',X2,'-x');

  mc=[ numel(X1.HK), numel(X2.HK) ];
  if diff(mc), wbdie('unexpected HK data (%dx%d)',mc); end
  mc=mc(1); k4=zeros(2,4); HD=QSpace;

  for j=1:mc
  for i=1:2
     if i==1, HD(i,j)=X1.HK(j); else HD(i,j)=X2.HK(j); end
     q=get_kidx(HD(i,j)); l=numel(q);
     if l~=4, wbdie('got rank-%d QSpace X%d.HK(%d)',l,i,j); end

     if j==1, k4(i,:)=q;
     elseif ~isequal(q.*[1 1 -1 1],k4(i,:))
        wbdie('inconsistent hconj data X%d.HK(j)',i,j);
     end
  end
  end

  Em=getIdentity(HAM.mpo(k1),2,'-0'); j=[]; HD_=HD;

  if diff(sign(k4(1,[2 3])))
               j=1; % insert Em'*Em at j==1 if mpo is incoming [k4(1,3)>0]
  elseif mc>1, j=2; % insert Em*Em' at j==2
  end
  if ~isempty(j), p=[1 2 4 3];
     ic={{'*'},{}}; if xor(j==1,k4(1,3)>0), ic=ic([2 1]); end
     for i=1:2
        HD(i,j)=contract(HD(i,j),Em,ic{i}{:}, p);
        if mc>1
           HD(i,j)=markitags(HD(i,j),3);
        end
     end
  end

  if mc==2
     for i=1:2, HD(i,1)=oplus(HD(i,1),HD(i,2),3); end
     HD=builtin('transpose',HD(:,1));
  end

  for i=1:2
     if i==1, AK=A1(end); else AK=A2(end); end
     [HD(i),nrm]=project_DD(HD(i),AK,3-i);
     if nrm(end)/nrm(1)<0.01
        HD(i)=project_DD(HD(i),AK,3-i);
     end
  end

  Ndisc=Ntot-dPsi(1,1:2);

  Q=contract({{HD(1),'14*',HD(1)},Psi},{HD(2),HD(2),'24*'},'12');
  H2L=contract(Q,Psi,'!1*');
  ELm=getIdentity(H2L,1,H2L,2,sprintf('Km%02d',k1));
  H2L=contract(ELm,'!3*',{H2L,'34',ELm});
  [U,S,Vd]=svdQS(H2L,2,'Nkeep',Ndisc(1));

  PL=contract(ELm,U);

  AD1=QSpace(svdQS(contract(HD(1),PL,[1 3 2]),2));

  e=normQS(contract(AD1,A1(end),'*'));
  if e>1E-12, wblog('WRN','got AD/AK overlap @ %.3g',e); end

  if NPsi<=1
     AD2=QSpace(svdQS(contract({PL,'*',Psi},HD(2)),1));
  else
     H2R=contract(Psi,'!2*',Q); ZR=getIdentity(H2R,4,'-0');
     ER=getIdentity(ZR,2,H2R,3,sprintf('Km%02d',k1));
     ER=contract(ZR,'*',ER);
     H2=contract(ER,'!3',{H2R,'34',ER,'*'});
     [U,S,Vd]=svdQS(H2,2,'Nkeep',Ndisc(2));

     PR=contract(ER,U,'*');

     AD2=QSpace(svdQS(contract(HD(2),PR,[3 1 2]),1));
  end
  e=normQS(contract(AD2,A2(end),'*'));
  if e>1E-12, wblog('WRN','got AD/AK overlap @ %.3g',e); end

  X1.AK=oplus(A1(end),AD1,2);
  X2.AK=oplus(A2(end),AD2,1);

  Psi=contract({X1.AK,'!2*',A1(1)},{X2.AK,'!1*',A2(1)},p3);

  if iflag
     Iout=add2struct(Iout,nk,dPsi);
     Iout.nrmPsi(2)=normQS(Psi);
     Iout.nrmDD=nrm;

     e=Iout.nrmPsi-1; if any(abs(e)>1E-12)
     wblog('WRN','got |Psi| = 1 %+.3g',min(e)); end
  end

  X1=updateHK(HAM,k1,'>>',X1);
  X2=updateHK(HAM,k2,'<<',X2);

end

% -------------------------------------------------------------------- %
% project into discarded space (by projecting out kept space)

function [HPsi,nrm]=project_DD(HPsi,A1,A2)

   q=(nargout>1); if q, nrm=norm(HPsi); end

   one_site=isnumeric(A2); ic='*'; p=[];
   if one_site, odir=A2;
      ic=sprintf(['!%d' ic],odir);
      if     odir==2, p=[1 3 4 2];
      elseif odir==1, p=[3 1 4 2];
      else wbdie('invalid odir=%d',odir); end
   end

      HPsi=HPsi-contract(A1,{A1,ic,HPsi},p); if q, nrm(2)=norm(HPsi); end
   if ~one_site
      HPsi=HPsi-contract({HPsi,A2,'*'},A2); if q, nrm(3)=norm(HPsi); end
   end

end

% -------------------------------------------------------------------- %
