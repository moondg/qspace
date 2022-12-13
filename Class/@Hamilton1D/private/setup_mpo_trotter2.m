function HAM=setup_mpo_trotter2(HAM,varargin)
% function HAM=setup_mpo_trotter2(HAM [,dt,opts])
%
%    This generates 2nd-order Trotter time step exp(-dt*H)
%    thus accepting nearest-neighbor terms in the Hamilonian only.
%    The result is stored in the structure HAM.user.trotter
%    with fields
%
%     .mpo   actual MPO of dimension dloc^2
%     .info  such as time step dt etc.
%
%     .data  if time evoultion is performed in memory
%     .mat   file tag if time evoultion is stored in file system
%            (NB! Trotter uses files separate from the DMRG data,
%            with file tag <DMRG>-dt*)
% Options
%
%    dt      time step tau to tak (0.005)
%   'e0',..  energy/site reference to take
%            default: ground state energy/site from DMRG info data
%
% Wb,May30,19

% Adapted from $HAM/Lab/tst_unitary_mpo.m
% See also MATH/trotter.tex

% for last version for translationally invariant Hamiltonian
% with ML, MC, MR see Archive/setup_mpo_trotter2_210522.m

  getopt('init',varargin);
     e0=getopt('e0',[]);
  dt=getopt('get_last',0.005);

  L=length(HAM.mpo);

  if isempty(e0)
     q=load_dmrg_info(HAM,'E0');
     q=reshape(q,[],1); q=[min(q),mean(q),std(q)];

     e0=q(1); 

     x=q(end-1)/q(end); if x<1E-2
          wblog(' * ','using GS energy E0/L=%g',e0);
     else wblog('WRN','using GS energy E0/L=%g (@ %.3g !?)',e0,x);
     end
  elseif numel(e0)~=1 || ~isnumeric(e0)
     wberr('invalid usage (scalar required for e0)');
  end

  ph0=exp(+1i*dt*e0);

  q=[HAM.mpo.stype];
  if any(diff(q)) || numel(q)~=numel(HAM.mpo), disp(unique(q));
     wberr('got invalid / unexpected stype !?');
  end

  HH=HAM.info.HH; Hloc=[]; Hcpl=[];

  [dx,Ix,Dx]=uniquerows(abs(HH(:,3)-HH(:,1)));
  if any(dx~=0 & dx~=1)
     wberr('only up to nearest-neighbor terms allowed with TR2');
  end

  Hloc=cell(1,L); nloc=0;
  Hcpl=cell(1,L-1);

  for i=1:numel(dx), Hx=HH(Ix{i},:); l=dx(i);
     [k,Ik,Dk]=uniquerows(Hx(:,1));
     if l==0, nloc=nloc+1;
       for j=1:numel(k), q=Hx(Ik{j},:);
          if ~isequal(uniquerows(q(:,[1 3]))-k(j),[0 0]) || ...
             ~isequal(q(:,2),q(:,4)), q, wberr('unexpected data');
          end
          Hloc{k(j)}=q(:,[2 5:end]);
       end
     elseif l==1
       for j=1:numel(k), q=Hx(Ik{j},:);
          if ~isequal(uniquerows(q(:,[1 3]))-k(j),[0 1]), q
             wberr('unexpected data'); end
          Hcpl{k(j)}=q(:,[2 4:end]);
       end
     else wberr('invalid usage (dx=%g)',l); end
  end

  E1=HAM.oez(1).op;
  A2=getIdentityQS(E1,E1,[1 3 2]); A2.info.itags={'s1','K2*','s2'};

  uform=1;
  if nloc, for i=2:numel(Hloc)
     if ~isequal(Hloc{i},Hloc{1}), uform=0; break; end
  end, end
  if uform, for i=2:numel(Hcpl)
     if ~isequal(Hcpl{i},Hcpl{1}), uform=0; break; end
  end, end

  Z=getIdentityQS(A2,3,'-0');
  Q=getIdentityQS(Z,1,Z,2,[1 3 2]);
  X2=untag(QSpace(contractQS(Q,3,Z,'2*',[1 3 2])));

  q=getDimQS(A2); dloc=q(end,1);

  HAM.user.trotter.dt=dt;
  HAM.user.trotter.mpo=QSpace(1,L);

  for k=[1 L, 2:L-1]
     bflag=0; if k<=1, bflag=-1; elseif k>=L, bflag=+1; end

     if k==1
        M=get_mpo_k(HAM, 0.,[], 0.5,Hloc{k}, [], bflag,dt,E1,A2,X2);
     else
       if bflag>0, w2=1; else w2=0.5; end
       [M,H2,U4]=get_mpo_k(HAM, ...
          0.5,Hloc{k-1}, w2,Hloc{k}, Hcpl{k-1}, bflag,dt,E1,A2,X2);
     end

     if 1
        q=normQS(M)^2; q=[q, dloc^3, q-round(q)];
        if bflag, q(2)=dloc^2; end
        if abs(diff(q(1:2)))>1E-12
        wblog('WRN','|M%02g|^2 =%4g/%3g  @ %8.2g',k,q); end
     end

     M.info.itags=set_itags(M.info.itags,k);
     HAM.user.trotter.mpo(k)=ph0*M;
  end

% -------------------------------------------------------------------- %
  if uform, Hloc=Hloc{1}; Hcpl=Hcpl{1}; end
  if ~nloc, Hloc={}; end

  HAM.user.trotter.info=add2struct('-','ttot=0',dt,Hloc,Hcpl,H2,U4);

end

% -------------------------------------------------------------------- %
function [M,H2,U4]=get_mpo_k(HAM, w1,h1,w2,h2, Hcpl, bflag,dt,E1,A2,X2)

  persistent Mlast

  args={w1,h1,w2,h2,Hcpl,bflag};
  if numel(Mlast)==2 && isequal(Mlast{1},args)
       M=Mlast{2}; return
  else Mlast=args; end

  if isempty(h1), w1=0; end
  if isempty(h2), w2=0; end; H2=QSpace;
  zero=mean(abs(HAM.info.HH(:,end)))*1E-32;

  if bflag>0 && w2 && w2~=1, q=[w2 1];
     wblog('WRN','for right boundary, setting w2=%g -> %g',q);  w2=q(2);
  end

  for iloc=1:2
     if iloc==1, w=w1; else w=w2; end;  if ~w, continue; end
     if iloc==1, Hloc=h1; ic='-op:s1'; 
     else        Hloc=h2; ic='-op:s2'; end
     if size(Hloc,2)~=2, wberr('invalid usage'); end
     Hloc(:,2) = w * Hloc(:,2);

     H1=QSpace; err=0;

     for i=1:size(Hloc,1), X=HAM.ops(Hloc(i,1));
        h=Hloc(i,2)*X.op;
           if X.hconj, h=h+h'; end
           if X.fermionic, err=err+1; end
        H1=H1+h;
     end

     r=rank(H1); e=normQS(H1-H1');
     if r~=2, wberr('invalid rank-%g for local ops',r); end
     if err || e>1E-12, wberr('got non-hermitian local ops (e=%g)',e);
     elseif norm(H1)>1E-14
          H2=H2+contractQS(A2,'!2*',{H1,ic,A2});
     else H1=QSpace; if iloc==1, w1=0; else w2=0; end
     end
  end

  if isempty(Hcpl)
     if w1 || bflag>=0, wberr('invalid usage'); end
     M=permute(setitags(X2,'-op:s1','mpoR'),[3 1 2]);
     if w2
        U1 = H1 + zero*E1;
        for i=1:numel(U1.data)
           U1.data{i}=expm(-(0.5i*dt)*U1.data{i});
        end
        M=contract({M,U1,'-op:s1','!1'},U1,'-op:s1','!2');
     end
     M=skipzeros(M); return

  elseif size(Hcpl,2)~=3, wberr('invalid usage'); end

  for i=1:size(Hcpl,1), X=HAM.ops(Hcpl(i,[1 2]));

     isferm=[X.fermionic]; hconj=[X.hconj];
     if norm(diff(isferm)) || norm(diff(hconj)), isferm, hconj
        wberr('got fermionic/hconj mitmatch');
     else isferm=isferm(1); hconj=hconj(1); end

     if isferm, X(2).op=HAM.oez(2).op*X(2).op; end

     Q=Hcpl(i,3)*QSpace(...
        contractQS(A2,'!2*',{{X(1).op','-op:s1',A2},X(2).op,'-op:s2'}));
     if hconj, Q=Q+Q'; end

     H2=H2+Q;
  end

  e=norm(H2-H2'); if e>1E-12
     wberr('got non-hermitian Hamiltonian (e=%.2g)',e); end

  U2 = H2 + zero*QSpace(getIdentityQS(A2,2));
  for i=1:numel(U2.data)
     q=U2.data{i};
       e=norm(q-q')/max(1,norm(q));
       if e>1E-14, wblog('ERR','got non-hermitian H2 @ %.3g',e); end

     U2.data{i}=expm(-(0.25i*dt)*(q+q')); % = (dt/2)*(q+q')/2
  end
  U4=skipzeros(QSpace(contractQS({A2,U2},A2,'!13*')));

  if bflag>0
     if w2 && w2~=1, wberr('invalid usage'); end
     M=contract(X2,'-op:s1:mpoL','*',{U4,'!12',U4});
  else
     M=contract(X2,'-op:s1:mpoL','*',{U4,'!12',...
     {X2,'-op:s2:mpoR',U4,'!4'}}, [1 3 2 4]);
  end

  M=skipzeros(M); Mlast={M,Mlast};

end

% -------------------------------------------------------------------- %
function tt=set_itags(tt,k)

   tt=regexprep(tt,'^mpo\w*$',  sprintf('m%02g',k-1),'ignorecase');
   tt=regexprep(tt,'^mpo\w*\*$',sprintf('m%02g*', k),'ignorecase');
   tt=regexprep(tt,'^s\d*',     sprintf('s%02g',  k),'ignorecase');

end

% -------------------------------------------------------------------- %

