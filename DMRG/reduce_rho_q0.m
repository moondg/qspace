function [Rho,sout]=reduce_rho_q0(AK,varargin)
% Usage #1: [Rho,sout]=reduce_rho_q0(AK,sc,H0)
%
%    This contracts AK with itself using sc (= string for ic).
%    With the result being Rho, this proceeds with usage #2.
%
% Usage #2: [Rho,sout]=reduce_rho_q0(Rho,H0)
%
%    Reduce density matrix out of AK'*AK (e.g. see usage #1)
%    assuming index order {KXX,Psi,KXX,Psi} within AK,
%    having global symmetry sectors specified with Psi.
%    H0 (required) represents Hamiltonian in the Psi-basis.
%    It is used to project to low-energy (mixed)
%    space within each symmetry sector in Psi to reduce the object
%    size of Rho in the presence of multiple global multiplets.
%    Thus Rho becomes block-diagonal in Psi and only carries
%    a single value for each symmetry.
%
%    Output string sout only contains warning if any
%    in which case if sout is used in the caller,
%    a warning log output is suppressed.
%
% Examples
%
%    rho=reduce_rho_q0(Psi,'12',Eg);
%    rho=reduce_rho_q0(Rho,Eg);
%
% Wb,Apr06,23

% tags: getRho, getrhoQS

  sout='';
  if nargin<2 || nargin>3, helpthis
     if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if nargin==3
     [sc,H0]=deal(varargin{:}); Rho=QSpace;
     if ischar(sc) 
        if ~isempty(regexp(sc,'^\d+$'))
           Rho=QSpace(contractQS(AK,sc,AK,[sc '*']));
        elseif ~isempty(regexp(sc,'^!\d+$'))
           Rho=QSpace(contractQS(AK,sc,AK,'*'));
        end
        rk=numel(Rho.Q); if rk==4
         % assuming OC for AK, all indices are incoming except for Psi
         % Xg'X'g -> XgX'g' with X \in LRs
         % Rho=permute(Rho,'1432'); % Wb,Apr16,24
         % actually leave as is => reductionb elow generates XX'g'(g) order
         % with the last index just an auxiliary singleton
        end
     end
     if ~Rho, wbdie('failed to contract Psi -> Rho'); end
  else
     Rho=AK; H0=varargin{1};
     rk=numel(Rho.Q);
  end

  if isempty(Rho), return; end

  if rk==2
     sz=sizeof(Rho); if sz>2E5
        [~,IR]=eigQS(Rho); Rho=QSpace(IR.EK);
     end
     return
  elseif rk~=4 || ~isfield(H0,'Q') || ~isfield(H0,'data')
     return
  end

  if ~isa(H0,'QSpace'), H0=QSpace(H0); end

  rk=rank(H0); t=getitags(H0,1);
  if rk~=2 || isempty(regexp(t,'Psi')) 
     wbdie('invalid input H0 (rank-%d with itag %s)',rk,t);
  end

  E2=getIdentity(Rho,2);
  E4=getIdentity(Rho,4);
  if ~isequal(E2,E4), E3=getIdentity(Rho,3);
     q=getDimQS(Rho); q=sprintf(' x %d',q(end,:));
     s={ q(4:end), sizeof(Rho,'-s') };
     if isequal(E3,E4)
        sout=sprintf('already got Rho in reduced format (%s @ %s)',s{:});
        if nargout<2, wblog(' * ',sout); end
     else wblog('WRN','got Rho of dimension %s @ %s',s{:});
     end
     return
  end

  i=find(sum((Rho.Q{2}-Rho.Q{4}).^2,2)==0);
  Rho=permute(getsub(Rho,i),[1 3 2 4]);

% NB! previously projected to lowest eigenstate in each symmetry sector
% --> reduced to ground state only if multiple states within a single
%     symmetry sector were targeted
% --> keep all (do not project to lowest in sector!) // Wb,Apr16,24
% ==> just keep diagonal contribution for ALL eigenstates
%     storing it in compact format as 3rd index (by keeping its symmetry
%     label, also 4th index must be kept still, even if singleton)

  nQ=numel(H0.data); nfull=0;
  for k=1:nQ, ek=H0.data{k};
     if numel(ek)>1 && ~diff(size(ek)), ek=eig(Hk+Hk'); end
     if any(diff(sort(ek))<1E-6), nfull=nfull+1; continue; end

     for i=matchIndex(Rho.Q{3},H0.Q{1}(k,:))
        Ri=Rho.data{i};
        for j=1:size(Ri,3), Ri(:,:,j,1,:)=Ri(:,:,j,j,:); end
        Rho.data{i}=Ri(:,:,:,1,:);
     end
  end

  if ~nfull, q=getDimQS(E2);
     q=[ q(end), trace(contract(Rho,getIdentity(Rho,[1 2]))) ];
     e=norm(diff(q));
     if e>1E-8, save2tmp
        wbdie('failed to reduce Rho (%g / %g @ %.2g)',q,e);
     end
  end

end

