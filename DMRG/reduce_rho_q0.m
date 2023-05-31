function [Rho,sout]=reduce_rho_q0(AK,varargin)
% Usage #1: [Rho,sout]=reduce_rho_q0(AK,sc,H0)
%
%    This contracts AK with itself using sc.
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
     end
     if ~Rho, wbdie('failed to contract Psi -> Rho'); end
  else
     Rho=AK; H0=varargin{1};
  end

  if isempty(Rho), return; end

  if numel(Rho.Q)==2, sz=sizeof(Rho); 
     if sz>2E5
        [~,IR]=eigQS(Rho); Rho=QSpace(IR.EK);
     end
     return
  elseif numel(Rho.Q)~=4 || ~isfield(H0,'Q') || ~isfield(H0,'data')
     return
  end

  if ~isa(H0,'QSpace'), H0=QSpace(H0); end

  r=rank(H0); t=getitags(H0,1);
  if r~=2 || isempty(regexp(t,'Psi')) 
     wbdie('invalid input H0 (rank-%d with itag %s)',r,t);
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

  E=diag(H0,'-d'); E=[min(E), max(E)];
  H0=H0-E(1);
  beta=1E-3*diff(E);

  [q,~,dc]=getQDimQS(H0,1); dc=prod(dc,2); Q1=H0.Q{1};
  [i1,i2,Im]=matchIndex(Q1,q,'-s');
  if ~isequal(i1,1:size(Q1,1)), wbdie('got Q-label mismatch !?'); end
  dc=dc(i2);

  RQ=H0;
  for i=1:numel(RQ.data)
     E=H0.data{i};
     if ~isvector(E), [U,E]=eig(E+E'); E=diag(E)/2; end
     R=exp(-beta*E); RQ.data{i}=diag(R/(dc(i)*sum(R)));
  end

  [iR,iQ,Im]=matchIndex(Rho.Q{end},H0.Q{1});
  if isempty(Im.ix1)
     for i=1:numel(iR), ir=iR(i);
        Ri=Rho.data{ir};
        Ri=contract(Ri,RQ.data{iQ(i)},[3 4],[1 2]);
        Rho.data{ir}=permute(Ri,[1 2 4 5 3]);
     end
  else wblog('WRN','got missing symmetry sectors in H0'); end

end

