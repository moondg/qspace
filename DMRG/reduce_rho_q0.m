function [Rho,sout]=reduce_rho_q0(AK,varargin)
% function Rho=reduce_rho_q0(AK,sc,HPsi)
% function Rho=reduce_rho_q0(Rho,  HPsi)
%
%    Reduce density matrix out of AK'*AK
%    assuming index order {KXX,Psi,KXX,Psi} within AK,
%    having global symmetry sectors specified with Psi.
%    HPsi is used to project to low-energy (mixed) space
%    within each symmetry sector in Psi to reduce the object
%    size of Rho in the presence of multiple global multiplets.
%    Thus Rho becomes block-diagonal in Psi and only carries
%    a single value for each symmetry.
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
     [sc,HPsi]=deal(varargin{:}); Rho=QSpace;
     if ischar(sc) 
        if ~isempty(regexp(sc,'^\d+$'))
           Rho=QSpace(contractQS(AK,sc,AK,[sc '*']));
        elseif ~isempty(regexp(sc,'^!\d+$'))
           Rho=QSpace(contractQS(AK,sc,AK,'*'));
        end
     end
     if ~Rho, wbdie('failed to contract Psi -> Rho'); end
  else
     Rho=AK; HPsi=varargin{1};
  end

  if isempty(Rho), return; end

  if numel(Rho.Q)==2, sz=sizeof(Rho); 
     if sz>2E5
        [~,IR]=eigQS(Rho); Rho=QSpace(IR.EK);
     end
  elseif numel(Rho.Q)==4 && isa(HPsi,'QSpace')
     if ~isfield(HPsi,'Q') || ~isfield(HPsi,'data')
        wbdie('invalid input HPsi (%s)',class(HPsi));
     end
     r=numel(HPsi.Q); t=getitags(HPsi,1);
     if r~=2 || isempty(regexp(t,'Psi')) 
        wbdie('invalid input HPsi (rank-%d with itag %s)',r,t);
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

     RQ=HPsi;
     for i=1:numel(RQ.data)
        Hi=HPsi.data{i};
        [U,E]=eig(Hi+Hi'); E=diag(E)/2; beta=1E-3*diff(E([1 end]));
        R=exp(-beta*E); RQ.data{i}=diag(R/sum(R));
     end

     [q,d,dc]=getQDimQS(RQ,1); dc=prod(dc,2);
     [i1,i2,Im]=matchIndex(RQ.Q{1},q);
     for i=1:numel(i1)
        RQ.data{i1(i)}=RQ.data{i1(i)}/dc(i2(i));
     end

     [iR,iQ,Im]=matchIndex(Rho.Q{end},RQ.Q{1});
     if isempty(Im.ix1)
        for i=1:numel(iR), ir=iR(i);
           Ri=Rho.data{ir};
           Ri=contract(Ri,RQ.data{iQ(i)},[3 4],[1 2]);
           Rho.data{ir}=permute(Ri,[1 2 4 5 3]);
        end
     else wblog('WRN','got missing symmetry sectors in HPsi'); end
  end
end

