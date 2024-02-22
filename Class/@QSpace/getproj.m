function P=getproj(H,varargin)
% function P=getproj(H [,eps,opts])
%
%    get projector with respect to H into ground state space
%    (diagonalize if necessary); eps determines the window wrt.
%    ground state energy, within which state are kept (1E-8).
%
% Options
%
%    '-q'  quiet mode
%
% Wb,Aug30,12

% adapated from getrhoQS.m

  getopt('init',varargin);
     qflag=getopt('-q');
  eps=getopt('get_last',1E-8);

  isd=isdiag(H);
  if isd
     if isd==1
        P=struct(diag(H));
     else
        if isd~=2, wbdie('invalid isd=%g',isd); end
        P=struct(H);
     end
     ee=cat(1,P.data{:});
  else
     if ~qflag
     wblog(' * ','diagonalizing input H'); end
     [ee,Is]=eigQS(H); P=Is.EK; ee=ee(:,1);
  end

  eref=min(ee)+eps; d=0;

  for i=1:numel(P.data)
     x=P.data{i};
       if isempty(x), wbdie('got empty data'); end
     s=size(x); s=s(find(s>1)); if isempty(s), s=1; end
       if numel(s)>1, wbdie('got matrix !??'); end

     j=find(x<eref); n=numel(j);
     if n, d=d+n;
        P.data{i}=full(sparse(j,1:n,ones(1,n),s,n));
     else
        P.data{i}=0;
     end
  end

  P=skipZerosQS(P);
  if ~isd
       P=QSpace(contractQS(Is.AK,2,P,1));
  else P=QSpace(P); end

end

