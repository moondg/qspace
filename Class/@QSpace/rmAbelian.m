function A=rmAbelian(A,ix)
% function A=rmAbelian(A,ix)
%
%    Remove Abelian symmetry at position ix from existing set
%    of symmetries. Note that this requires that only a single
%    symmetry sector is present for any leg.
%    This shall be used only when tweaking the output of
%    getLocalSpace, e.g., by projecting to unique symmetry subspace.
%
% Wb,Mar15,24

% adapted from addSymmetry.m

  nA=numel(A); if ~nA, return; end
  if numel(ix)~=1, wbdie('invalid usage (index ix)'); end

  for k=1:nA, Ak=A(k);
     [r,sym]=getsym(Ak,'-r'); n=numel(r);
     if ix>n, wbdie('index for symmetry out of bounds (ix=%d/%d)',ix,n); end
     if r(ix), wbdie(...
       'invalid usage (got non-abelian %s for ix=%d)',sym{ix},ix);
     end

     r1=r; r1(r1<1)=1; j=sum(r1(1:ix));
     for l=1:numel(Ak.Q)
        q=Ak.Q{l}(:,j); if any(diff(q)), wbdie(['invalid usage ' ... 
           '(got %d sectors for symmetry %d/%d on leg %d)'], ix,n,l); end
        Ak.Q{l}(:,j)=[];
     end

     if ~isempty(Ak.info.cgr), Ak.info.cgr(:,ix)=[]; end
     if ~isempty(sym)
        sym(ix)=[];
        Ak.info.qtype=strjoin(sym,',');
     end
     A(k)=Ak;
  end
end

