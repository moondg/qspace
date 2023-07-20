function [rr,tt] = gotITags(A)
% function [r,t] = gotITags(A)
%
%    return number of itags (if any).
%    second return argument contains listing of ilables as cell array.
%
% Wb,Aug04,12

  rr=0; tt={};
  if ~isfield(struct(A),'info'), return; end

  nA=numel(A); rr=repmat(-1,size(A)); tt=cell(size(A));
  for k=1:nA
     Ak=A(k); if ~isfield(Ak.info,'itags'), continue; end
     tk=Ak.info.itags; if isempty(tk), continue; end

     if ischar(tk)
        tk=strread(tk,'%s','delimiter',',;');
     end
     r=[numel(tk), numel(Ak.Q)];
     if diff(r)==0, rr(k)=r(1); tt{k}=tk;
     else wbdie('invalid number of itags (%d/%d)',r); end
  end
  if nA==1, tt=tt{1}; end

end

