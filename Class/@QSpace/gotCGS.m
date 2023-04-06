function [i,rsym] = gotCGS(A)
% function [i,rsym] = gotCGS(A)
% Wb,May16,10

  i=0; rsym=[];
  if ~isfield(A,'info') || ~isfield(A.info,'cgr') || ...
      isempty(A.info.cgr)
      if ~isempty(A.Q), rsym=zeros(1,size(A.Q{1},2)); end
      return
  end

  i=1;
  if numel(A.Q)>2
     s=regexprep(A.info.qtype,'\<(A|Z\d+|SU2)\>','');
     if ~isempty(regexp(s,'[\d\w]')), i=2; end
  end

  s=size(A.info.cgr);
  n=length(find(A.info.qtype==','))+1;
  if ~isempty(A.Q) && ~isempty(A.Q{1})
     if s(1)~=numel(A.data) || length(s)>2 ...
        || ~isempty(A.Q) && s(2)>size(A.Q{1},2) || s(2)~=n ...
        || isempty(A.Q) && s(2)~=0, s
        wbwrn('invalid info.cgr data');
     end
  elseif s(1)~=numel(A.data) || s(1)>1
     wbwrn('invalid info.cgr data');
  end

end

