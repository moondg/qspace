function i = gotCGS(A)
% function i = gotCGS(A)
% Wb,May16,10

  i=0;

  if nargin~=1, wblog('ERR','%s() invalid usage',mfilename); end
  if isfield(struct(A),'info') && isfield(A.info,'cgr') && ~isempty(A.info.cgr)
     i=1;
     if numel(A.Q)>2
        s=A.info.qtype;
        s=strrep(s,'A,',''); s=strrep(s,'SU2,','');
        if ~isempty(s), i=2; end
     end

     s=size(A.info.cgr);
     n=length(find(A.info.qtype==','))+1;
     if ~isempty(A.Q) && ~isempty(A.Q{1})
        if s(1)~=numel(A.data) || length(s)>2 ...
           || ~isempty(A.Q) && s(2)>size(A.Q{1},2) || s(2)~=n ...
           || isempty(A.Q) && s(2)~=0, s
           wblog('WRN','QSpace::%s() invalid info.cgr data',mfilename);
        end
     else
        if s(1)~=numel(A.data) || s(1)>1
           wblog('WRN','QSpace::%s() invalid info.cgr data',mfilename);
        end
     end
  end

end

