function A=join(A,B)
% function A=join(A,idx)
% Wb,May19,10

  cgflag=gotCGS(A);

  if numel(A.Q)~=numel(B.Q) || ...
    ~isempty(A.Q) && size(A.Q{1},2)~=size(B.Q{1},2)

     if ~isempty(A.Q)
          i=[size(A.Q{1},2), size(B.Q{1},2)];
     else i=[-1,-1]; end

     wblog('ERR','rank/size mismatch (%d/%d; %d/%d)', ...
     numel(A.Q), numel(B.Q), i);
  end

  if xor(cgflag,gotCGS(B)) || cgflag && ~isequal(A.info.qtype,B.info.qtype)
     wblog('ERR','qtype mismatch (%d/%d)', cgflag, gotCGS(B));
  end

  if cgflag, A.info.cgr=[ A.info.cgr; B.info.cgr ]; end

  for p=1:length(A.Q), A.Q{p}=[ A.Q{p}; B.Q{p} ]; end

  A.data=[ reshape(A.data,1,[]), reshape(B.data,1,[]) ];

end

