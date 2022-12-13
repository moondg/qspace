function M=blkdiag(A,varargin)
% function M=blkdiag(A,varargin)
% Wb,Aug10,16

  if numel(A)~=1, error('Wb:ERR','invalid usage'), end

  if ~isempty(A.Q)
     if numel(A.Q)~=2 || ~isdual_(A.Q{1},A.Q{2})
        wblog('WRN','got rank-%g object with off-diagonal blocks !?',...
        numel(A.Q));
     end
  end

  M=A.data; M=blkdiag(M{:});

end

