function varargout=mpsGetMaxDim(A)
% mpsGetMaxDim - get maximum D for every dimension of given MPS state
% Usage: D=mpsGetMaxDim(A)
%
%    where D is a vector of length rank(A)
%
% Wb,Aug08,06

  varargout=cell(1,nargout);

  if ~isa(A,'QSpace')
     if isempty(A), D=[]; DD=[]; return; end
     if ~isstruct(A) || ~isfield(A,'Q') || ~isfield(A,'data')
     eval(['help ' mfilename]); return; end

     A=QSpace(A);
  end

  [varargout{:}]=dsize(A);

end

