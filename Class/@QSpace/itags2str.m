function ss=itags2str(A,varargin)
% function s=itags2str(A [,use_col])
%
%    return combined itags as string
%
% Wb,Jul13,23

% tags: itag2str
  ss=cell(size(A)); n=numel(A);

  if ~nargout && nargin<2, varargin={'QS:info'}; end
  for i=1:n
     ss{i}=itags_to_str(A(i).info.itags,varargin{:});
  end
  if ~nargout
     fprintf(1,'\n');
     fprintf(1,'    %s\n',ss{:}); fprintf(1,'\n');
     clear ss
  elseif n==1, ss=ss{1}; end

end

