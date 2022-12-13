function fprintf2(fid,varargin)
% function fprintf2(fid,varargin)
%
%    duplicate output to all fid's specfied
%
% Wb,Aug18,16

  for i=1:numel(fid)
  fprintf(fid(i),varargin{:}); end

end

