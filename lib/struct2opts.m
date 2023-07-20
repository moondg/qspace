function f=struct2opts(S)
% Function: struct2opts(S)
%
%    convert given structure to option set
%    by pairing field and value.
%
% Wb,Dec10,07

  if ~isstruct(S)
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'); end
     return
  end

  f=fieldnames(S);

  for i=1:length(f)
  f{i,2}=getfield(S,f{i}); end

  f=reshape(permute(f,[2 1]),1,[]);

end

