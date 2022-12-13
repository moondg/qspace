function ies=isemtpystruct(s)
% Function i=isemtpystruct(s)
%
%    checks whether given input structure has the form
%    of a structure with all fields being empty.
%
% Wb,Aug07,07

  ies=0;

  if nargin~=1
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  if ~isstruct(s), return; end

  f=fieldnames(s); nf=length(f);
  for i=1:nf
     if ~isempty(getfield(s,f{i})), return; end
  end

  ies=1;

end

