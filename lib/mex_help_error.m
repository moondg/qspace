function mex_help_error()
% function mex_help_error()
% Wb,Jan16,19

  S=dbstack;
  l=length(S); if l<2, return; end
  n=S(2).name;

  wbdie([ 10 ...
  '   ERR ' n '.m is the help script for the mex-function ''' n ''''   10  ...
  '   ERR Ensure that the mex-routine has predence on the matlab path, i.e.,' 10  ...
  '   ERR that it resides in the same path as the mex-file or later.' ...
  ])

end

