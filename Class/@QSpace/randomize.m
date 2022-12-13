function A=randomize(A,varargin)
% function A=randomize(A [,nrm,opts])
%
%    Randomize A.data using randn()
%    while preserving overall norm (unless explicitly specified).
%
% Options
%
%    -z   use complex random numbers (default: real)
%
% Wb,Nov11,09

  getopt('init',varargin);
     zflag=getopt('-z');
  nrm=getopt('get_last',-1);

  if nrm<=0, nrm=normQS(A); end

  n=numel(A.data);
  for i=1:n
     A.data{i}=randn_z(size(A.data{i}),zflag);
  end

  nrm = [nrm, normQS(A)];
  if ~all(isfinite(nrm) & nrm), wblog('WRN','got norm %g/%g !?',nrm); end
  nrm=nrm(1)/nrm(2);

  for i=1:n, A.data{i}=A.data{i}*nrm; end

end

% -------------------------------------------------------------------- %
function X=randn_z(sz,zflag)

   X=randn(sz);
   if zflag, X=complex(X,randn(sz)); end

end

% -------------------------------------------------------------------- %

