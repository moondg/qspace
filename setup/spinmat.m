function S=spinmat(varargin)
% Function: S=spinmat([d, OPTS])
%
% Get generalized Pauli spin matrizes for any dimension d
% Options
%
%      d       dimension of spin matrizes (e.g. 2 for s=1/2)
%     '-pauli' get pauli matrizes (fac=2)
%     '-sp'    keep matrizes sparse
%     '-pm'    return Sp,Sm,Sz instead of Sx,Sy,Sz
%     '-sym'   return spin-1 IROP ~ (Sp,Sz,Sm)
%     'fac',.. overall factor to multiply (1.)
%
% Ref. Griffiths p.195
% Wb,Jun15,07

  d=2; fac=[];

  getopt('init',varargin);
     spflag =getopt('-sp');
     pmflag =getopt('-pm');
     pauli  =getopt({'-pauli','pauli'}); if pauli, fac=2; end
     symflag=getopt('-sym'); if ~symflag
     fac    =getopt('fac',fac); end
  d=getopt('get_last',d);

  if pauli && (d~=2 || pmflag || symflag) || pmflag && symflag
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'); else return; end
  end

  s=(d-1)/2; m=(s:-1:-s)';
  Sz=spdiags(m,0,d,d);

  Sm=spdiags(sqrt(s*(s+1) - m.*(m-1)),-1, d, d);
  Sp=Sm';

  if symflag
     S = {-sqrt(0.5)*Sp,Sz,+sqrt(0.5)*Sm};
  elseif pmflag, S = {Sp,Sm,Sz};
  else
     Sx=0.5 *(Sm+Sp);
     Sy=0.5i*(Sm-Sp); S = {Sx,Sy,Sz};

   % check commutator relations (safeguard)
   % e = [ norm(full(comm(Sx,Sy)-1i*Sz))
   %       norm(full(comm(Sy,Sz)-1i*Sx))
   %       norm(full(comm(Sz,Sx)-1i*Sy))
   %       norm(full(Sx+1i*Sy-Sp))
   %       norm(full(Sx-1i*Sy-Sm))
   % ]; if any(e), e, wblog('ERR','invalid spin operators'); end
  end

  if ~isempty(fac)
  for p=1:length(S), S{p}=fac*S{p}; end, end

  if ~spflag
  for p=1:length(S), S{p}=full(S{p}); end, end

end

