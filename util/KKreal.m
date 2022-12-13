function [gr,Ir] = KKreal(om,gi0,varargin)
% function [gr,Ir] = KKreal(om,gi0 [,opts])
%
%    Use Kramers Kronig relation to obtain Re(G).
%    See kramers (mex routine) for options.
%    Gi is assumed in a continuous (non-discrete) form.
%    (and keep it this way! ;)
%
% Wb,Aug25,05, Wb,Jun04,06 Wb,Dec12,06

% latest matlab version: KKreal_061213.m
% gr = KKreal_061213(om,gi0,varargin{:});
% return

% kramers only accepts ascencding order for om
  [ox,i]=sort(om);
  if any(diff(i)~=1)
     om=ox; iflag=1;
     if size(gi0,1)==length(ox), gi0=gi0(i,:); else gi0=gi0(:,i); end
  else iflag=0;
  end

  [gr,Ir]=kramers(om,gi0,varargin{:});

  if iflag
     if size(gi0,1)==length(ox), gr(i,:)=gr; else gr(:,i)=gr; end
  end
end

