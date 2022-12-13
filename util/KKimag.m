function [gi,Ii]=KKimag(om,gr0,varargin)
% function [gi,Ii] = KKimag(om,gr0 [,opts])
%
%    Use Kramers Kronig relation to obtain Im(G).
%    See kramers (mex routine) for options.
%
% Wb,Aug25,05, Wb,Jun04,06, Wb,Dec12,06

% latest matlab version: KKreal_061213.m

% kramers only accepts ascencding order for om
  [ox,i]=sort(om);
  if any(diff(i)~=1)
     om=ox; iflag=1;
     if size(gr0,1)==length(ox), gr0=gr0(i,:); else gr0=gr0(:,i); end
  else iflag=0;
  end

  [gi,Ii]=kramers(om,gr0,'r2i',varargin{:});

  if iflag
     if size(gr0,1)==length(ox), gr(i,:)=gr; else gr(:,i)=gr; end
  end

end

