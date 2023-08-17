function [Rho,Iout,IL,E3]=initRho(HAM,tau,varargin)
% function [Rho,Iout]=initRho(HAM,tau [,p])
%
%    Obtain Rho := exp(-tau*H) with tau << 1
%    based on a power expansion up to the optional power p specified.
%    without(!) normalization, i.e., without 1/Z factor.
%    This starts from the full MPO for the HAM (thus required in input).
%
% Wb,Jul27,23

% for testing purposes:
% wsys='Spin'; L=16; use_mpo=1; tflag=2; dberr; tst_Hamilton1D
% [Rho,Iout]=getRho(HAM,0.001);

  if nargin<2, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  t0=tic();

  getopt('init',varargin);
     vflag=getopt('-v');
     tflag=getopt('-t');
     kflag=getopt('-k');

     rtol =getopt('rtol', []); % 1E-24 // defaults in MPO_add.m
     Nkeep=getopt('Nkeep',[]); % -1
     nsw  =getopt('nsw',  []); %  3

   % (estimate for) ground state energy to subtract (e0 = E0/L)
     e0=getopt('e0',0);
     pwr=getopt('pwr',[]);

  if isempty(pwr)
       pwr=getopt('get_last',[]);
  else getopt('check_error'); end

  if ~using_full_MPO(HAM)
     wbdie('invalid usage (full MPO required for HAM)'); end

  [Dtot,dloc]=MPO_get_Dtot(HAM); % requires/ensures full MPO

  osw={};
% NB! instead of normalizing HAM.mpo, adjusting rtol in MPO_add()
% see also comments in MPO_get_Dtot.m
% WRN! Dtot=dloc^L quickly becomes astronomically large,
% and so does |Rho| ~ |Id| for tau<<1!
  if ~isempty(rtol) && rtol>0, setopts(osw,rtol); end % rtol*Dtot
  if ~isempty(Nkeep) && Nkeep>0, setopts(osw,Nkeep); end
  if ~isempty(nsw), setopts(osw,nsw); end

  L=numel(HAM.mpo);
  use_hconj=HAM.info.mpo.use_hconj;

% NB! eventually for large beta (larger n in XTRG):
% beta_n = (2^n)*tau >> 1, such that [rho(beta_n)] -> exp(-beta_n E_0)
% For example, for Heisenberg, with E0 ~ L*(-0.44) for beta=50 and L=16:
% exp(-beta*L*|e0|) = exp(+342) = 7.44E+152 (!)
% Similarly, if E0 is positive, like e0 -> -e0 above,
% exp(-beta*L*e0) = exp(-342) = 2.96E-149 (!)
% ==> best to subtract (good estimate for) ground state energy
%     even though for tau<<1 initially this has near negligible effect
%     since rho(tau<<1) ~ 1 - tau*H is dominated by Id.
% Wb,Jul29,23

  Eref=L*e0;
  if Eref
     wblog(' * ','subtracting energy reference e0=%.4g (Eref=L*e0=%.4g)',e0,Eref); 
  end

  if use_hconj
     fac=[ 1 +1; 1 -1]; % H_full = H+H'
     if Eref, fac=[-Eref 0; fac]; end % H_full -> (H_full-Eref) // e0
     [Ham,IH,IL,E3]=MPO_add(HAM.mpo,fac,osw);
  elseif Eref
       [Ham,IH,IL,E3]=MPO_add(HAM.mpo,[-Eref 0; 1 +1],osw); % H -> H - Eref
  else [Ham,~, ~, E3]=MPO_add(HAM.mpo,osw);
  end

% obtain energy fluctuation to check bounds on tau
  E1=MPO_trace(Ham,E3);
  E2=MPO_trace(Ham,2,E3);

  q2=[ E2, MPO_norm2(Ham) ]; % check hermiticity of Ham // safeguard
  e2=abs(diff(q2))/norm(q2);
  if e2>1E-12
     if e2<1E-9
          wblog('WRN','Ham hermitian @ %.3g',e2); 
     else wbdie('Ham not hermitian @ %.3g',e2); end
  end

  dE=sqrt(abs(E2)/(Dtot*L));
  % e.g. Heisenberg (J=1, L=16) => dE=0.4193 with E1=0
  % to be compared to e1 = -0.64342365 => dE ~ (2/3) |e1|

  if tflag  % run consistency checks
     q0=[ MPO_trace(HAM.mpo,2), MPO_norm2(HAM.mpo) ]; % MPO_trace(HAM.mpo),  
     e1=norm(diff(q0))/norm(q0); e=0;
     if use_hconj
         if e1<1E-12, e=1; wblog('WRN','already got hermitian HAM @ %.3g',e1); end
     elseif e1>1E-12, e=2; wblog('ERR','got non-hermitian HAM @ %.3g',e1); end

     if use_hconj
      % | H+H' |^2 = 2*|H|^2 + 2*real(tr(H^2))
        q = [ 2*(q0(1) + real(q0(2))), E2 ];
        e3=norm(diff(q))/norm(q);
        if e3>1E-12, wbdie('got inconsistent H+H'' @ %.3g',e3); e=4; end
     else
        e3=norm(q2-q0)/norm(q2);
        if e3>1E-12, wbdie('got inconsistent HAM/Ham @ %.3g',e3); e=5; end
     end

     if ~e
        if use_hconj
             wblog('ok.','got hermitian Ham @ %.3g [HAM @ %.3g]',e2,e1);
        else wblog('ok.','got hermitian HAM @ %.3g / %.3g',e1,e2); end
     end
  end
  tt=get_time(t0);

  p_=round(log(1E-14)/log(tau*dE));
  if isempty(pwr), pwr=min(4,p_);
     if pwr>4, wblog(' * ','expanding RHO up to power p=%d',pwr); end
  elseif pwr>min(4,p_)
     wblog(' * ','expanding RHO up to power p=%d / %d',pwr,p_);
  end

  n=max(6,pwr);
  for p=n:-1:1
     fac(p+1,:)= [ (-tau)^p/factorial(p), p ];
  end
  fac(1,:)=[1  0];

% add 3rd column to estimate numerical range of H^p (info only)
  fac(1:3,3)=[Dtot, E1, E2]; % assuming (Ham^0 = Id)
  for p=3:n, fac(p+1,3)=MPO_trace(Ham,p,E3); end

  tt(end+1)=get_time(t0);

% HAM.mpo may not be hermitian => need Ham here
  [Rho,Iout,IL]=MPO_add(Ham, fac(1:pwr+1,1:2), E3,osw);

  tt(end+1)=get_time(t0);

  add2struct(Iout,fac,pwr,'timing=tt',dloc,Dtot,Eref,tau); % e0
  % Iout already contains Nkeep, rtol, etc.

  Iout.Z =MPO_trace(Rho,E3); % normalization = partition function
  Iout.Z2=MPO_norm2(Rho);  % partition function at 2*tau
  % since Rho ~ Id => Z, Z2 ~ Dtot

  if kflag, wbstop; end

end

% -------------------------------------------------------------------- %
% see also tensor/MPO_add.m

function S=get_time(t0)

   S=dbstack();
   S=struct('line',S(2).line,'time',toc(t0));

end

% -------------------------------------------------------------------- %

