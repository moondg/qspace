function [g0,a0,g1]=getConductance0(om,a0,varargin)
% function [g0,a0]=getConductance0(om,a0 [,T,gfac,...])
%
%    Get conductance by folding discrete a0 data with
%    derivative of Fermi distribution (-df/dw).
%
% Options
%
%   '-f'  got functional data (rather than discrete data set)
%
% Return values
%
%    g0 regular conductance (Meir-Wingreen)
%    g1 only includes discrete contributions strictly within |omega|<=1
%
% outsourced from MEX//rnrg_gg.m ; tags: conductivity
% Wb,Jul08,13

  if nargin<2 || ~isvector(om)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  e=0; fflag=0; T=[]; param=[]; Gamma=[]; gfac=[]; 

  while ~isempty(varargin) && isstruct(varargin{1}), q=0;
     if isfield(varargin{1},'T') && isempty(T), q=q+1;
        T=varargin{1}.T; end
     if isfield(varargin{1},'param') && isempty(param), q=q+1;
        param=varargin{1}.param; end
     if isfield(varargin{1},'Gamma') && isempty(param), q=q+1;
        param=varargin{1}; end
     if q, varargin=varargin(2:end); else break; end
  end

  while ~isempty(varargin) && isnumber(varargin{1})
     if isempty(T), T=varargin{1};
     elseif isempty(gfac), gfac=varargin{1};
     else break; end
     varargin=varargin(2:end);
  end

  if isempty(T)
     getbase('Idma');
     if isfield(Idma,'T'), T=Idma.T; else e=e+1; end
  end

  if isempty(gfac)
     if isempty(param),
        param=getuser(0,'param');
        if isempty(param)
           getbase('Inrg');
           if isfield(Inrg,'param'), param=Inrg.param; end
        end
     end
     if isfield(param,'Gamma'), gfac=pi*param.Gamma;
     else e=e+10; end
  end

  if e, error('Wb:ERR','\n   ERR invalid usage (e=%g)',e); end

  if numel(varargin)
     getopt('init',varargin);
        fflag=getopt('-f');
     getopt('check_error');
  end

  if size(om,2)~=1, om=om'; end

  beta=1/T;
  fw=(beta/2)./(1+cosh(beta*om));

  if fflag
     dw=[ om(2)-om(1); 0.5*(om(3:end)-om(1:end-2)); om(end)-om(end-1) ];
     fw=fw.*dw;
  end

  g0=gfac*(fw'*a0);

  if nargout>1
     if nargout>2
        i=find(abs(om)<=1);
        g1=gfac*(fw(i)'*a0(i,:));
     end

     dE=T/10; fw=exp(-(om/dE).^2) / (sqrt(pi)*dE);
     a0=fw'*a0;
  end

end

