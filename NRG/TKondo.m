function [TK,nK,Iout]=TKondo(varargin)
% function [TK,nK,Iout]=TKondo(varargin)
%
%    This routine takes param structure or 'var',value,... list
%    as input and returns Kondo temperature.
%    
%    If Lambda is specified with the input
%    nK gives the Wilson chain position that represents TK
%
% Examples
%
%    TK=TKondo(param,'epsd',0);  % further args overwrite values in param
%
% Wb,Mar01,05 ; Wb,Aug18,17

  global param g_TKondo

  if ~isempty(g_TKondo), TK=g_TKondo;
     if ~isnan(g_TKondo)
        wblog(' * ','using global TK=%.3g', g_TKondo);
        TK=check_value(TK);
     end
     return
  end

  getopt('init',varargin);
     vflag = getopt('-v');
     wTK   = getopt('wtk',[]);
  varargin=getopt('get_remaining'); narg=numel(varargin);

  if narg && isstruct(varargin{1})
       p=varargin{1}; varargin(1)=[]; narg=numel(varargin);
  else p=param; end

  ff={'Lambda','J','fJ','Jz','U','epsd','Gamma'};
  n=numel(ff); i1=2;

  if narg
     getopt('init',varargin);
     if nargout>1, i1=1; end
     for i=i1:n, f=ff{i};
        x=[]; if isfield(p,f), x=getfield(p,f); end
        x=getopt(f,x); p=setfield(p,f,x);
        if (i==2 || i==3) && ~isempty(x), break; end
     end
     getopt('check_error');
  else
     for i=i1:n
        if isfield(p,ff{i}), break; end
     end
  end

  if i==3, p.J=p.fJ;
  elseif i==2 && isfield(p,ff{3})
     wblog('WRN','got values for J and fJ (ignore)');
  end

  if  isfield(p,'epsd'), [TK,wTK,iTK]=TKondo_SIAM(p,wTK);
  elseif isfield(p,'J'), [TK,wTK,iTK]=TKondo_J(p,wTK);
  else wbdie('invalid usage (expecting Kondo or SIAM parameters)');
  end

  TK=check_value(TK);

  if nargout>1, nK=nan;
     if ~isempty(p.Lambda)
          nK=-2*log(TK)/log(param.Lambda);
     else wblog('WRN','got missing option / field ''Lambda''');
     end
     if vflag, wblogx(1,'<i> %s (nK=%.1f)',iTK,nK); end
  elseif vflag
     wblogx(1,'<i> %s',iTK);
  end

  if nargout>2
     Iout=add2struct(p,wTK,iTK);
  end

end

% -------------------------------------------------------------------- %

function TK=check_value(TK,p)

  if numel(TK)~=1
     if nargin>1, disp(p); end
     wblog('WRN','got TK=[%s] -> %g !?',TK,vec2str(TK(:)),nan);
     TK=nan;
  elseif isinf(TK) || TK<0
     if nargin>1, disp(p); end
     wblog('WRN','got TK=%g -> %g',TK,nan);
     TK=nan;
  elseif isnan(TK)
     if nargin>1, disp(p); end
     wblog('WRN','got TK=NaN');
  elseif TK>10
     if nargin>1, disp(p); end
     wblog('WRN','got TK=%g',TK);
  end

end

% -------------------------------------------------------------------- %
function varargout=deal_pfield(param,varargin)

   nargs=numel(varargin);
   varargout=repmat({nan},1,nargs);

   for i=1:nargs, v=varargin{i};
      if ~isfield(param,v)
         wblog('WRN','missing field or variable ''%s''',v);
      else
         q=getfield(param,v); n=numel(q);
         if isempty(q), TK=nan;
            wblog('WRN','got empty value for ''%s''',v);
         else
            if n>1, q=unique(q(:)); n=numel(q); end
            if n>1
               varargout{i}=mean(q); wblog('WRN',...
               'got %g values for ''%s'' (taking mean %g)',n,v,varargout{i});
            else varargout{i}=q; end
         end
      end
   end

end

% -------------------------------------------------------------------- %

function [TK,wTK,istr]=TKondo_J(p,wTK)

   istr='';

   J=p.J; n=numel(J); q=[ n==3, isfield(p,'Jz') ];
   if any(q)
      if q(2), Jz=p.Jz;
         if q(1) && J(3)~=Jz, disp({ J, Jz })
            wbdie('invalid usage [J(3) ~= Jz]'); 
         end
      else 
         if diff(J(1:2)), J, wbdie('invalid usage (Jx ~= Jy !?)'); end
         Jz=J(3);
      end
      if Jz==abs(J(1)) && Jz>0, J=Jz; n=1;
      else 
        Jx=abs(J(1));
        a2=Jz^2 - Jx^2;
        if Jz>Jx && Jx
           a=sqrt(a2);
           TK = ((Jz-a)/(Jz+a)) .^ (1/(2*a));
           wTK='zKondo';
        elseif a2<0
           a=sqrt(-a2);
           TK = exp( (atan(Jz/a) - pi/2)/a );
           wTK='xKondo';
        else
           TK=0;
           wTK='noKondo';
        end
        TK=TK*sqrt(norm([Jx,Jx,Jz])/sqrt(3));
        return
      end
   end

   if n~=1, wbdie('invalid usage (got %g values for J)',n); end

   if isnan(J) || J<0, TK=nan; return
   elseif J==0, TK=0; return; end

 % in Rosch et al. (PRB 2003) the matrix elements are J/2
 % same as in M.Sindel (thesis p. 36)
 % equation TK is offset from energy flow convergence regime
 % => take J->J/2 for definition of Kondo temperature
 % NB!  Rosch et al. (PRB 2003) footnote [33]
 %   TK =~ sqrt(J)*exp(-1/J)   - exact within 2-loop perturbation
 %   TK =~         exp(-1/J)   - exact within 1-loop perturbation

   if isfield(p,'wsys') && ~isempty(regexp(p.wsys,'\dCK'))
        wTK='Kondo2';
   else wTK='Kondo'; end

   switch wTK
      case 'Kondo'
       % using J = fJ = J_calc = 2 \nu JK = 2nJ
       % ==> TK \equiv sqrt(2nJ) exp(-1/2nJ) (!!)
       % equivalently: using (2D nuJ) Sd.S0 = (nuJ) Sd.tau0 as Kondo Hamiltonian
       % => J_NRG = nuJ (see setupKondo*.m)
       % => still: TK =~ sqrt(nuJ)*exp(-1/nuJ) // Wb,Sep08,13
         TK=sqrt(J)*exp(-1/J);
      case 'Kondo2'
         TK=exp(-2/J);
      otherwise, wTK, error('Wb:ERR','\n   ERR invalid switch'); 
   end

   if nargout>1
      istr=sprintf('J=%g => TK=%g',J,TK);
   end

end

% -------------------------------------------------------------------- %
% * Hewson, p.66
%   jk =abs(e*(e+u))/(g*u); TK = sqrt(1/jk)*exp(-pi*jk/2);
%   TKondo = sqrt(U*Gamma) * exp(pi*epsd*(epsd+U)/(U*Gamma));
% * Thesis, M. Sindel, p.29
%   NB! requires Gamma -> Gamma/2, otherwise TK too small!)
% * Haldane, PRL, 6, 416 (1978) - Delta \equiv myGamma
%   TK = sqrt(u*g)/2 * exp(-abs(pi*e*(e+u)/(u*g)));
% * for symmetric Anderson model (see Merker12 + Costi) // Wb,Jun20,12
%   TK = sqrt(u*g/2) exp(-pi(u/8g - g/2u)) \equiv 1/4chi0
%   TK = sqrt(u*g/2) * exp(pi*e*(e+u)/(2*u*g));
% -------------------------------------------------------------------- %

function [TK,wTK,istr]=TKondo_SIAM(p,wTK)

   [u,g,e]=deal_pfield(p,'U','Gamma','epsd'); istr='';

   if any(isnan([u,g,e])) || any([u,g]<0), TK=nan; return; end
   if g==0, TK=0; return; end
   if u==0, TK=g; return; end

   if isempty(wTK), wTK='SIAM'; end

 % -------------------------------------------------------- %
 % NB! the prefactor of sqrt(u*g/2) is somewhat problematic
 % as it goes to zero for U->0 (resonant-level model) (!?!)
 % Wb,Mar28,13
 % see also Filippone17, Eq.(22) (!) // Wb,Aug18,17
 % including higher order correction for SIAM!
 % -------------------------------------------------------- %
   switch wTK
     case 'bethe'

        a = pi*(-u/(8*g) + g/(2*u));
        x = (e+u/2)*sqrt(pi/(2*u*g));

        TK=min(1,sqrt(u*g/2))*exp(a+x^2);

     case 'SIAM', a=pi*e*(e+u)/(2*u*g);
       TK=min(0.575,sqrt(u*g/2))*exp(a);

     case 2, a=pi*e*(e+u)/(2*u*g); TK=sqrt(u*g/4)*exp(2*a);

     case 3, a=pi*e*(e+u)/(2*u*g); TK=sqrt(u*g/(2*3.00))*exp(1.050*a);

     otherwise wTK, error('Wb:ERR','invalid switch');
   end

   if TK>10, TK=10; end

   if nargout>1
      istr=sprintf('U=%g, Gamma=%g, epsd/U=%g => TK=%g',u,g,e/u,TK);
   end
end

% -------------------------------------------------------------------- %

