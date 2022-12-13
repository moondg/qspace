function [b,a,I] = oliveira(Gamma,Lambda,z,varargin)
% function [beta,alpha,I] = oliveira(Gamma,Lambda,z [,opts])
%
%    Logarithmic discretization using Oliveira z-shift (PRB '90)
%
% Wb,May16,06

% -------------------------------------------------------------------- %
% for testing purposes [args directly handed over to oliveira()]:
%  '-pos'  just positive or
%  '-neg'  negative energies only (for testing purposes)
%  Wb,Mar10,16
% -------------------------------------------------------------------- %

  N=0; if nargout>2, I=struct; end
  getopt ('INIT', varargin);
     if getopt({'nolog','-q'}), qflag=1;
     elseif getopt('-Q'), qflag=2; else qflag=0; end
     N      = getopt('N',N  );
     xi     = getopt('xi'   );
     ALambda= getopt('-AL'  );
     Nmax   = getopt('Nmax',[]);

     if getopt('-pos'),     pnflag=+1;
     elseif getopt('-neg'), pnflag=-1;
     else pnflag=0; end

  getopt('check_error');

  if numel(N)~=1 || isinf(N) || isnan(N)
     error('Wb:ERR','\n   ERR invalid N (%g)',N); end
  if abs(N-round(N))>0
     q=N; N=ceil(N);
     wblog('WRN','got N=%g -> N=%g',q,N);
  end

  if isempty(Nmax)
     betaMin=1E-16;
     Nmax = ceil( 2* abs(log(betaMin) / log(Lambda)) );
  else
     betaMin=0;
     if Nmax<0, Nmax = 2*N; end
  end

  if Gamma<0, error('Wb:ERR','invalid Gamma (%g)',Gamma); end

  if z==0 && ~ALambda && ~pnflag
     if qflag==1, wblog('<i>',... 
     'plain Wilson (z=%g; Lambda=%g; Gamma=%g)',z,Lambda,Gamma); end

     ii=0:N-3; Delta=1/Lambda;
     xx=(1-Delta.^(ii+1))./sqrt((1-Delta.^(2*ii+1)).*(1-Delta.^(2*ii+3)));

     if xi
          b = [ sqrt(2*Gamma/pi), (1+Delta)/2                  * xx ];
     else b = [ sqrt(2*Gamma/pi), (1+Delta)/2*(Delta.^(ii/2)) .* xx ]; end
     a = zeros(1,length(b)+1);

     I.pnflag=pnflag;
     return

  elseif N>Nmax
     wblog('WRN','N too large; set N=%g -> %g (z=%g)', N, Nmax, z);
     N=Nmax;
  end

  z_=mod(z,1);
  z =mod(1-z,1);

  iz=0:Nmax; if z, iz=[0,iz+z]; end

  Edisc = (1/Lambda).^iz;
  n=length(Edisc)-1;

  if ALambda
   % NB! there appears to be systematic error / offset for Lambda<2
   % e.g. standard Friedel sum rule for SIAM well-obeyed for Lambda>=2,
   % but A(0)=~0.98 becomes < 1 for Lambda=1.7 (!!) the latter case
   % seems to work better with ALambda turned off! (systematic error
   % compensates overshoot?) -- Wb,Jun16,09

     d0=((1-1/Lambda)/log(Lambda)) * Edisc(1:end-1);
     if z==0
          d0(1)=      (1-1/Lambda  )/log(Lambda);
     else d0(1)=(1-z)+(1-1/Lambda^z)/log(Lambda);
     end

  else
     d0=0.5*(Edisc(1:end-1) + Edisc(2:end));
  end

  d1 = sqrt(abs(diff(Edisc)));

  hdd=sparse(diag(d0)); Z=sparse(n,n);

  if pnflag==0
     H0=[ 0     d1    d1
          d1'  +hdd   Z
          d1'   Z    -hdd ];
  elseif pnflag>0
       H0=[ 0,  d1; d1', +hdd ];
  else H0=[ 0,  d1; d1', -hdd ];
  end

  D=size(H0,1);
  niter=D;

  if ALambda
     s='ALambda; ';
     global param; param.ALambda='using shifted discrete energies';
  else s=''; end

  if qflag
     if qflag==1, wblog('<i>',...
     'Oliveira (%sz=%g; Lambda=%g; Gamma=%g)',s,z_,Lambda,Gamma); end
  else
     wblog('<i>','Oliveira (%sz=%g)',s,z_);
     wblog('Lambda = %g', Lambda);
     wblog('Gamma  = %g', Gamma );
     wblog('D      = %d -> %d', n, D);
  end

  U=zeros(D,1); U(1)=1;

  alpha=zeros(1,niter); beta=zeros(1,niter);

  for i=1:niter
     v=H0*U(:,i);
     alpha(i)=U(:,i)'*v;

     v=v-U*(U'*v); v=v-U*(U'*v);
     beta(i)=norm(v);

     if beta(i)<=betaMin
        if ~qflag, wblog('<i>',...
        'Lanczos converged at iteration %d/%d (%.3g).',i,D,beta(i)); end
        break
     elseif i>=N, break;
     else U(:,i+1)=v/beta(i); end
  end

  alpha=alpha(1:i); beta=beta(1:i-1); beta_=beta;

  beta(1)=beta(1)*sqrt(Gamma/pi);

  if xi
       a=alpha; b=beta .* [1, Lambda.^(((1:length(beta)-1)-1)/2)];
  else a=alpha; b=beta; end

  k=length(beta);
  br=abs(beta(3:k)./beta(2:k-1)); e=max(br);
  if e>=1
       kx=find(br==max(br),1);
  elseif beta(2)>beta(1), kx=2; else kx=1; end

  if kx>3
     s=sprintf(['got increasing couplings: ' ... 
       'ff(%g)=%.3g @ %.3g (z=%.3g)'],kx+1,beta(kx+1),br(kx),z);
     if e>2, error('Wb:ERR',['\n   ERR ' s]);
     else wblog(iff(kx<5,'WRN','ERR'),s); end
  end

  if nargout>2,
     I=add2struct('-',...
     Lambda,z,Gamma,Edisc,d0,alpha,beta,H0,U,kx,ALambda);

     if ALambda
        I.ALfac = Lambda^z*(Lambda-1)/log(Lambda);
     end
  end

  if qflag
     a=a(1:N); b=b(1:N-1);
     return
  end

% -------------------------------------------------------------------- %
ah=smaxis(2,2,'tag',mfilename);

  header('%M :: NRG coefficients (z=%g, \\Gamma=%g, \\Lambda=%g)',...
   z_, Gamma, Lambda)

setax(ah(1,1));

  semilogy(abs(alpha),'ko-','Disp','\alpha'); hold on
  semilogy(abs(beta), 'ro-','Disp','\beta'); sms(2)
  xtight(1.1); ytight(1.1);

  if ~pnflag
     h=plot(sqrt(Edisc),'b','Disp','E_{disc}^{1/2}');
  else
     h=plot(Edisc,'b','Disp','E_{disc}');
     postext('NE','NB! got ''%s'' flag',iff(pnflag>0,'-pos','-neg'));
  end
  blurl(h); mv2back(h);

  legdisp('Location','W');

setax(ah(1,2));

  if ~pnflag && abs(z)<1E-3
     beta0=zeros(size(beta));
     beta0(1)=sqrt(2*Gamma/pi);

     ii=2:length(beta);
     beta0(2:end) = ...
        (1+1/Lambda) ./ (2*Lambda.^((ii-2)/2)) .* ...
        (1-Lambda.^(-(ii-2)-1)) ./ (sqrt(1-Lambda.^(-2*(ii-2)-1)) .* ...
        sqrt(1-Lambda.^(-2*(ii-2)-3)));

     semilogy(abs(beta-beta0),'r'); hold on
     semilogy(abs(alpha),'k'); axis tight

     title('difference to Wilson''s \xi_{ i}')
  end

setax(ah(2,1));

  T=diag(alpha)+diag(beta_,1)+diag(beta_,-1);

  Q=U'*H0*U-T;

  mp(Q,'cmap','-gca'); 
  postext(0.65,0.9,'U''HU-T',{'BackgroundC','w'});

  e=norm(Q(:));
  if e>10*beta(end)
     wblog('WRN','got e=%.3g (%.3g) !?',e,beta(end)); 
  end

setax(ah(2,2));

  mp(U'*U-eye(size(U,2)),'cmap','-gca');
  h=postext(0.75,0.9,'U''U-E',{'BackgroundC','w'});

  a=a(1:N); b=b(1:N-1);

end

% -------------------------------------------------------------------- %

