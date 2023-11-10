function [fo,I]=openWilson(Gamma,Lambda,N,varargin)
% function [fo,I]=openWilson(opts)
% Wb,Nov10,16

  getopt('init',varargin);
     npts =getopt('npts',2048);
     z    =getopt('z',0);
     sigma=getopt('sigma',1.28);
     pflag=getopt('--plot');
     kflag=getopt('--kb');
     vflag=getopt('-v');
  getopt('check_error');

  if nargin<3 || nargout>2 || ...
    ~isnumber(Gamma )   || Gamma <=0 || ...
    ~isnumber(Lambda)   || Lambda<=1 || ...
    ~isnumber(N) || N<3 || N~=round(N) || ~isnumber(z)

     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  z=mod(z,1);
  I=add2struct('-',Gamma,Lambda,N,z,npts);

  wblog('<i>','openWilson (z=%g, Lambda=%g, Gamma=%g)',z,Lambda,Gamma); 

  Ex=[1, Lambda.^(-(0.010+2*z:N+4)/2)];
  Ex=Ex(1:N+1);

  ek=linspace(-Ex(1),Ex(1),2*npts)';

  u2=repmat(1/2,size(ek));

  EU=repmat({{ek,u2}},1,2);

  if kflag || pflag, vflag=1; end

  for k=1:N
     if vflag, fprintf(1,'\r  %4g/%g ...\r',k,N); end
     de=diff2(ek,'len');

     oz=complex(ek,sigma*de);

     G=zeros(size(u2));
     for i=1:numel(u2)
        G(i)=de'*(u2./(oz(i)-ek));
     end

     aa=-(1/pi)*imag(G);
     nrm(k,1:2)=[ intxy(ek,u2), intxy(ek,aa) ];

     alpha(k)=intxy(ek,ek.*aa);

     Delta=ek-alpha(k)-1./G;
     Gamma=-(1/pi)*imag(Delta);

     EU{k+1,1}={ek,Gamma/intxy(ek,Gamma)};

     i=find(ek> Ex(k+1)); t2(1)=intxy(ek(i),Gamma(i));
     i=find(ek<-Ex(k+1)); t2(2)=intxy(ek(i),Gamma(i)); gamma(k)=sqrt(sum(t2));
     i=find(abs(ek)<=Ex(k+1)); t2=intxy(ek(i),Gamma(i)); beta(k)=sqrt(t2);

     u2=Gamma/t2;

     e2=linspace(-Ex(k+1),Ex(k+1),2*npts)';
     u2=spline(ek,u2,e2); ek=e2;

     e=sum(u2(find(u2<0))); if e
        wblog('WRN','k=%g got negative u2 @ %.3g !?',k,e); end

     nrm(k,3)=intxy(ek,u2);
     EU{k+1,2}={ek,u2};

  end
  if vflag, fprintf(1,'\r60s\r',''); end

  I.alpha=alpha;
  I.beta =beta;
  I.gamma=gamma;

  if nargout
     fo=[sqrt(2*I.Gamma/pi), beta];
  end

  if ~pflag, return; end

ah=smaxis(2,2,'tag',mfilename,'dx',0.12,'dy',0.12); addt2fig wb
header('%M :: \Gamma=%g, \Lambda=%g, N=%g, npts=%g, z=%g',I.Gamma,Lambda,N,npts,z);

  Escale=Lambda.^(0:-.5:-(N-1)/2);

setax(ah(1,1))

  [ff,q,Ic]=getNRGcoupling(I.Gamma,Lambda,N+2,'-x','z',z);
  [fw,q,Ic]=getNRGcoupling(I.Gamma,Lambda,N+2,'-x','z',z,'-w');

  ff=ff(2:end);
  fw=fw(2:end);

  plot(beta,'.-','Disp','\beta'); hold on
  plot(ff,'.-','Disp','Zitko');
  plot(fw,'.-','Disp','Wilson');
  plot(abs(alpha),'.-','Disp',sprintf('|\\alpha|\\leq%.3g',max(abs(alpha))));
  plot(gamma,'.-','Disp','\gamma');

  set(gca,'YScale','log');

  legdisp('Location','W');

  label('Wilson shell n','energy [D]');

setax(ah(1,2))

  h1=plot(beta./Escale,'.-','Disp','\beta'); hold on
  plot(ff./Escale,'.-','Disp','Zitko');
  plot(fw./Escale,'.-','Disp','Wilson');

  h0=plot(alpha*sqrt(Lambda)./Escale,'.-','Disp','\alpha');
  h2=plot(gamma*sqrt(Lambda)./Escale,'.-','Disp','\gamma');

  if z==0, k=1; elseif z==0.5, k=2; else k=[]; end
  if Lambda==2
     q=load('openWilson_SeungSup_Lambda=2');
  elseif Lambda==4
     q=load('openWilson_SeungSup_Lambda=4');
  else q=[]; end

  if ~isempty(q) && ~isempty(k)
     plot(q.ff(2:end,k) ./ Escale(1:length(q.ff)-1)','--',...
     'LineW',2,'Color',get(h1,'Color'),'Disp','\beta Seung-Sup');
  end

  legdisp('Location','SE','-erase');
  label('Wilson shell n','energy / \omega_n');

setax(ah(2,1))

  for k=1:size(EU,1)
     h=plot(EU{k,2}{1}, EU{k,2}{2}); hold on
     if k==2, set(h,'Color',[0 .5 0]);
     elseif k==3, set(h,'Color','r');
     end
  end

  label('|\omega|','u^2_n(\omega)');

  togglelogx; togglelogy

  h=xmark(Ex,'k'); set(h,'Color',[.9 .9 .9]);
  q=diff(log(Ex)); s=sprintf(', %.3g',q(1:4)/q(end));
  title(sprintf('{\\delta}k = [ %s ...]',s(3:end)));

setax(ah(2,2))

  for k=1:size(EU,1)
     h=plot(EU{k,1}{1}/Ex(k), EU{k,1}{2}*Ex(k)); hold on
     if k==2, set(h,'Color',[0 .5 0]);
     elseif k==3, set(h,'Color','r');
     end
  end
  xtight(1.1); ytight(1.1,'y1',0);
  xmark([-1 1],'k--');

  q=nrm(:); q=[mean(q),std(q)]; q(2)=q(2)/q(1);
  s=sprintf('normalization = %.4g @ %s',q(1),num2tex(q(2),'%.2E'));
  xlabel('\omega / \omega_n');
  ylabel('u^2_n(\omega) \cdot \omega_n [same as left; rescaled]');
  title(s);

  e=norm([q(1)-1,q(2)]);
  if e>0.05, wblog('ERR','got %s !?',s); end

  if kflag, keyboard; end

end

