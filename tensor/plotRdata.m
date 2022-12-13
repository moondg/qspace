function [Iout,ES]=plotRdata(RR,varargin)
% function Iout=plotRdata(RR [,opts])
% Options
% 
%   'kQ',k  at which k to extract q-labels
%           (this also decides the order of qlabels)
%   'kQ',{k,nmax} show the lowest nmax qlabels (default: 5)
% 
% Examples
% 
%    plotRdata(RR,'k2x',{[-2 6],'J_{2}'},'dk',3)
% 
%    o={'k2x',{[ min(J(:,2)), max(J(:,2)) ], 'J_{\perp}'},'kQ',{10,8},'-x3'};
%    plotRdata(RR(1:2:end,1),o{:},'HSS',HAM.info.HSS{end})
%
%    plotRdata(RR,'k2x',{TT,'T','--log'})
% 
% Wb,Dec13,16

% adapted from @Hamilton1D/plotRdata.m 

  getopt('init',varargin);
     k2x    = getopt('k2x',{}); % { x(k), 'xlabel' }
     semin  = getopt('semin',[]);
     yl12   = getopt('yl12',[0 8]);
     kQ     = getopt('kQ',[]);
     ksplit = getopt('dk',[]);
     xstr   = getopt('xstr','');
     dyl    = getopt('dyl',[]);
     HSS    = getopt('HSS',[]);
     istr   = getopt('istr','');
     xlb    = getopt('xlb','bond index k');
     vflag  = getopt('-v');
     x3flag = getopt('-x3');

     if getopt('-gca'), ah=gca;
     else ah=getopt('ah',[]); end

  getopt('check_error');

  if numel(find(size(RR)>1))>1
     error('Wb:ERR','\n   ERR invalid usage (got matrix RR !?)');
  end

  if numel(yl12)==1 && yl12>0
     yl12=[0 yl12];
  end

  L=numel(RR); dd=nan(L,2); ES=QSpace;

  if ~isempty(k2x)
     if ~iscell(k2x) || numel(k2x)<2 || ~ischar(k2x{2})
        error('Wb:ERR','\n   ERR invalid usage (k2x)');
     end

     s=k2x{2};
     if numel(k2x{1})==2
        if numel(k2x)>2, k1L=k2x{3}; else k1L=[1 L]; end
        k2x={'lin',k1L,k2x{1},s};
     end
  end

if isempty(ah), pflag=1;
ah=smaxis(2,2+x3flag,'landscape','tag',mfilename); % ,'dy',0.10
  mvaxis(ah(2,:),[0 0.025]);
  if isempty(ksplit), ksplit=1; end
addt2fig Wb
else
  pflag=0;
  if isempty(ksplit), ksplit=numel(ah); end
end

% NB! k<kx vs. k>kx only flips q-labels with dual labels
% but leaves the eigenspectrum unchanged
% ==> for bond entanglement spectra for global singlets
%     q -> q_bar leads to same entanglement spectrum *by value*
%     while q and q_bar may very well have different 
%     spectra themselves when considering k<kx or k>kx individually!
% Wb,Mar22,15
% kx=L/2; => can lead to discontinuities in entanglement flow diagram!
  kx=0;

  nah=numel(ah);

if pflag
   setax(ah(2,1)); m=max(1,round(L/100));
end

  if ~isempty(HSS), rflag=0;
  else rflag=1; end

  for k=1:L
     if isdiag(RR(k))<=1
        [~,I]=eigQS(RR(k)); RR(k)=I.EK;
     end
  end

  for k=L:-1:1, if isempty(RR(k)), continue; end
     d=getDimQS(RR(k)); dd(k,:)=d(:,end);

     [se(k),rr]=SEntropy(QSpace(RR(k)),1E-8);

     ES(k)=RR(k); q=ES(k).data;
     for j=1:numel(q)
        q{j}=flip(-log10(q{j}(find(q{j}>1E-14))));
     end
     ES(k).data=q;

     if nah==1, continue; end

     if rflag==1 || x3flag
        n=sum(rr);
        if abs(n-1)>1E-8, wblog('WRN','got sum(rr)=1%+.3g !?',n-1); end
        j=find(rr>1E-16);
        wr(k)=numel(j)/sum(1./rr(j));
     elseif rflag==2
        if k<12 || k>L-12 || m==1 || mod(k,m)==1
           h(k)=semilogy(rr,'.-'); hold on
        end
     end
  end

  Iout=add2struct('-',se,ah);

  if ksplit, ok={'dk',ksplit}; else ok={'-1'}; end

  if pflag==0
     oe={'yl',yl12,'ah',ah,'-ldisp',ok{:}};

     [Iout.ES,Iout.hh,~,IEE]=getEEdata(ES,oe{:});
        if ~isempty(k2x), 
           for h=reshape(IEE.ah,1,[])
              setax(h), map_xdata(k2x{:});
           end
        end
        legdisp -erase
     return
  end

  if rflag==0
     if ~isnumeric(HSS)
        if ~isstruct(HSS) || numel(HSS)~=1
           error('Wb:ERR','\n   ERR invalid HSS data !?'); end
        q=HSS;
        HSS=q.HSS{end}; s=size(HSS);
        H=sparse(q.HH(:,1),q.HH(:,3),q.HH(:,end),s(1),s(2));
     else H=[];
     end

     q=diag(HSS,1); gotlad=(norm(q(2:2:end))==0);
     J=diag(H,1); i=find(J); q(i)=q(i)./J(i);

     lo={'Color',getcolor(1)};
     if gotlad
          h=plot(q(1:2:end),'.-',lo{:},'Disp','rung');
     else h=plot(q,'.-',lo{:},'Disp','<S{\cdot}S>_{k,k+1}');
     end
     blurl(h,'--line');
     hold on

     q=diag(HSS,2);
     J=diag(H,2); i=find(J); q(i)=q(i)./J(i);

     if gotlad
        for i=1:2
           h=plot(q(i:2:end),'.-','Color',getcolor(4-i));
           if i==1, set(h,'LineSt','--');
           else set(h,'Disp','leg');
           end
           blurl(h,'--line');
        end
     else 
        h=plot(q,'.-','Color',getcolor(2),'Disp','<S{\cdot}S>_{k,k+2}');
        blurl(h,'--line');
     end

     label(xlb,'<S_i{\cdot}S_{i+x}>');
     legdisp
  end

  if x3flag && ~rflag, setax(ah(2,3)); end
  if x3flag || rflag==1
     semilogy(wr,'.-');
     label(xlb,'smallest eigenvalues in \rho_k');
     Iout.wr=wr;
  elseif rflag==2
     hx=plot(nan,nan,'-');
     h(find(h==0))=hx;
     ylim([1E-16,1]); xtight(1.1,'view','x1',0); 
     recolor
  end

setax(ah(1,1));

  kk=1:length(se);

  h=plot(kk,se,'.-'); blurl(h,'--line'); hold on
  xtight(1.05,'x1',0);

  label(xlb,'block entropy S(k)');

if ~isempty(semin), kse_=semin;

  dse=diff(se(1:2:end));
  k=2*find(dse(1:end-1).*dse(2:end)<0);
  q=abs(k-kse_); i=find(q==min(q)); kse=k(i(1));

  i=max(1,kse-2):min(L,kse+2); m=min(5,length(i))-1;
  [p,s,mu]=polyfit(kk(i),se(i),m);
  x=linspace(i(1),i(end),1E4);

  y=polyval(p,x,s,mu); j=find(y==min(y),1); kmin=x(j);
  plot(x,y,'r');

  p1=p(1:m).*(m:-1:1);
  z=roots(p1); xmin=real(z(find(abs(imag(z))<1E-12)));
  if numel(xmin)>1
     q=abs(xmin-y(j)); xmin=xmin(find(q==min(q)));
  end

  kmin=xmin*mu(2)+mu(1);
  semin=polyval(p,kmin,s,mu); q=exp(semin);

  text(x(j),y(j),...
    sprintf('\nk_{min}=%.4g\n\\Rightarrow S=log(%.8g)',kmin,q),...
    'Color','r','VerticalAl','top','HorizontalAl','center','FontSize',8)
  xmark(kmin,'r:');
  if abs((q-round(q))/q)<1e-3, ymark(log(round(q)),'r:'); end

end

  setuser(gcf,'amin',[]); xd=[];

  if iscell(xstr) && numel(xstr)>1
     iso=xstr(2:end); if numel(iso)<2
        if isequal(iso{1},'SW')
           iso={iso{:},'scale',[1.6 0.7],'dx',[0.12 0.08]};
        end
     end
     xstr=xstr{1};
  else
     iso={'NW','scale',[1.6 0.7],'dx',[0.12 -.03]};
  end

  if ~isempty(xd)
     pk2x=polyfit(1:length(xd)-2, xd(2:end-1)', 3);

     x=0:length(xd); hold on
     h=plot(x,polyval(pk2x,x),'r--');

     i=find(abs(pk2x)<1E-12); pk2x(i)=0;

     s=poly2str2(pk2x,'x');
     title(s);

  else pk2x=[]; end

  setuser(gcf,'xd',xd);
  setuser(gcf,'xstr',xstr);
  setuser(gcf,'pk2x',pk2x);

  if ksplit>1
     setax(ah(1,1));
     splitline(ksplit); sms(4);
  end

  oe={'yl',yl12,'ah',ok{:}};

  [Iout.ES,Iout.hh]=getEEdata(ES,oe{:},ah(1,2));

  label(xlb,'entanglement spectrum (using log_{10})','entanglement flow');

  [EE,hh]=getEEdata(ES,oe{:},ah(2,2),'-rel');

  ylabel('[same as (c) above, yet with E_0^{ES} subtracted]');

  xtight('x1',0);

  if isset('kmin')
     xmark(kmin,'r:');
     setuser(gcf,'kmin',kmin);
  else
     setuser(gcf,'kmin',[]);
  end

setax(ah(1,2));

  if isempty(kQ), kQ=xlim; kQ=kQ(1); end
  if iscell(kQ), nmax=kQ{2}; kQ=kQ{1}; else nmax=4; end

  lh=findall(gca,'type','line');
  nh=numel(lh); yh=zeros(nh,1); qh=[];

  for i=1:nh, h=lh(i);
     q=getuser(h,'Q'); if isempty(q), continue; else qh(i,:)=q; end
     y=get(h,'YData'); j=find(~isnan(y));
     x=get(h,'XData'); j=j(findmel(x(j),kQ,-1));
     yh(i)=mean(y(j));
  end

  [qh,Iq,Dq]=uniquerows(qh); nq=size(qh,1); is=zeros(1,nq);
  for i=1:nq
     x=yh(Iq{i});
     is(i)=Iq{i}(find(x==min(x),1));
  end

  yh=yh(is); lh=lh(is);
  [yh,is]=sort(yh); ih=is(1:min(end,nmax));
  hl=lh(ih);

  for i=1:numel(hl), h=hl(i); u=get(h,'UserData');
     if isfield(u,'sym') && isfield(u,'Q')
        if numel(u.Q>1)
           if regexp(u.sym,',')
                s=sprintf('[%s ]',sprintf(' %g',u.Q));
           else s=sprintf('(%s)',sprintf('%g',u.Q)); end
        else s=sprintf('%g',u.Q);
        end
        set(h,'Disp',s);
     end
  end

  [~,i1,i2]=intersect(hl,get(gca,'Children')); [i1,is]=sort(i1); i2=i2(is);
  [~,is]=sort(i2); p(is)=1:length(is);
  legdisp('-erase','pick',fliplr(numel(p)+1-p));

if x3flag
   ah(1,3)=repax(ah(1,3),ah(1,2));
   h=findall(gca,'type','line'); set(h,'Marker','o'); sms(h,4);
   ylabel('[zoom into (c) with markers]');

   h=hl(1);
   y=get(h,'YData'); j=find(~isnan(y));
   x=get(h,'XData'); xj=x(j); dx=mean(diff(xj)); j=j(findmel(xj,kQ,-1)); j=j(1);

   xlim(x(j)+[-min(x(j),3*dx) 3*dx]); ylim([0 4]);
   if dx==round(dx)
      set(gca,'XTick',x(j)+(-5:5)*dx); grid on
   end
end

  q=max(dd);
  if diff(q)
       header(['%M :: L=%g, D^\ast \leq %g multiplets (%g states)'],L,q);
  else header(['%M :: L=%g, D \leq %g states'],L,q(1)); end
  setuser(gcf,'ah',ah);

  for i=1:numel(ah), poslabel(ah(i),i,'dx',[0 -0.02],'FontSize',13); end

  if ~isempty(k2x)
     for ih=1:4
        setax(ah(ih))
        map_xdata(k2x{:});
     end
  end

  if ~nargout, clear Iout; end

end

% -------------------------------------------------------------------- %

