
% simple plot routine for the data produced by tst_Hamiltonian1D
% Wb,Aug28,15

  L=length(HAM.mpo);

  if ~isset('isw'), isw=HAM.info.isw; end
  if ~isset('nsw'), nsw=nan; end
  if ~isset('IL'), IL=Il; end

  if isemptystruct(Il(1,2))
     Il(1:end-1,2)=Il(2:end,2);
     Il=Il(1:end-1,:);

     IL(1:end-1,2,:)=IL(2:end,2,:);
     IL=IL(1:end-1,:,:);
  end

  e0=getdatafield(Il,'E0','-q');
  r2=getdatafield(Il,'r2','-q'); i1=ceil(L/3); i2=L-i1+1;
  q=reshape(e0(i1:i2,:)',[],1);
  r2=mean(reshape(r2(i1:i2,:)',[],1));

  istr=sprintf('sweep %g/%g : <e_0>=%s (@ %.3g)',...
    isw,nsw,num2str2(mean(q),std(q)),r2);

ah=smaxis(2,2,'tag',mfilename,'landscape');
header('%M :: %s\newline%s',regexprep(HAM.info.istr,'\n',' '),istr); addt2fig Wb
mat_to_headerSE

  if isset('NPsi')
     d=getDimQS(HAM(1).Psi); % Wb,Feb22,17
     if size(d,2)==4
          npsi=d(:,end)';
     else npsi=0; end

     if npsi(1)~=NPsi
        wblog('WRN','got NPsi = %g / %g !?',NPsi,npsi(1));
     end
  else npsi=0; end

setax(ah(1,1))

  se=getdatafield(IL,'se','-q');
  hh=plot(exp(se(:,:)));

  if npsi
     if numel(npsi)>1
          s=sprintf('N^\\ast_\\Psi= %g (%g)',npsi);
     else s=sprintf('N_\\Psi= %g',npsi);
     end
     postext('N',s,{'Color','r'});
  end

% NB! for NPsi>1, got L/R asymmetric SEntropy data! // Wb,Feb22,17
  set(hh(2:2:end),'LineSt','--');

% for i=1:2:numel(hh), h=hh(i);
%    x=get(h,'XData'); set(h,'XData',x+1);
% end

  xtight(1.01); xl=get(gca,'XLim');
% scalecolor(hh,'alpha',0.99);

  smax=max(se(:));
  ymark(exp(smax),'r--',...
    'istr',{sprintf('S=%.5g',smax),0.02,'BackgroundC','w','FontSize',8})

  label('site k','exp( entanglement entropy )');

setax(ah(1,2))

  e0=getdatafield(IL,'E0','-q'); eg=min(e0(:));
  hh=semilogy(e0(:,:)-eg);

% NB! for NPsi>1, got L/R asymmetric SEntropy -> E0 data! // Wb,Feb22,17
  set(hh(2:2:end),'LineSt','--');

% for i=1:2:numel(hh), h=hh(i);
%    x=get(h,'XData'); set(h,'XData',x+1);
% end

  label('site k','ground state convergence');
  title2('e_0^{ref}=%.8g',eg);
  xlim(xl);

setax(ah(2,1))

  r2=getdatafield(IL,'r2','-q');
  hh=semilogy(r2(:,:));

% NB! for NPsi>1, got L/R asymmetric SEntropy -> rdisc data! // Wb,Feb22,17
  set(hh(2:2:end),'LineSt','--');

% for i=1:2:numel(hh), h=hh(i);
%    x=get(h,'XData'); set(h,'XData',x+1);
% end

  label('site k','discarded weight \delta\rho');
  xlim(xl);

  y=ytight(1.1,'y2',1); if y(1)<1E-15, y(1)=1E-15; ylim(y); end

setax(ah(2,2))

  if isfield(IL,'Dk') && ~isempty(IL) && size(IL(2).Dk,1)>1
    nk=getdatafield(IL,'Dk',1,'-q');
    Dk=getdatafield(IL,'Dk',size(IL(2).Dk,1),'-q');
  else
    nk=getdatafield(IL,'NK','-q');
    Dk=[];
  end

  if ~isempty(Dk)
    h2=semilogy(Dk(:,:)); blurl(h2,'LineW',1); hold on
    Dmax=max(Dk(:));
    title2('N^\ast_{kept} \leq %g (%s)',max(nk(:)),int2str2(Dmax));
  % ymark(Dmax,'k:');
  else
    title2('N_{kept} \leq %g',max(nk(:)));
  end
 
  hh=semilogy(nk(:,:));

% NB! for NPsi>1, got L/R asymmetric SEntropy -> Dk data! // Wb,Feb22,17
  set(hh(2:2:end),'LineSt','--');

% for i=1:2:numel(hh), h=hh(i);
%    x=get(h,'XData'); set(h,'XData',x+1);
% end

  label('site k','D{^\ast} (= kept number multiplets)');
  xlim(xl);
 
  set(gca,'YTick',unique(NKEEP(1:size(IL,3))));
  ytight(1.1);

setax(ah(1,2));
ih12=inset('SE','scale',[1.5 1.0],'dx',[0 0.1]);

  ns=size(se,3); s1=zeros(ns,1); r1=s1; e1=s1; n1=s1;
  for k=1:ns
     s1(k)=max(reshape(se(:,:,k),[],1));
     r1(k)=max(reshape(r2(:,:,k),[],1));
     n1(k)=max(reshape(nk(:,:,k),[],1));
     e1(k)=min(reshape(e0(:,:,k),[],1));
  end

  if 1
     x=1./n1; y=e1-eg;
     h=loglog(x,y,'o-'); hold on
     label('1/N','e_0 - e_0^{ref}');
  else
     h=semilogx(r1,e1,'o-'); hold on
     label('r^2','e_0');
  end

  xtight(1.05,'x1',0);
  mvlabel('x',0.2)

setax(ah(2,2));
ih22=inset('SE','scale',[1.5 1.0],'dx',[0 0.15]);

  x=1./n1; y=s1;
  h=plot(x,y,'o-'); hold on

  label('1/N','max(S)');
  xtight(1.05,'x1',0);

  i=find(n1>60);
  if numel(i)>=3
   % i=ns-3:ns;
     [p,mu]=polyfit(x(i),y(i),2);

     xx=linspace(0,max(x(i)));

     set(h,'Marker','none');
     h=plot(xx,polyval(p,xx,mu),'r');
     h=plot(x,y,'o'); sms(h,4);

     title2('S_{\infty}=%.4g',polyval(p,0,mu));
  end

  togglelogx

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% Wb,Jun22,17

if isvar('Jp1') && isvar('Jp2') && Jp1<Jp2, setax(ah(1,1));
   setax(ah(1,1));
   map_xdata(gca,'lin',[1 L],[Jp1 Jp2],'xlb','J_{\perp}');
%  h=get(gca,'XLabel');
%  set(h,'String',sprintf('%s \\rightarrow %s',get(h,'String'),x1opts{2}));
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

