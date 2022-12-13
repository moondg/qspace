
% simple plot routine for the data produced by tst_Hamiltonian1D
% Wb,Aug28,15

  L=length(HAM.mpo);
  param=HAM.info.param;

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
header('%M :: %s\n%s',regexprep(HAM.info.istr,'\n',' '),istr); addt2fig wb
mat_to_headerSE

  if isset('NPsi')
     npsi=getfield2(Inrg,'NPsi','--def',[]);
     if ~isempty(npsi), NPsi=npsi;
     else
        d=getDimQS(HAM(1).Psi);
        if size(d,2)==4
             npsi=d(:,end)';
        else npsi=0; end

        if npsi(1)~=NPsi
           wblog('WRN','got NPsi = %g / %g !?',NPsi,npsi(1));
        end
     end
  else npsi=0; end

setax(ah(1,1))

  se=permute(getdatafield(IL,'se','-q'),[1 3 2]);
  smax=max(reshape(se(:,end,:),[],1));

  h1=plot(exp(se(:,:,1))); set(gca,'ColorOrderIndex',1); hold on
  h2=plot(exp(se(:,:,2)));

  set(h2(end),'Color','k','LineW',1.5);

  if npsi
     if numel(npsi)>1
          s=sprintf('N^\\ast_\\Psi= %g (%g)',npsi);
     else s=sprintf('N_\\Psi= %g',npsi);
     end
     postext('N',s,{'Color','r'});
  end

  blurl(h1,'LineW',1,'cfac',0.5); % ,'LineSt',':');

  xtight(1.01); xl=get(gca,'XLim');

  o={'FontSize',8,'VerticalAl','bottom'}; % ,'BackgroundC','w'
  ymark(exp(smax),'r--','istr',{sprintf('S=%.5g',smax),0.02,o{:}});

  label('site k','exp( entanglement entropy )');

setax(ah(1,2))

  e0=getdatafield(IL,'E0','-q');
  eg=min(e0(:))-std(reshape(e0(:,:,end),[],1));

  e0=permute(e0,[1 3 2]);
  h1=semilogy(e0(:,:,1)-eg); set(gca,'ColorOrderIndex',1); hold on
  h2=semilogy(e0(:,:,2)-eg);

  blurl(h1,'LineW',1,'cfac',0.5); % ,'LineSt',':');

  label('site k','ground state convergence');
  title2('e_0^{ref}=%.8g',eg);
  ytight(1.1); xlim(xl);

setax(ah(2,1))

  r2=permute(getdatafield(IL,'r2','-q'),[1 3 2]);
  h1=semilogy(r2(:,:,1)); set(gca,'ColorOrderIndex',1); hold on
  h2=semilogy(r2(:,:,2)); set(h2(end),'Color','k','LineW',1.5);

  blurl(h1,'LineW',1,'cfac',0.5); % ,'LineSt',':');

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
  nk=permute(nk,[1 3 2]);

  if ~isempty(Dk)
    h2=semilogy(Dk(:,:)); blurl(h2,'LineW',1); hold on
    Dmax=max(Dk(:));
    title2('N^\ast_{kept} \leq %g (%s)',max(nk(:)),int2str2(Dmax));
  else
    title2('N_{kept} \leq %g',max(nk(:)));
  end

  h1=semilogy(nk(:,:,1)); set(gca,'ColorOrderIndex',1); hold on
  h2=semilogy(nk(:,:,2));

  blurl(h1,'LineW',1,'cfac',0.5); % ,'LineSt',':');

  label('site k','D{^\ast} (= kept number multiplets)');
  xlim(xl);

  set(gca,'YTick',unique(NKEEP(1:min(numel(NKEEP),size(IL,3)))));
  ytight(1.1);

setax(ah(1,2));
ih12=inset({'SE',[0 0.15]},'scale',[1.2 1.2]);

  ns=size(se,2); s1=zeros(ns,1); r1=s1; e1=s1; n1=s1;
  for k=1:ns
     s1(k)=max(reshape(se(:,k,:),[],1));
     e1(k)=min(reshape(e0(:,k,:),[],1));
     r1(k)=max(reshape(r2(:,k,:),[],1));
     n1(k)=max(reshape(nk(:,k,:),[],1));
  end

  setdef('r2flag',1); r_=0;
  for i=1:numel(r1)
     if r1(i), r_=r1(i); else r1(i)=r_; end
  end

  if r2flag
       xd={r1,   'r^2'};
  else xd={1./n1,'1/N'}; end

  i=find(xd{1}>1E-16,1,'last');
  xd{3}=xd{1}(i)/8;

  x=xd{1}; y=e1-eg;
  h=loglog(x,y,'o-'); sms(4); hold on
  label(xd{2},'e_0 - e_0^{ref}');

  semilogx(x(end),y(end),'*','Color',getcolor(1));

  xtight(1.05,'x1',xd{3});
  q=logspace(-10,0,6); set(gca,'XTick',q,'YTick',q);

setax(ah(2,2));
ih22=inset({'SE',[0 0.15]},'scale',[1.2 1.2]);

  y=exp(s1);
  h=semilogx(x,y,'o-'); hold on
  label(xd{2},'max(S)');

  xtight(1.05,'x1',xd{3}); ytight(1.2);
  set(gca,'XTick',q);

  i=find(n1>60); m=numel(unique(n1(i)));
  if m>=3
     [p,mu]=polyfit(x(i),y(i),2);

     xx=linspace(0,max(x(i)));

     set(h,'Marker','none');
     h=plot(xx,polyval(p,xx,mu),'r');
     h=plot(x(i),y(i),'o'); sms(h,4);

     title2('S_{\infty}=%.4g',polyval(p,0,mu));
  end

  semilogx(x(end),y(end),'*','Color',getcolor(1));

% -------------------------------------------------------------------- %

if isvar('Jp1') && isvar('Jp2') && Jp1<Jp2
   setax(ah(1,1));
   map_xdata(gca,'lin',[1 L],[Jp1 Jp2],'xlb','J_{\perp}');
elseif isvar('ph1') && isvar('ph2') && numel(ph1)==1 && numel(ph2)==1
   q=[ph1 ph2]/180;
   map_xdata(ah(:,1),'lin',[1 L],q,'xlb','\phi/\pi');
elseif isequal(wsys,'KNL')
   if isvar('JK1') && isvar('JK2')
      map_xdata(ah(:,1),'lin',[1 L],[JK1 JK2],'xlb','J_K');
   end
elseif isequal(wsys,'Spin'), s='';
   q=param.J;
   if size(q,1)>=L-1
      for i=1:size(q,2), x=unique(q(:,i));
         if numel(x)>1
            s=sprintf('J_%d',i);
            break
         end
      end
   end
   if isempty(s), q=param.B;
   if size(q,1)==L
      for i=1:size(q,2), x=unique(q(:,i));
         if numel(x)>1
            if i==1, s='B_z'; else s='D_z'; end
            break
         end
      end
   end
   end

   if ~isempty(s)
      map_xdata(ah(:,1),'lin',[1 L],x([1 end]),'xlb',s);
   end
elseif isequal(wsys,'Kitaev') && isfield(param,'oJ2'), q=param.oJ2;
   if iscell(q) && numel(q)==2
      if isequal(q{1},'--lin')
         map_xdata(ah(:,1),'lin',[1 L],q{2},'xlb','J_2');
      elseif isequal(q{1},'--const')
         map_xdata(ah(:,1),'lin',[1 L],[0 L/param.W],'xlb','x');
      end
   end
end

% -------------------------------------------------------------------- %

