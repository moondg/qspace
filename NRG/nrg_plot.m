
% -------------------------------------------------------------------- %
  if ~usejava('jvm')
   % see i ml => NOJVM_ERR_IN_DEBUG_MODE // Wb,May23,17
	 banner(3,'WRN option -nojvm no longer compatible with debug mode');
     return
  end
  if ~isset('Inrg')
     m=chooseFile([getenv('LMA') '/NRG/*_info.mat']);
     if isempty(m), return; end
     m=regexprep(m,'_info.mat.*','');

     global param
     wblog(' * ','loading %s*',repHome(m));
     Inrg=load([m '_info']); q=load([m '_00']);
     H0=QSpace(q.HK);

     if isfield(Inrg,'om') && isfield(Inrg,'a0') && isfield(Inrg,'rhoNorm')
        Idma=getfields(Inrg,'om','a0','a4','ISpec','rhoT','rhoNorm','ops');
        om=Idma.om;
        a0=Idma.a0;
     end
  elseif ~isvar('H0'), H0=QSpace(Inrg.HK(1));
  end

  if isfield(Inrg,'Itr') && all(Inrg.Itr.Etrunc>0)
     param.D=max(Inrg.NK(:,1));
  end

% --------------------------------------------------------------------- %
ah=smaxis(3,1,'tag',mfilename,'dy',0); addt2fig;

  setdef('param',struct); s={'',''};
  if isfield(param,'istr'), s{1}=[' :: '    param.istr ]; end

  if isfield(param,'sym'), s{2}=param.sym;
  elseif exist('A0','var') && isfield(A0.info,'qtype'), s{2}=A0.info.qtype;
     if exist('IS','var') && ~isempty(s{2}) && ...
        isfield(IS,'SOP') && isfield(IS,'E')
        if ~isequal(IS.E.info.qtype,s{2}), wblog('WRN',...
          'got sym-mismatch (IS.E,A): %s <> %s',IS.E.info.qtype,s{2});
        else s{2}=['{ ' strhcat(IS.SOP.info,'sep',' * ') ' }'];
        end
     end
  end

  if ~isempty(s{2}), s{2}=[' using ' s{2}  ]; end

  nrg_header('-x','sym');
  header(['NRGWilsonQS' cat(2,s{:})]);

  ap=get(ah(3),'Pos'); set(ah(3),'Pos',[ap(1) ap(2)+0.03 ap(3) 0.22]);
  ah=[ ah(1); ah(2); splitax(ah(3),[3,1],'dy',2) ]; setuser(gcf,'ah',ah);
  mvaxis(ah(1),[0 0.03]);

% --------------------------------------------------------------------- %

  yl=[0 4];
  if isset('yl_nrg'), 
     if numel(yl_nrg)==1, if yl_nrg>0, yl(2)=yl_nrg; end
     else yl=yl_nrg; end
  end

  if ~exist('Inrg','var')
     Inrg=add2struct('-','phE0?','EScale?','E0?','EE?','NK?');
  end
  if ~exist('E0','var'), E0=Inrg.E0; end

  l=[];
     if exist('Lambda','var'), l(end+1)=Lambda; end
     if exist('param','var') && isfield(param,'Lambda'), l(end+1)=param.Lambda; end
  if ~isempty(l)
     if numel(l)>1 && norm(diff(l))>1E-12,
        wblog('ERR','inconsistent (param.)Lambda'); disp(l); l=[];
     else l=l(1); end
     if l==2, l=[]; end
  end
  if ~isempty(l) && isset('NRG_LSCALE')
     xsc=[log(l)/log(2), 2*log(3/(l+1))/log(l)];
  else xsc=[]; end

  ylb={'even spectrum (k=0,2,4,...)','odd spectrum (k=1,3,5,...)'}; 

  if isset('EO_FLAG'), ah(1:2)=ah([2 1]); end
  eo_flag=0;

  if isfield(Inrg,'HK') && ~isset('EEflag'), EEflag=0;
     o={'E0',Inrg.E0};
     if isset('qsel_nrg'), o(end+(1:2))={'qsel',qsel_nrg}; end

     if exist('A0','var') && isa(A0,'QSpace')
        q=getDimQS(A0); if max(q(:,1))>1, eo_flag=1;
           ah(1:2)=ah([2 1]);
           ylb(1:2)=ylb([2 1]);
        end
     end

     if exist('QMAP') && ~isempty(QMAP), o={o{:},'Qmap',QMAP}; end
     EE_=getEEdata(Inrg.HK,'ah',ah([2 1]),'EK',yl(2)+2,'-red',o{:});
     EE_=EE_';

     if ~isempty(xsc)
        for i=1:2
           for h=findall(ah(i),'type','line')'
           set(h,'XData',xsc(1)*get(h,'XData')+xsc(2)); end
        end
     end
  else
     if ~isfield(Inrg,'EE'), EE_=EE;
     else
        if ~isequalwithequalnans(Inrg.EE,EE)
            wblog('WRN','got different Inrg.EE / EE !?'); end
        EE_=Inrg.EE;
     end
     EEflag=1;

     if exist('A0','var') && isa(A0,'QSpace')
        q=getDimQS(A0); if max(q(:,1))>1, eo_flag=1;
           ah(1:2)=ah([2 1]);
           ylb(1:2)=ylb([2 1]);
        end
     end

     if size(EE,1)>4000 && exist('Inrg','var') && isfield(Inrg,'EE') && ...
        isequal(size(EE),size(Inrg.EE))
        EE_=reduceEE(EE); % ,'eps',1E-2,'yl',yl(2));
     end
  end

  if ~exist('savemat','var') && ...
    exist('Imain','var') && ~isfield(Imain,'mat') && exist('mat','var')
    Imain.mat=mat;
  end

  k=[1 2]; if eo_flag, k=[2 1]; end

  i1=k(1):2:size(EE_,2); if ~isempty(xsc), x1=i1*xsc(1)+xsc(2); else x1=i1; end
  i2=k(2):2:size(EE_,2); if ~isempty(xsc), x2=i2*xsc(1)+xsc(2); else x2=i2; end
  ii=  1 :  size(EE_,2); if ~isempty(xsc), xx=ii*xsc(1)+xsc(2); else xx=ii; end

setax(ah(1));

  if EEflag
     plot(x2, real(EE_(:,i2)')); hold on; if ~isreal(EE_)
   h=plot(x2, imag(EE_(:,i2)')); blurl(h); mv2back(h); end
  end

  ylabel(ylb{1}); ylim(yl)
  set(gca,'tag',regexprep(ylb{1},'(even|odd).*','$1'));

setax(ah(2));

  if EEflag
     plot(x1, real(EE_(:,i1)')); hold on; if ~isreal(EE_)
   h=plot(x1, imag(EE_(:,i1)')); blurl(h); mv2back(h); end
  end

  ylabel(ylb{2}); ylim(yl)
  set(gca,'tag',regexprep(ylb{2},'(even|odd).*','$1'));

  xtight; ylim(yl); hold on

  if exist('U','var') && isscalar(U) && U==1E-3 && ...
     exist('Gamma','var') && Gamma==U/(12.66*pi) && epsd==-U/2
     eta=[ 1.520483, 3.952550, 9.882118 ]';
     h=plot(xlim, [eta, eta ; 2*eta, 2*eta; 3*eta, 3*eta], 'k--');
  end

  l1=i1(end-1); l2=i2(end-1);
  if E0(l2)<E0(l1)
       k=l2; j=l1; ih=[1 2];
  else k=l1; j=l2; ih=[2 1]; end
  if eo_flag, ih=ih([2 1]); end

  if isfield(Inrg,'HK')
     Ek=cat(2,Inrg.HK(k).data{:});
     E2=cat(2,Inrg.HK(j).data{:}); s=getQ0str(Inrg.HK(k));
  else
     Ek=EE_(:,k); E2=EE_(:,j); s='';
  end

  show_deg(Ek,ah(ih(1)));

  o={ 0.98, 0.95,'VerticalAl','top','HorizontalAlign','right', ...
     'Units','normalized','EdgeColor','k','Margin',5,...
     'BackGroundColor','w','FontSize',10};
  text(o{1:2}, sprintf('ground state %s\\newlineE_g=%.8g',s,Inrg.phE0),o{3:end});

  show_deg(E2,ah(ih(2)));
  if isfield(Inrg,'HK')
     h=text(o{1:2}, sprintf('lowest %s',getQ0str(Inrg.HK(k-1))),o{3:end},'FontSize',8);
     set(h,'EdgeColor',[.99 .99 .99]);
  end

setax(ah(3));

  plot(x2,real(E0(i2)),'bo-','Disp','\DeltaE even'); hold on
  plot(x1,real(E0(i1)),'ro-','Disp','\DeltaE odd'); sms(2); clear q

     q=ylim; l=1.2*max(abs(E0(5:end))); if max(abs(q))>l
        q=q* (l/max(abs(q))); ylim(q);
     end

  if ~isreal(E0)
     h=plot(x2,imag(E0(i2)),'b','Disp','imag. part'); blurl(h); mv2back(h);
     h=plot(x1,imag(E0(i1)),'r'); blurl(h); mv2back(h);

     if isfield(param,'gamma0'), q=ylim;
        dd=imag(E0)-param.gamma0*length(Inrg.E0)./Inrg.EScale;
        h=plot(x2,dd(i2),'b','Disp','imag. part');
        h=plot(x1,dd(i1),'r'); ylim(q);
     end
  end

  xtight; ytight(1.1); % set(gca,'XTickLabel',[]);
  if exist('q','var'), ylim(q); end

  ylabel('\DeltaE_0')
  legdisp('Location','NorthWest','orientation','horiz','-detach',...
    'dx',[0.03 0.01],'fs',8); % ,'xsc',0.8

  set(gca,'tag','deltaE0');

setax(ah(4));

  if isfield(Inrg,'Itr') && all(Inrg.Itr.Etrunc>0)
     h=plot(xx,Inrg.EK(:,1:2),'o-'); set(h(2),'Color','r'); blurl(h,0.5,'LineW',1);
     set(h(1),'Disp','max(E_K)');
     set(h(2),'Disp','min(E_D)'); mv2back(h(2));

     if isfield(Inrg,'EK')
        et=Inrg.EK(10:end-1,:); if isempty(et), et=Inrg.EK(end,:); end
        et=reshape(et(:,1:2),[],1); et(find(isnan(et)))=[];
        if numel(et)>1
           ylim(mean(et)+[-1 1]*5*std(et));
        end
        et=Inrg.Itr.Etrunc;
     else
        et=Inrg.Itr.Etrunc;
        y=ylim; if y(2)>2*et, ylim([0 2.5*et]); end
     end
     h=ymark(et,'Color',[.8 .8 1],'LineW',2);

     h=ylabel('E_{trunc}','Rot',90); % ,'Color','b'
     mvlabel('y',-.02);
  else
     h=plot(xx,max(EE_),'bo-','LineW',2,'Disp','E_{max}');
     blurl(h,0.5,'LineW',1);
     h=ylabel('E_K ; E_D','Rot',90); % ,'Color','b'
  end
  clear EE_

  sms(4); hold on

  h=legdisp('Location','SouthWest','Orientation','horiz','-detach',...
   'dx',[ 0.03 -0.01],'dy',-0.15,'fs',8,'xsc',0.8);

  set(gca,'tag','Etrunc');

setax(ah(5));

  q=max(Inrg.NK,[],1);
  n=size(Inrg.NK,2); if n==4, j=[1 3]; else j=1:n; end

  h1=plot(xx,Inrg.NK(:,j(1))); y=ytight(1.1); hold on
  h2=plot(xx,Inrg.NK(:,j(2:end))); h=[h1;h2];
  set(h(1),'Disp','N_K');  if n>1
  set(h(2),'Disp','N_{tot}'); end

  if n==4
     set(h(2),'Color',[0 .5 0]); set(h(2),'LineSt','--');
     s=num2cell(q); for i=1:numel(s), s{i}=int2str2(s{i}); end
     postext({'NE',[0.035 -0.2]},...
     'N_K^{max}=%s (%s);  N_{tot}^{max}=%s (%s)', ...
       s{[1 2]}, s{[3 4]},{'FontSize',8});
  elseif n>1
     postext({'NE',[0.02 -0.05]},'N_K^{max}=%d',q(1));
  else
     set(h(1),'Disp',sprintf('N_K^{max}=%d',q(1)));
  end

  if n>1
     i=[find(Inrg.NK(:,j(1))==Inrg.NK(:,j(2))); Inrg.N];
     i=i(diff(i)>1);
     if ~isempty(i)
        h2=plot(xx(i),Inrg.NK(i),'o','Color',[1 .5 .5]);
        mv2back(h2); sms(h2,4);
        xmark(xx(i+1),'LineW',4,'Color',[1 .7 .7]);
     end
  end

  if 0 && max(q)/min(q)>32
     set(gca,'YLimMode','auto','YScale','log')
     ytight(1.1,'y1',10);
     set(gca,'YTick',10.^(1:2:9)); %,'FontSize',8 grid on
  else ytight; end

  if isfield(param,'Nkeep')
     h=ymark(param.Nkeep,'LineW',2,'Color',[.8 .9 .8]);
     mv2back(h);
  end

  ylim(y);

  xlabel('Wilson shell (NRG iteration) k') %,'max(E_k)');
  if isAbelian(H0),l='N_{kept}'; else l=['N_{kept}' 10 '(multiplets)']; end
  ylabel(l); % 'Color',[0 .5 0] % mvlabel('y',-0.02);

  legdisp('Location','SouthWest','orientation','horiz','-detach', ...
    'dx',[.04 -.005],'dy',-0.2,'xsc',0.8,'fs',8);

  set(gca,'tag','Nkept');

  xl=[1 max([x1,x2])];
  set(ah,'XLim',xl,'FontSize',12);

  set(ah(3:end-1),'XTickLabel',[]);
  set(ah(3:end),'TickLength',[0.005 0.01],'FontSize',12)

  set(ah(1:2),'TickLength',[0.005 0.01]); % ,'XColor',[1 1 1]*0.6)

if ~isempty(xsc)
   p=[ get(ah(end),'Pos'); get(ah(1),'Pos') ]; d=0.053;
   p=[ p(1,1), p(1,2)-d, p(1,3), p(2,2)+p(2,4)-p(1,2)+d ];

ah(end+1)=axes('Position',p);
   set(gca,'XLim',(xl-xsc(2))/xsc(1),'Color','none','YColor',[1 1 1]*0.999,'YTick',[]);
   set(gca,'FontSize',10,'TickLength',[0.002 0.002]); grid on

for i=1:5, setax(ah(i)); if i<4, continue; end
   h=leg2front; leg_mvtext(h,-0.15);
   if i==4
      xlabel(sprintf(['NRG iteration @ \\Lambda=2\n' ...
      'actual NRG iteration (\\Lambda=%g)'],Lambda));
   end
end
end

% --------------------------------------------------------------------- %

  if exist('wsys','var')
     switch wsys
        case {'KondoJH'}
        header(istr)
        header('fleft',sprintf(['NRG N=%g, \\Lambda=%g, D=%g; ', ...
         'AM: U=%g, \\epsilon_d=%g, \\Gamma_{(1)}=%g, \\Gamma_2=%g, J_H=%g' ],...
          N,Lambda,Nkeep,U,epsd,Gamma,Gamma2,JH));
     end
  end

  drawnow

% --------------------------------------------------------------------- %

