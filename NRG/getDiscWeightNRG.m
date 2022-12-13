function [dw,Iout,IR]=getDiscWeightNRG(varargin)
% function [dw,I,IR]=getDiscWeightNRG(NRGtag [,opts])
%
%    get estimate for discarded weight (dw) from reduced density
%    matrices, considering backward update from overall ground
%    state space (default NRGtag: './NRG/NRG').
%
% Options
%
%   'n0',...   how many times to back-propagate rho until min(rho) is obtained
%   'k0',...   consider first NRG_<it> iterations with it<k0 exact
%   'beta',..  local 1/T to initialized rho to ground state space (100)
%   'chi',..   lower range of eig(rho) to use (default: 0.05 = last 5%)
%
%   '-E'       use energy eigenbasis (E-basis s, rather than
%              eigenbasis of rho, R-basis r)
%   'Etr',..   truncation energy (as used in histogram for energy)
%              (default: Inrg.Itr.Etrunc(1)) 
%
%   '-q'       quiet (do not show iteration k while calculating rho)
%   '-fig'     show summarizing figure
%
% The returned (estimated) discarded weight is given by dw = mean(I.e1).
% where the returned info structure I contains (amongst others)
% the following fields:
%
%   .chi   lowest range of eigenvalues used (0.05 = 5%)
%   .e1    discarded weight based on the lowest weight states
%          in the reduced density matrices (RDM) R(n-n0:n);
%          (lowest chi percent in terms of number of states)
%   .e2    discarded weight based on the top energy states
%          in the reduced density matrices (RDM) R(n-n0:n);
%          however, since rho_r and E_r are strongly related,
%          e2 is typically only slightly larger than e1 (R-basis).
%   .dw    estimated discarded weight for each Wilson shell
%   .se    Shannon entropy
%   .beta  1/T (in rescaled units) used when initializing the density
%          matrix in the kept energy states (projection into
%          the ground state space; default: beta=100).
%
% Wb,Jan24,11

% adapted from: getrhoNRG_resall.m ; Wb,Aug17,12

  getopt('init',varargin);

     n0    =getopt('n0',-1);
     k0    =getopt('k0', 0);
     beta  =getopt('beta',200);

     smooth=getopt('smooth',0);

     if smooth<=0
        chi =getopt('chi',0.05);
     else
        chi =getopt('chi',0.01);
     end

     vflag =getopt('-v');
     qflag =getopt('-q');
     Eflag =getopt('-E');
     Etr   =getopt('Etr',[]);
     showfig=getopt('-fig');

     yl    =getopt('yl',[]);

     aflag =getopt('-a');
     xl31  =getopt('xl31',[]);
     yl31  =getopt('yl31',[]);

     nochi =getopt('~chi');
     nx    =getopt('nx',[]);

     fflag =getopt('-ff');

  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg
     if narg>1 || ~ischar(varargin{1}), varargin
     error('Wb:ERR','invalid NRGdata tag'); end
     nrg=varargin{1};
  else nrg='./NRG/NRG'; end

  if numel(chi)~=1 || (chi<=0 || chi>=1)
     error('Wb:ERR','\n   ERR invalid chi (%g)',chi);
  end

  if smooth>0 && smooth<2, error('Wb:ERR',...
    '\n   ERR invalid smoothing parameter for disc. weight (%g)',smooth);
  end

  if nrg(1)~='/' && isempty(findstr(pwd,'Data')), cto lma, end
  ff=getnrgfiles(nrg,'-f');
  N=numel(ff); 

  if ~N, error('Wb:ERR',...
  '\n   ERR no NRG files ''%s'' found',nrg); end

  Inrg=load([nrg '_info']);
  Lambda=Inrg.Lambda; Nkeep=Inrg.Nkeep;

  if isempty(Etr)
     Etr=Inrg.Itr.Etrunc(1);
  end
  if numel(Etr)~=1 || Etr<=0, Etr
     error('Wb:ERR','\n   ERR invalid usage (Etr)');
  end

  mat=ff(N).name; load(mat);
     if ~isempty(QSpace(HK)), error('Wb:ERR',...
     'HK must be empty at last iteration (%s) !??',mat); end

  AK=AD; AD=QSpace;
  HK=HD; HD=QSpace;

  if fflag, getbase FC
      if isempty(FC), wblog('ERR',...
        'failed to find FC operators (got -ff flag)'); fflag=0;
      else
         for i=1:numel(FC)
            x(i)=normQS(FC(i))^2;
         end
      end
   end
   e3=[];

 % -------------------------------------------------------------------- %
  if aflag
     if isequal(n0,-1), n0=[]; end
     n1=1; n2=N;
  else
     if isequal(n0,-1),
        q=max(2, 4/log(Lambda));

        d=getDimQS(Inrg.ops.A0); d=d(end);
        k=max(Inrg.NK); if numel(k==4) k=k(2); else k=k(1); end
        q=max(q, log(k)/log(d));

        n0=floor(q)+[1 2];
     end

     if ~isempty(n0), n1=min(n0); n2=max(n0);
     else, n1=4; n2=4;
     end
  end

  for k=N:-1:2

     [R,I0]=getrhoQS(HK,beta);
     IR(k)=struct('k1',[],'k2',k,'R',R,'R0',R,'I0',I0);

     if fflag && ~isempty(HD) && ~isempty(HD.Q), for i=1:numel(FC)
        R=getrhoQS(HK,beta,'-sqrt');
        q=contractQS(contractQS(AK,2,R,1),2,FC(i),1);
        q=contractQS(q,'13',AD,'13*'); %  R(o),R'
        l=numel(q.Q);

        a=HD; for j=1:numel(HD.data), a.data{j}=1./a.data{j}; end
        q=contractQS(q,l,diag(QSpace(a)),1);
        q=normQS(q)^2;

        e3(k,i)=( x(i)*q*Inrg.ops.ff(k)^2 )^2;
     end, end

     if k>1
     load(ff(k-1).name,'HK','HD'); end

     kk=unique([ (k-1)+(1:n2), N ]); kk(find(kk>N))=[];

     if ~qflag, fprintf(1,...
       '\r   %s %4g/%g (%g)... \r',mfilename,k,N,numel(kk)); end

     for k2=kk
        R=IR(k2).R;

        q=contractQS(AK,'2*',IR(k2).R,2);
        IR(k2).R=QSpace(contractQS(AK,'23',q,'32'));
        IR(k2).k1=k-1;

        if k2-k+1<n1
           continue
        end

        s.I0=IR(k2).I0;
        s.k1=IR(k2).k1;
        s.k2=IR(k2).k2; if k2==N, RR(k)=IR(k2); end

        [rr,Ie]=eigQS(IR(k2).R); q=Ie.EK; n=numel(q.data);
        for j=1:n, d=q.data{j}; d(find(d<1E-16))=nan;
            q.data{j}=fliplr(-log(d));
        end
        s.ES=q;

        [s.se,s.rr,s.Is]=SEntropy(QSpace(IR(k2).R),'-xm');

        l=ceil((1-chi)*length(s.rr));
        if smooth>1
           s.e1=(1./((s.rr/s.rr(l)).^smooth+1))*s.rr';
         % this is the Fermi function on a logscale using beta:=smooth(=4)
         % [with the variable beta already used in a different context above]
         % (NB! smooth>1 required to suppress contribution of rr>>rr(l) data)
         % In contrast to mean() below, here e1 is the SUM of all discarded
         % weights together with a smooth transition weighting function.
         % Wb,May12,14 // tags: SMOOTH
        else
           ix=find(s.rr<=1.001*s.rr(l));
           s.e1=mean(s.rr(ix));
        end

        if s.e1<0
           if abs(s.e1)>1E-12, wblog('ERR','got s.e1=%g',s.e1); end
        end

        if isempty(HD.data)
           s.erd=[]; s.e1=[]; s.e2=[];
        else
           R=struct(IR(k2).R);
           [i1,i2,I]=matchIndex(R.Q{1},HK.Q{1});
           dz=getzdim(QSpace(R),2,'-p');

           for i=1:numel(i1)
              j1=i1(i); ri=R.data{j1};
              j2=i2(i); ei=HK.data{j2};

              di=repmat(dz(j1),size(ri,1),1);

              if Eflag
                 ri=diag(ri); ei=ei';
              else
                 [u,ri]=eig(ri+ri'); ri=0.5*diag(ri);
                 ei=real(diag(u'*diag(ei)*u));
              end
              R.data{j1}=[ei,ri,di];
           end
           s.erd=cat(1,R.data{:});

           ex=(1-chi)*max(s.erd(:,1));

           if smooth>1
              s.e2=(1./((s.erd(:,1)/ex).^smooth+1))'*(s.erd(:,2).*s.erd(:,3));
           else
              ix=find(s.erd(:,1)>ex);

              s.e2=(s.erd(ix,2)'*s.erd(ix,3))/sum(s.erd(ix,3));;
           end

        end

        r=struct(diag(IR(k2).R));
        [s.rd,i]=sort(cat(1,r.data{:}));
        s.ed=min(s.rd);

        if k2-k<n2, II(k-1,k2-k+1)=s; end
        if k2==N, II(k-1,n2+1)=s; end

        if abs(trace(IR(k2).R)-1)>1E-8, IR(k2).R
        error('Wb:ERR','invalid density matrix R'); end
     end

     if k>1
        if fflag, v={'AD'}; else v={}; end
        load(ff(k-1).name,'AK',v{:},'HD');
     end

  end
  if ~qflag, fprintf(1,'\r%50s\r',''); end

  Iout=add2struct('-',nrg,Nkeep,Etr,beta,smooth,e3,n0,chi,fflag,aflag);
  if aflag && nargout<3, Iout.IA=II; end

  II=II(:,[n0 end]);

  Iout.param=Inrg.param;

  se=getdatafield(II,'se','-q'); se=seteps2nan(se); Iout.se=se;
  e1=getdatafield(II,'e1','-q'); [e1,gotneg_e1]=seteps2nan(e1);
  e2=getdatafield(II,'e2','-q'); e2=seteps2nan(e2);

% -------------------------------------------------------------------- %
% NB! e2 is preferential in particular also when using Eflag
% due to the large spread of lowest-weight states!
% NB! Eflag leads to markedly smooter data when plotting
% discarded weight vs. Wilson shell!
% Wb,Nov19,15
% -------------------------------------------------------------------- %

  e=e1(k0+1:end,:); e(find(isnan(e)))=0; Iout.e1=max(e, [],1);
  e=e2(k0+1:end,:); e(find(isnan(e)))=0; Iout.e2=max(e, [],1);

  dw=mean(Iout.e2); Iout.dw=max(e2,[],2);

  if gotneg_e1~=0 && (dw<=0 || ~qflag || vflag)
     wblog('WRN','got e1<0 !! (%.3g%%)',100*gotneg_e1);
  end

  if vflag
  wblog('<i>','max. discarded weight: %.3g',dw); end

% -------------------------------------------------------------------- %
  if ~showfig, return; end
  EE=Inrg.EE; nE=min(512,size(EE,1));

  k=round(size(II,1)/2);
  k1=cat(2,II(k,:).k1); if norm(diff(k1))>1E-12
     error('Wb:ERR','inconsistent k1''s !??'); end
  k2=cat(2,II(k,:).k2);
  n0_=k2(1)-k1(1);

  t=mfilename;

  if Eflag, bs='_s'; else bs='_r'; end

ah=smaxis(8,1,'tag',t,'dx',0.10,'dy',0.007,'DY',0.02); 

  s=repmat({''},1,4); getbase param

  if isfield(param,'istr'),s{1}=param.istr; end
  if isfield(param,'sym'),s{2}=param.sym;
  elseif isfield(param,'wsys'),s{2}=param.wsys; end

  if k0>0, s{3}=sprintf('k_0=%g, ',k0+1); end
  if Eflag
       s{4}='{\color{red}using energy eigenbasis}';
  else s{4}='using eigenbasis of \rho';
  end

  evalin('caller','nrg_header'); header(['%M :: %s (%s)' 10 ...
    '\\chi=%.3g%% (%sn_0=%g, m=%g, \\beta=%g, smooth=%g, %s)'],...
     s{1:2},100*chi,s{3},max(n0),numel(n0),beta,smooth,s{4});
  set(gcf,'DefaultAxesLineW',1.5); addt2fig Wb

  ah=mat2cell(ah,[2,3,3],1); dy=0.012;

  mvaxis(ah{1}(1),[0  dy+.003]);
  mvaxis(ah{1}(2),[0  dy+.006]);
  mvaxis(ah{2}(1),[0  dy+.006]);
  mvaxis(ah{2}(2),[0  dy+.009]);
  mvaxis(ah{2}(3),[0  dy+.012]);

  mvaxis(ah{3},[0  -.006]);

  p=get(ah{3,1}(3),'Pos'); p(3)=0.6*p(3); p(4)=2.2*p(4);
  set(ah{3,1}(1),'Position',p);

  set(ah{3,1}(2),'Position',[p(1), p(2)+p(4)+0.005, p(3), p(4)/3]);

  set(ah{3,1}(3),'Position',[p(1)+p(3)+0.005, p(2), p(3)/6, p(4)]);

setax(ah{1,1}(1))

  cg=[.8 .6 .4];

  lo={{'Color','k'},{'Color',cg}};

  if ~isempty(EE)
     n=size(EE,2);
     i=2:2:n; h=plot(i,EE(1:nE,i),lo{2}{:}); set(h(1),'Disp','odd'); hold on
     i=1:2:n; h=plot(i,EE(1:nE,i),lo{1}{:}); set(h(1),'Disp','even');
     ylim([0 2*Inrg.ops.ff(end)]);
  end

  rr=cat2(1,II(:,end).rr); rr(find(rr<=0))=nan;
  rr=-log(rr); nr=min(512,size(rr,2));

  legdisp({'NE',[0 0]},'Orientation','horizontal','-eraset','pick',[2 1]);
  ylabel(['energy' char(10) 'flow diagram']);
  set(gca,'XTickLabel',[]);

setax(ah{1,1}(2))

  lo={{'Color','k'},{'Color',cg}};

  n=size(rr,1);
  i=1:2:n; h=plot(i,rr(i,1:nr),lo{1}{:}); set(h(1),'Disp','even'); hold on
  i=2:2:n; h=plot(i,rr(i,1:nr),lo{2}{:}); set(h(1),'Disp','odd');

  sms(2); xtight; ylim([-0.05 5]);
  legdisp({'NE',[0 0]},'Orientation','horizontal','-eraset');

  ylabel(['entanglement' char(10) 'flow diagram']);
  set(gca,'XTickLabel',[]);

setax(ah{2,1}(end))

  h1=semilogy(e1); hold on
  set(h1(1),'Disp','e1: <eig(\rho_n)> using \chi_N');

  h2=semilogy(e2);
  set(h2(1),'Disp','e2: <eig(\rho_n)> using \chi_E');

  blurl(h1,'cfac',0.85);

  if ~isempty(yl), ylim(yl);
  else ytight(1.2); y=ylim; if y(2)<1E-4, y(2)=1E-4; ylim(y); end
  end
  set(ah{2,1}(end),'XLim',get(ah{1,1}(1),'XLim'),'YTick',logspace(-15,0,6));

  label('Wilson iteration n','\epsilon_n^D');

  legdisp({'SE',[ 0.01  0.01]},... %'Orientation','horizontal',...
    'dy',[-0.2],'ysc',0.6,'fs',10,'-detach');
  h=ymark(dw,'r--','-fg','istr',{'\epsilon^D',0.10,'VerticalAl','bottom'});

setax(ah{3,1}(1)); j=1; n=size(II,1); h=zeros(n,1);

  ex=zeros(n,1);
  for i=1:n, setax(ah{3,1}(1)); j=1;
     erd=II(i).erd; if isempty(erd) || i<k0+1, continue; end
     j=find(erd(:,2)>0);

     ex(i)=mean(erd(find(erd(:,1)>0.8*max(erd(:,1))),1));

     h(i)=semilogy(erd(j,1),erd(j,2),'.','Color',getcolor(i));
     hold on
  end

setax(ah{3,1}(1));

  k=find(h); h=h(k); ex=ex(k);

  if isempty(nx)
     ix=[ find(ex==min(ex)), find(ex==max(ex)) ];
     sx={'n_{max}', 'n_{min}' };
  else
     ix=zeros(size(nx));
     for i=1:numel(nx), ix(i)=find(k==nx(i)); end
     sx={ 'n_1', 'n_2' };
  end

  om={'Color',[1 .3 .3],'LineW',2}; o={om{1:2},'Marker','o'}; % 'ro'
  i=ix(1);
    set(h(i),o{:},'Disp',sprintf('%s=%g',sx{1},k(i)));
    sms(h(i),3);

  om={'Color','k','LineW',2}; o={om{1:2},'Marker','d'}; % 'kd'
  i=ix(2);
    set(h(i),o{:},'Disp',sprintf('%s=%g (N=%g)',sx{2},k(i),Inrg.N));
    sms(h(i),3);

  mv2front(h(ix));
  hx=h; hx(ix)=[]; blurl(hx,'cfac',0.7);

  label(['rescaled energy E^{[n]}' bs], ['\rho^{[n;n_0]}' bs]);
  set(gca,'YTick',logspace(-20,0,5));

setax(ah{2,1}(1));

  n=size(II,1); lo={'Color',[1 1 1]*.65};
  i=2:2:n; h1=plot(i,(se(i,:)),lo{:});   set(h1(1),'Disp','even'); hold on
  i=1:2:n; h2=plot(i,(se(i,:)),'--',lo{:}); set(h2(1),'Disp','odd');

  h=[h1(end),h2(end)];
  ds='n_0\rightarrow \infty';
  set(h,'Color','b','LineW',2); set(h(1),'Disp',ds); mv2front(h);

  legdisp({'NE',[0 0]},'Orientation','horizontal',...
    'dy',[0 0 -0.2]); %,'pick',[1 3 2] 
  y=ytight(1.1,'y1',1); ylabel('entropy S'); % mvlabel('y', 0.04)
  if y(2)<3,
       set(gca,'YTick',0:.5:y(2));
  else set(gca,'YTick',0:y(2)); end

setax(ah{2,1}(2));

  plot(Inrg.NK(:,1),'o-'); sms(4);
  ylabel('N_{keep}'); mvlabel('y',0.02);

  ytight(1.1,'y1',0);

  y=get(gca,'YTick'); if numel(y)<=2,
  set(gca,'YTick',linspace(0,y(end),3)); end

  h=ah{1,1}(1);
  set([ah{1,1}; ah{2,1}],'XLim',get(h,'XLim'),'TickLength',[0.005 0.005],...
    'XTick',get(h,'XTick'), 'XTickLabel',get(h,'XTickLabel'),'FontSize',10);

  setax(ah{2,1}(end));
  set(gca,'XTickLabel',get(gca,'XTick'));

setax(ah{3,1}(1));

  ERD=cat(1,II.erd);
  n=ceil(22/sqrt(Lambda));
  xb=linspace(0,Etr,16);

  wbhist([],[],xb); n=size(II,1); Np=0;
  for i=1:n
     erd=II(i).erd; if isempty(erd), continue; end
     wbhist(erd(:,1),erd(:,2)); Np=Np+1;
  end

  [x,y]=wbhist; y=(y/Np)./diff2(x,'len');

  s='\rho(E)'; dwc=nan;

  if 1
    i=2:length(x)-1;

    i=i(find(y(i))); % plot(x,y,'k');
    p=polyfit(x(i),log(abs(y(i))),1); xx=xlim; xx=linspace(xx(1),xx(2));

    h1=plot(x(i),abs(y(i)),'r','Disp','cumulative weight');
    h2=plot(xx,exp(polyval(p,xx)),'k');
    mv2back(h2); blurl([h2]);

    q=p; q(2)=log(abs(q(1)));
    h=plot(xx, exp(polyval(q,xx)),'-.','Color',[.3 .3 .3],...
      'Disp',sprintf('%s \\sim e^{%.4g{\\cdot}E}',s,p(1)));

    dwc=exp(p(1)*Etr);

  else
    i=find(ERD(:,2)>1E-12);
    i=find(ERD(:,2)>exp(mean(log(ERD(i,2)))));

    p=polyfit(ERD(i,1),log(ERD(i,2)),1);

    xx=xlim; xx=linspace(xx(1),xx(2));
    h=plot(xx,exp(polyval(p,xx)),'k','Disp',...
      sprintf('%s \\approx %.3g e^{%+.4g{\\cdot}E}',s,exp(p(2)),p(1)));
  end
  header('fleft','-append',' \rho(E) fit: %.4g',p(1))

  Iout.p=p;

  if ~isempty(xl31)
     if numel(xl31)==2, xlim(xl31); else axis(xl31); end
  else xtight(1.05); %'x1',-0.25);
  end

  if ~isempty(yl31)
     if numel(yl31)==2, ylim(yl31); else axis(yl31); end
  else ytight(1.1,'y1>',1E-15); end

  o={'-detach','fs',10};
  if nochi
     legdisp('pick',[3 4],o{:},'dx',[0.02 0.005],'dy',[0 0.05]); y=0.65;
  else
     ymark(dw,'r--','-fix','Color',[1 .4 .2],'-fg');
     xmark(Etr,'r--');

     legdisp('pick',[3 4],'dy',[0 0.10],o{:},'ysc',0.7,'-detach');
  end

  if ~isfield(Inrg,'Itr') || isempty(Inrg.Itr) || all(Inrg.Itr.Etrunc<=0)
     s=sprintf('\\Lambda=%g, N_K=%g\nn_0=%g',Inrg.Lambda,Inrg.Nkeep,n0_);
  else
     q=Inrg.Itr.Etrunc(1);
     e=norm(diff(Inrg.ops.ff(end-3:end,:)))/abs(Inrg.ops.ff(end));
     if e>1E-2, wblog('WRN','couplings still vary (e=%.3g)',e); end
     s=sprintf('\\Lambda=%g,E_K=%.2g\nn_0=%g',Inrg.Lambda,q,n0_);
  end
  postext({'E',[0.02 0.12]},s,{'FontSize',10});

  legdisp({'SW',[-0.02 -0.005]},'pick',[2 1],'fs',10,'-eraset',...
    'dy',-0.05,'xsc',{0.8 [0 1.5]},'-detach','ysc',0.8);

setax(ah{3,1}(2));

  xb=linspace(0,ceil(max(ERD(:,1))),32);
  wbhist([],[],xb); n=size(II,1);

  for i=1:n
     erd=II(i).erd; if isempty(erd), continue; end
     wbhist(erd(:,1),[]);
  end

  [ee,nn,xb]=wbhist; nn=nn/Np;

  plot(ee,nn,'o-'); sms(3); hold on
  ylabel('\nu(E)');

  ytight(1.2,'y1',0);
  set(gca,'XLim',get(ah{3,1}(1),'XLim'),...
    'XTick',get(ah{3,1}(2),'XTick'),'XTickLabel',[]);

setax(ah{3,1}(3));

  wbhist([],[],linspace(-60,0,48)); n=size(II,1);

  for i=1:n
     erd=II(i).erd; if isempty(erd), continue; end
     j=find(erd(:,2)>0);
     wbhist(log(erd(j,2)),[]);
  end
  [logr,nn,xb]=wbhist; nn=nn/Np;

  semilogy(nn,exp(logr),'o-'); hold on
  sms(3); % scalecolor('-flip','alpha',2);

  if ~nochi
     h=ymark(dw,'r--','-fg'); x=xlim;

     if smooth>0
          s=sprintf('\n    [using smooth=%g]',smooth);
     else s=''; end

     text(x(2)+0.3*diff(x),2*dw, ['\rightarrow  ' sprintf(... 
     '\\epsilon^D_{\\chi=%g} = %s%s',chi,num2tex(dw,'%.2E'),s)], ...
     'Color','r','VerticalAl','top');

     h=ymark(dwc,'k--');
     text(x(2)+0.3*diff(x),2*dwc, ['\rightarrow  ' sprintf(... 
     '\\epsilon^D_{\\chi=%g} = %s%s',chi,num2tex(dwc,'%.2E'),s)], ...
     'VerticalAl','top');
  end

  xtight(1.2,'x1',0); 
  xlabel('\nu(\rho)');

  set(gca,'YLim',get(ah{3,1}(1),'YLim'),'YTick',...
  get(ah{3,1}(1),'YTick'),'YTickLabel',[]); %,'YAxisLocation','right');

  setax(ah{2,1}(end,1));

end

% -------------------------------------------------------------------- %

function [dd,neg]=seteps2nan(dd,eps)

   if nargin<2, eps=0; end

   if nargout>1
      i=find(dd<0); neg=numel(i); if neg
        neg=neg/numel(find(~isnan(dd)));
        dd(i)=abs(dd(i));
     end
   end

   i=find(dd<=eps); e=norm(dd(i));
   if e>1E-12
      s=sprintf('WRN skipping small weight e=%g !??',e);
      if e<1E-6, wblog('WRN',s'); else
      error('Wb:ERR','\n   ERR %s',s); end
   end

   dd(i)=nan;

end

% -------------------------------------------------------------------- %

