function [Iout,EK,EH]=plotRdata(HAM,varargin)
% function Iout=plotRdata(HAM [,opts])
% Options
% 
%   'kQ',k  at which k to indicate q-labels
%           default: ymax=4, otherwise 'kQ',{k,ymax}  
% 
% Examples
% plotRdata(HAM,'k2x',{[-2 6],'J_{2}'},'dk',3)
% [Iout,EK,EH]=plotRdata(HAM,'k2x',{[min(ph),max(ph)]/pi,'\phi/pi'},'yl12',[0 6],'ylH',[-0.02 0.5])
%
% SEE ALSO MEX/plotRdata(RR) using RR as returned to tst_Hamilton1D!
% 
% Wb,Apr18,14

  getopt('init',varargin);
     k2x    = getopt('k2x',{}); % { x(k), 'xlabel' }
     semin  = getopt('semin',[]);
     yl12   = getopt('yl12',[0 8]);
     ylH    = getopt('ylH', [-0.05 1]);
     kQ     = getopt('kQ',[]);
     isw    = getopt('isw',[]);
     eoflag = getopt('-eo');
     ksplit = getopt('dk',1);
     xstr   = getopt('xstr','');
     dyl    = getopt('dyl',[]);
     vflag  = getopt('-v');
  getopt('check_error');

  L=numel(HAM.mpo); h=zeros(L,1);
  D=-1; % D=max(getPsiDim(HAM,'-a','-q'),[],1);

ah=smaxis(2,2,'tag',mfilename,'dy',0.10);
addt2fig Wb

  odir=1; kc=L; EH=QSpace; EK=QSpace; R=QSpace; Rl=QSpace; Ql=[]; qloc=[];
  r2=nan(1,L); dsw={};
  ek=nan(L,1);
  dd=nan(L,2);
  DX=nan(L,2);

% NB! k<kx vs. k>kx only flips q-labels with dual labels
% but leaves the eigenspectrum unchanged
% ==> for bond entanglement spectra for global singlets
%     q -> q_bar leads to same entanglement spectrum *by value*
%     while q and q_bar may very well have different 
%     spectra themselves when considering k<kx or k>kx individually!
% Wb,Mar22,15
% kx=L/2; => can lead to discontinuities in entanglement flow diagram!
  kx=0;

setax(ah(2,1)); h=gobjects(0); largf=0; b=0;

  I1=load_dmrg_data(HAM,2,'info');
  I2=load_dmrg_data(HAM,L-1,'info');

  if isempty(isw), isw=0; end
  if isw<0
     i1=[I1.isw, I2.isw];
     i2=floor([size(I1.dsw,1),size(I2.dsw,1)]/2);

     if ~isempty(i1) && ~isempty(i2), j1=min(i1); j2=min(i2);
        if norm(j1-j2)>1
           wbdie('unexpected sweep setting !?'); end
        isw=max(0,min(j1,j2)+(1+isw));
     else isw=0;
     end
  end

  for k=L:-1:1
     if 1 || largf, fprintf(1,...
       '\r   loading %3g/%g (%.3gM) ...  \r',k,L,b/2^20);
     end

     I =load_dmrg_data(HAM,k,'info');
     HK=load_dmrg_data(HAM,k,'HK');

     q=whos('I','HK'); b=sum([q.bytes]); if b*L>1E6, largf=1; end

     if isempty(HK.Q) || isempty(find(HK.info.itags{1}=='E'))
        d=getDimQS(HK); if ~isempty(d), dd(k,:)=d(:,1); end
     else
        AK=load_dmrg_data(HAM,k,'AK');
        d=getDimQS(AK); if ~isempty(d), dd(k,:)=d(:,2); end
     end

     if ~isempty(HK)
        [ee,Ih]=eigQS(HK); ee=ee(:,1);
        k0=get_kout(HK.info.itags);
        EH(k0)=QSpace(Ih.EK)-min(ee); EH(k0).info.itags=HK.info.itags;
     end

     k1=get_kout(I.itags); if isempty(k1), k1=k; end
     D=max([D; I.Dk(:)]);

     if isfield(I,'r2'), r2(k1)=I.r2; end
     if isfield(I,'dsw'), dsw{k1}=I.dsw; end

     if isfield(I,'ek') && ~isempty(I.ek), q=I.ek;
        if iscell(q), q=q{end}; end
        if ~isempty(q)
           ek(k1,1:size(q,1))=q(:,end);
        end
     end

     if isfield(I,'rho') && ~isempty(I.rho)
        Rl(k)=I.rho; ql=I.rho.Q;
        Ql{k}=uniquerows(cat(1,ql{:}));
     end

     if isfield(I,'Rho')
        if isw
           if k>1 && k<L, j=2*isw; else j=isw; end
        else j=numel(I.Rho);
        end

        Rho=QSpace(I.Rho(j)); if isdiag(Rho)>1, Rho=diag(Rho); end
        k2=get_kout(Rho.info.itags);
        if isfield(I,'DX'), DX(k2,:)=I.DX(:,1); end
        R(k2)=QSpace(Rho);
     else
        k2=get_kin(I.Psi.info.itags);
        if isfield(I,'DX'), DX(k2,:)=I.DX(:,1); end

        if isfield(I,'Rho'), R(k2)=I.Rho;
        elseif isfield(I,'Psi')
           if k<kx
                R(k2)=QSpace(contractQS(I.Psi,I.Psi,'!2*'));
           else R(k2)=QSpace(contractQS(I.Psi,I.Psi,'!1*')); end
        else continue; end
     end

     [rr,Ie]=eigQS(R(k2));

     if isfield(I,'rr'), rr=I.rr; end
     if size(rr,2)==2, rr=prod(rr,2); end

     if rr(end)>rr(1), rr=flipud(rr); end
     rr(find(rr<=0))=nan;

     EK(k2)=QSpace(Ie.EK); EK(k2).info.itags=R(k2).info.itags;
     for j=1:numel(EK(k2).data), q=EK(k2).data{j};
        EK(k2).data{j}=fliplr(-log10(q(find(q>1E-14))));
     end

     se(k2)=SEntropy(R(k2));

     h(k)=semilogy(rr,'o-'); hold on

     odir(k)=itags2odir(QSpace(I.itags));
  end
  if largf, fprintf(1,'\r%60s\r',''); end

  setuser(gcf,'Rl',Rl);
  setuser(gcf,'Ql',Ql);
  setuser(gcf,'EK',EK);
  setuser(gcf,'EH',EH);

  se(find(~se))=nan;
  setuser(gcf,'se',se);

  setuser(gcf,'r2',r2);
  setuser(gcf,'dsw',dsw);
  setuser(gcf,'isw',isw);
  setuser(gcf,'ek',ek);
  setuser(gcf,'dd',dd);
  setuser(gcf,'DX',DX);

  kc=find(odir==0); if isempty(kc)
  kc=find(diff(odir)); end
  if isempty(kc), kc=0; else kc=kc(1); end

  hx=plot(nan,nan,'-');
  h(find(h==0))=hx;

  ylim([1E-16,1]); xtight(1.1,'view','x1',0); 
  recolor

  m=floor(L/3); k1=1+m; k2=L-m;
  if all(numel(h)>=[k1 k2])
     set(h(k1),'Disp',sprintf('k=%g',k1));
     set(h(k2),'Disp',sprintf('k=%g',k2));

     h0=h([1:k1-1, k2+1:end]);
     blurl(h0,'LineW',1); mv2back(h0); set(h0,'Marker','none','LineSt','--');

     set(h(k1+1:k2-1),'Marker','none');

     set(h(k1),'Color',[1 0.5 0]); sms(h(k1),4);
     set(h(k2),'Color',[0 0 1]); sms(h(k2),4);

     mv2front(h([k2 k1]));

     legdisp('Location','SouthWest');
  end

ih21=inset('NE','scale',[1.5 1.0]);

  semilogy(r2); xtight(1.1); ylim([1E-10 1E-3]); set(gca,'YTick',logspace(-10,0,6));
  label('bond k','discarded weight');
  postext('SW','\Delta\rho_{disc} \leq %.2e',max(r2));

setax(ah(2,1));

  label('','eig(R_k)','[data for left and right ends of chain blurred and dashed]')

setax(ah(1,1));

  kk=1:length(se);

  plot(kk,se,'o-'); sms(4);
  xtight(1.05,'x1',0);

  if kc>2 && kc<(L-2)
       o={'istr',{sprintf('k_c=%g',kc),0.4,'Rot',90,'FontSize',10,...
          'FontW','norm','BackgroundC','none'}};
  else o={}; end
  xmark(kc,o{:}); % 'k:',

  q=load_dmrg_data(HAM,1); d=getDimQS(q.AK);

  label('bond k','block entropy S(k)');

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

  if 0
  ih=inset(iso{:});

    x=0:size(HAM.hcpow,1)-1;
    h=plot(x, HAM.hcpow,'o'); sms(h,2);
    xtight(1.05); ytight(1.15); label('bond index k','\alpha');
    set(gca,'tag','k2x');

    xd=HAM.hcpow(:,2); 
    if isempty(xstr), xstr='x [hcpow(:,2)]'; end

    if isset('kmin'), m=size(HAM.hcpow,2); amin=zeros(1,m);
       xmark(kmin,'r:');
       k=[max(1,floor(kmin-2)) : min(L,floor(kmin+2)) ]';
       for j=1:m
          p=polyfit(k,HAM.hcpow(k,j),3);
          amin(j)=polyval(p,kmin);
       end
       postext({'N',[ 0.1 0.3]},'\\alpha_{min(S)}=[%s ]',sprintf(' %.5g',amin));
       if numel(amin)==3
          a0=[1, 13/108, 1/216];
          for i=1:numel(a0)
             wblog(' * ','alpha(%g) = %8.5f = %8.5f (=%6s) @ %.3g',...
             i, amin(i),a0(i),num2rat(a0(i)),(amin(i)-a0(i))/a0(i)); 
          end
       end
       setuser(gcf,'amin',amin);
       setuser(gca,'info','hcpow');
    end

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

  if nargout
     Iout=add2struct('-',odir,R,kc,'kmin?','amin?');
  end

  if ksplit>1
     setax(ah(1,1));
     splitline(ksplit); sms(4);
  end

setax(ah(2,2));
  if eoflag, o={'ah',ah(2,[2 2])}; else o={'ah',ah(2,2),'-1'}; end

  [EH_,hH]=getEEdata(EH,o{:},'yl',ylH);
  if ksplit>1, splitline(ksplit); end

  title('energy flow');
  ylabel('eig(H_{L/R}) - E_0');

setax(ah(1,2));
  if eoflag, o={'ah',ah(1,[2 2])}; else o={'ah',ah(1,2),'-1'}; end

  [EE,hh]=getEEdata(EK,o{:},'yl',yl12);
  if ksplit>1, splitline(ksplit); end

  title('entanglement flow');

  label('bond k','entanglement spectrum (using log_{10})');
  xtight(1.05,'x1',0);
  xmark(kc);

  if isset('kmin')
     xmark(kmin,'r:');
     setuser(gcf,'kmin',kmin);
  else
     setuser(gcf,'kmin',[]);
  end

  if ~isempty(kQ), setax(ah(1,2));

    if iscell(kQ), emax=kQ{2}; kQ=kQ{1}; else emax=4; end

    to={'VerticalAl','bottom','HorizontalAl','center','FontSize',8};
    lh=findall(gca,'type','line'); % hh=cat(1,hh{:})';
    eh=[]; qq=[]; qs={};

    dx=diff(xlim)/10; if isempty(dyl), dyl=emax/20; end
    for ih=1:numel(lh), h=lh(ih);
      qi=getuser(h,'Q'); if isempty(qi), continue; end
      x=get(h,'XData');
      y=get(h,'YData');
      j=find(~isnan(y)); if isempty(j), continue; end

         [dk,is]=sort(abs(x(j)-(kQ+0.1))); j=j(is);

         if numel(j)>2 && y(j(2))<y(j(1)), j=j(2); else j=j(1); end

      if y(j)>emax, continue; end

      if ~isempty(qq)
         [i1,i2]=matchIndex(qq,qi);

         if ~isempty(i1)
            if eh(i1,1)<y(j), continue; end
            delete(eh(i1,2:end)); eh(i1,:)=nan;
            k=i1;
         else
            k=size(qq,1)+1;
         end
         l=numel(find(abs(eh(:,1)-y(j))<dyl));
      else
         l=0; k=1;
      end

      qs{k}=['[' sprintf('%g',qi) ']'];

      tx={'Color',get(h,'Color'),'tag','qlabel'};
      h1=plot(x(j),y(j),'o',tx{:}); sms(h1,4);
      h2=text(kQ+l*dx,y(j),qs{k},to{:},tx{:});

      qq(k,:)=qi;
      eh(k,:)=[y(j),h1,h2];
    end

    [~,is]=sort(eh(:,1),'descend');
    eh=eh(is,:); qq=qq(is,:); qs=qs(is);

    if vflag
       for i=1:size(qq,1)
          fprintf(1,'   %10s : %20.16f\n',qs{i},eh(i,1));
       end
    end

  end

  if ~isempty(Ql)
     ql=uniquerows(cat(1,Ql{:})); nl=numel(Ql); sa=zeros(nl,nl);
     for i=1:numel(Ql)
     for j=i+1:numel(Ql)
        sa(i,j)=isequal(Ql{i},Ql{j});
     end
     end
     sa=sa+sa'; ql=sum(sa,2);

     i2=find(ql==max(ql)); qloc=Ql{i2(1)};
     ix=find(ql~=max(ql));

     if ~isempty(ix)
        s=sprintf(...
          'WRN got different local state space for %g site(s)!',numel(ix));
        printfc(['\n\e[31;1m>> ' s '\e[0m \n\n']);

        for i=1:numel(ix)
           fprintf(1,'site %g:\n',ix(i));
           disp(Ql{ix(i)});
        end

        fprintf(1,...
          'local state space otherwise (%g/%g sites):\n',numel(i2),numel(ql));
        disp(Ql{i2(1)});
     end
  end

  setuser(gcf,'ah',ah);

  nd=numel(dsw); d1=zeros(1,nd); d2=0;
  for i=1:nd
     d1(i)=size(dsw{i},1);
     q=max(dsw{i}(:,1)); if d2<q, d2=q; end
  end

  rtol=getopt(HAM.info,'sweep','rtol',[]);
  if ~isempty(rtol)
       sr=sprintf(' @ rtol=%g',rtol);
  else sr=''; end

  if isw<1, i=mean(d1(2:end-1))/2;
  else i=isw; end

  [~,h]=header(['%M :: %.3g sweeps at N_{kept} \leq %g multiplets ' ... 
    '({\\equiv}%g states)%s'],i,d2,D,sr);

  if isfield(HAM.info,'host')
     q=HAM.info.host;
       s=q.host; if iscell(s), s=s{1}; end
       t=datestr(datenum(q.date),'yymmdd');
     header('SE',[s ' :: ' repHome(q.pwd) ' :: Wb' t],{'FontSize',8});
  end

  setax(h);
    s=regexprep(HAM.info.istr,'\n',', ');
    s=regexprep(s,'qloc=\[[\d\s]*\]',...
      sprintf('qloc=[%s]',mat2str2(qloc,'fmt','%g','rowsep','; ')));

  postext([0.01 0.96], s);

if isempty(k2x), return, end
if ~iscell(k2x) || numel(k2x)<2 || ~ischar(k2x{2})
   wbdie('invalid usage (k2x)');
end

  s=k2x{2}; if numel(k2x)>2, k1L=k2x{3}; else k1L=[1 L]; end
  k2x=k2x{1};

  for ih=[1 3 4], setax(ah(ih))
     if numel(k2x)==2
      % interpreted has mapping k=[1 L] to k2x([1 2])
      % for h=hh
      %    x=get(h,'XData');
      %    x=k2x(1)+(x-k1L(1))*(diff(k2x)/diff(k1L));
      %    set(h,'XData',x);
      % end
        map_xdata(gca,'lin',[1 L]-1,k2x);
     else
        hh=findall(gca,'type','line')';
        for h=hh
           x=get(h,'XData'); j=find(~isnan(x));
           if ~isempty(j), x(j)=k2x(x(j)); set(h,'XData',x); end
        end
     end

     xlabel(['bond index k  \rightarrow  ' s]);
     xtight(1.05);
  end

end

