function plotSE(HAM,IL,varargin)
% function plotSE(HAM,IL,varargin)
% Wb,Nov22,17

  getopt('init',varargin);
     dk=getopt('dk',2);
     perBC=getopt('-pBC');
     yfac=getopt('yfac',1); if numel(yfac)==1, yfac=[-0.2 1]*yfac; end
     wsc =getopt('wsc',1);
     xflag=~getopt('-noex');
  varargin=getopt('get_remaining');

  [se,Ie,IC]=getSE(HAM,IL,varargin{:});

  SE=Ie.se;
  if xflag, SE(:,:,end+1)=se; end

  niter=size(SE,3); n0=ceil(niter/2);

  L=size(SE,1)+1;
  x=(L/pi)*sin(pi*(1:L-1)'/L);
  lx=log(x); xl=max(lx); x_=pi/L; x4=0.25;

  ell=getsymb('\ell');
  sll=getsymb('\ll');
  seq=getsymb('\simeq');

  if wsc<=0,
       xlb=[ 'log '  ell];  x_=-log(x_); x4=-log(x4);
  elseif wsc==2
       xlb=['1/log ' ell];  x_=-1/log(x_); x4=-1/log(x4);
  elseif wsc==3
       px=1/2; x_=x_^px; x4=x4^px;
       xlb=[ '1/' ell sprintf('^{%s}',rat2(px))];
  elseif wsc==4
       px=1/3; x_=x_^px; x4=x4^px;
       xlb=[ '1/' ell sprintf('^{%s}',rat2(px))];
  else xlb=['1/'     ell];
  end

  x3=3*x_; 

  if perBC, cfac=3; else cfac=6; end

ah=smaxis(dk,2,'tag',mfilename); addt2fig Wb

  s=['mat_to_headerSE; header(''' sprintf('%s :: ',mfilename) ' %s'',istr2)'];
  evalin('caller',s);

  h=header('hleft');
     Dk=max(max(getdatafield(IL(:,:,end),'Dk',1)));
     i=size(IL(1).Dk,1);
     if i>1, Dk(2)=max(max(getdatafield(IL(:,:,end),'Dk',i))); end
  set(h,'String',sprintf('%s\nL=%g, D^\\ast \\leq %s (%s)',...
     get(h,'String'), numel(HAM.mpo), ...
     int2str2(Dk(1)), int2str2(Dk(end)) ...
  ));

  h=header('fleft',[ ell ' \equiv (L/\pi) sin({\pi}x/L) ' seq ' x for x' sll 'L']);

  ox={'Color','r'};

  for k=1:dk
     ik=[ k:dk:(L/2-1E-8) fliplr(L-k:-dk:(L/2-1E-8)) ];

     yk=permute(SE(ik,:,:),[1 3 2]);
     li=lx(ik);

   setax(ah(k,1));

     h=plot(li,yk(:,:,1));
       if xflag, set(h(end),ox{:}); end; blurl(h,'LineW',1);
       set(gca,'ColorOrderIndex',1); hold on; 
     h=plot(li,yk(:,:,2));
       if xflag, set(h(end),ox{:}); end; blurl(h(1:n0),'LineW',1);

     label(['log ' ell],'entropy S',sprintf('k=%g,%g,\\ldots,L-%g',k,k+dk,k));
     xlim([0,xl]);

     n=numel(h);
     for i=n-3:n
        if i<n || ~xflag, s=sprintf('sweep %g',i); else s='extrap.'; end
        set(h(i),'Disp',s);
     end

     legdisp('Location','NW');

     if k==dk
        xm=0.20;
        lm=(L/pi)*sin(pi*xm);
        om={'k','Color',[.8 .8 .8],'istr',...
           {sprintf('L/%g',1./xm),'FontW','normal','Margin',0.01}};
        xmark(log(lm),om{:});
     else lm=[];
     end

   setax(ah(k,2));

     dx=diff(li); xi=x(ik);
     x1=0.5*(xi(2:end)+xi(1:end-1));
     dy=cfac*diff(yk,[],1) ./ repmat(dx,[1,size(yk,2),2]);

     if     wsc<=0, x1=log(x1);    lm=log(lm);
     elseif wsc==2, x1=1./log(x1); lm=1./log(lm);
     elseif wsc>=3, x1=x1.^(-px);  lm=lm.^(-px);
     else           x1=1./x1;      lm=1./lm;
     end

     ix=find(abs(dx)<1E-10);
     dy(ix,:,:)=nan;

     h=plot(x1,dy(:,:,1));
       if xflag, set(h(end),ox{:}); end; blurl(h,'LineW',1);
       set(gca,'ColorOrderIndex',1); hold on; 
     h=plot(x1,dy(:,:,2));
       if xflag, set(h(end),ox{:}); end; blurl(h(1:n0),'LineW',1);

     label(xlb,['c \sim ' num2str(cfac) ' dS / d ln ' ell]);

     if wsc>0
          ik=find(x1>2*x_);
     else ik=find(x1<x_/2);
     end
     q=reshape(dy(ik,:,:),[],1); q=[min(q),max(q),std(q)];
     ylim(q(1:2)+q(3)*[-0 1]);
     if wsc>0
          xlim([0.5*x_,x4]);
     else xlim([x4,1.2*x_]);
     end

     h=xmark([x_,x3],'k','Color',[.8 .8 .8]);
     mv2back(h);

     if k==dk, xmark(lm,om{1:end-2}); end

   ah(k,3)=repax(inset({'SE',[0 0.1]},'scale',[1.6 1.2]),ah(k,2));

     if wsc>0
          i=find(x1>1.2*x_ & x1<1.5*x_);
     else i=find(x1<0.8*x_ & x1>0.6*x_);
     end

     if norm(diff(sort(x1(i))))<0.2
        i=unique([i,i+1]); i(find(i>length(x1)))=[];
     end

     q=reshape(dy(i,end-3:end,:),[],1); q=[min(q),max(q),std(q)];

     xlim([0.*x_,x3]);
     ylim(q(1:2)+q(3)*yfac);
     xmark(x_,'k:');

     label('','','');
     set(gca,'FontSize',9);
  end

  setuser(gcf,'ah',ah);

end

