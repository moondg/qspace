function plotHSS(HAM,varargin)
% function plotHSS(HAM [,opts])
% Wb,Jan22,16

% plotHSS(HAM,'xl',[101.5 107.5],'--lubox','sfac',[1 2 1],'--rbox2')
% plotHSS(HAM,'xl',[101.5 107.5],'--lubox','-k','--rbox2')

% adapted from plotHAM.m

  getopt('init',varargin);
     isw  = getopt('isw',[]);
     sfac = getopt('sfac',1);
     rfac = getopt('rfac',0);
     tflag= getopt('-t');
     xdk  = getopt('xdk',[]);
     nofix= getopt('-nofix');
     ah   = getopt('ah',[]);
     xl   = getopt('xl',[]);

     io   = getopt('io',[]);

     bflag= getopt('-b');

     if getopt('--rbox')
          sh={'shape','rbox','--noedge'};
     elseif getopt('--rbox2')
          sh={'shape','rbox2','--noedge'};
     else sh={}; end

     lubox= getopt('--lubox');
  varargin=getopt('get_remaining');
  if isempty(varargin), varargin={'-pn'}; end

  if isempty(HAM) || ~isfield(HAM.info,'HSS')
     wberr('invalid usage (missing field HAM.info.HSS)'); end
  L=numel(HAM.mpo);

  sym=HAM.info.IS.E.info.qtype;
  if regexp(sym,'^SU\d+$')
       N=str2num(sym(3:end));
  else N=0; end

  if isfield(HAM.info,'XY')
     xy=HAM.info.XY;
  else xy=[]; end

  if isempty(xy)
     xy=[(1:L)',repmat(0.1,L,1)]; xy(1:2:end,2)=-0.1;
  end

  if isempty(isw)
  isw=numel(HAM.info.HSS); end

  HSS=HAM.info.HSS{isw};
  if ~isreal(HSS)
     e=[norm(imag(HSS),'fro'), norm(real(HSS),'fro')];
     if e(1)/e(2)>1E12, wblog('WRN','got imag(HSS) @ %.3g %.3g!?',e); end
     HSS=real(HSS);
  end

  if ~bflag
     plot_simple_line_plot(HAM,HSS,ah);
     return
  end

  if size(HSS,3)>1
     if isempty(io), io=1;
        wblog(' * ','using io=%g/%g',io,size(HSS,3));
     end
     HSS=HSS(:,:,io);
  end

  if ~nofix

     [k1,k2,ss]=find(HSS); kk=[k1,k2];
     i=find(k1<k2); kk1=kk(i,:); s1=ss(i);
     i=find(k1>k2); kk2=kk(i,:); s2=ss(i);

     [i1,i2,I]=matchIndex(kk1,fliplr(kk2));
     if ~isempty(I.ix2)
        kx=fliplr(kk2(I.ix2,:)); sx=s2(I.ix2);
        for i=1:size(kx,1), k1=kx(i,1); k2=kx(i,2);
           wblog(' * ','quick fix for bonds connecting leftward of site %g',k2);
           HSS(:,k2)=HSS(k2,:);
        end
     end
  end

  HSS=triu(HSS);

  [k1,k2,ss]=find(HSS);

  dk=abs(k1-k2);
  if ~isempty(xdk)
     ix=matchIndex(dk,unique(xdk(:)));
     k1(ix)=[]; k2(ix)=[]; ss(ix)=[]; dk(ix)=[];
  end

  if ~isequal(sfac,1)
     sfac(end+1:max(dk)+1)=sfac(end);
     ss=ss.*reshape(sfac(dk+1),[],1);
  end

  smean=mean(abs(ss)); 
  if N
       ss=ss*(3/(0.5+1/(2*N)));
  else ss=ss*(1/smean);
  end

  if isempty(ah)
     ah=smaxis(1,1,'tag',mfilename);
     header('%M :: %s',HAM.info.istr); addt2fig Wb
  end; setax(ah(1))

  for i=1:numel(ss)
     x=xy([k1(i),k2(i)],:);
     if rfac, x=x+rfac*randn(size(x)); end

     if ~tflag
        plotbond(x(:,1),x(:,2),ss(i),sh{:},varargin{:});
     else
        plot(x(:,1),x(:,2),'Color',getcolor(1+floor(100*rand)));
     end

     hold on
  end

  title(sprintf('<S.S> \\sim %.4g',smean));
  if ~isempty(xl), xlim(xl); else xl=xlim; end

if ~lubox, return; end

  axis equal
  xlim(xl);

  i=find(k2-k1==2 & xy(k2,2)>0); x=xy(k1(i),1);
  i=i(find(x>=xl(1)-1 & x<=xl(2)+1));

  s=ss(i); s=sign(s-mean(s));
  I={ find(s<0), find(s>0) };
  n=[ numel(I{1}), numel(I{2}) ];

  if n(1)>=n(2)
       j=find(s>0);
  else j=find(s<0);
  end

  o={'EdgeC','none','FaceC',[1 .8 .3]};

  xb=xy(k1(i(j)),1); % xmark(x+0.5,'k:');
  for i=2:numel(xb)
     x=[xb(i-1)+1, xb(i)]; x=x([1 2 2 1 1]);
     y=[0 0 1 1 0];
     patch(x,y,repmat(-1,size(x)),'y',o{:});
  end

end

% -------------------------------------------------------------------- %
function plot_simple_line_plot(HAM,HSS,ah)

   sH=size(HSS); sH(end+1:3)=1;

   if isempty(ah)
      ah=smaxis(sH(3),1,'tag',mfilename);
      header('%M :: %s',HAM.info.istr); addt2fig Wb
   end

   n=max(sH(1:2)); sm={'o','+'};
   for k=1:sH(3), setax(ah(k));
      for i=0:n, c=getcolor(i+1);
         for l=[2 1]
            if l==1, id=i; elseif l>0, id=-i; else continue; end
            yd=diag(HSS(:,:,k),id); if norm(yd)<1E-14, continue; end

            xd=(1:(n-i))+i/2;
            h=plot(xd,yd,'-','Color',c); blurl(h,'LineW',1); hold on
            h=plot(xd,yd,sm{l},'Color',c); sms(h,2);
            if l==2, blurl(h,'LineW',1,'cfac',0.8); end
            set(h,'Disp',sprintf('i=%g',id))
         end
      end

      mv2back(findall(gca,'Type','line','LineSt','-'));
      legdisp
   end

end

% -------------------------------------------------------------------- %

