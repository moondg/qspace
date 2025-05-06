function gg=figdata(varargin)
% Function gg=figdata([opts])
%
%   Get data from existing figure which exists as an image
%   e.g. extracted from some pdf / publication.
%
% Options
%
%   '-x'    (re)specify x-range
%   '-y'    (re)specify y-range
%
%   '--ax'  (re)specify x-range where ticks are located automatically (bottom)
%   '--axt' same as -ax, but look for ticks at top of axis
%
%   '--ay'  (re)specify y-range where ticks are located automatically (left)
%   '--ayr' same as -ay, but look for ticks at righ side of axis
%
% Wb,Nov30,08 ; Wb,Jun28,18

% tags: pdf, png, xcoord, ycoord, points

  getopt('init',varargin);

     if getopt('-xy'), xflag=1; yflag=1;
     elseif getopt('-axy'), xflag=2; yflag=2;
     else
        if getopt('-x'), xflag=1;
        elseif getopt('--ax'), xflag=2;
        elseif getopt('--axt'), xflag=3; else xflag=0; end

        if getopt('-y'), yflag=1;
        elseif getopt('--ay'), yflag=2;
        elseif getopt('--ayt'), yflag=3; else yflag=0; end
     end

  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg==1
    if ischar(varargin{1})
        f=varargin{1}; varargin(1)=[];
        if ~exist(f,'file'), wblog('ERR',...
        'no image file %s available',f); end

        fh=findall(groot,'type','figure','tag',mfilename);
        if ~isempty(fh), figure(fh); ah=getuser(fh,'ah'); else
           ah=smaxis(1,1,'tag',mfilename); setax(ah(1));
           A=imread(f); image(flipud(A)); hold on; axis xy
           xflag=2; yflag=2;
        end

     elseif isscalar(varargin{1}) && isaxis(varargin{1})
        ah=varargin{1}; varargin(1)=[];
        setax(ah);
     end
  end

  xsc=getuser(gcf,'xsc');
  if xflag || isempty(xsc) || ~isnumeric(xsc) || numel(xsc)~=3
     if xflag<=1
          xsc=get_scale('x');
     else xsc=get_scale_auto(0,xflag);
     end
     setuser(gcf,'xsc',xsc);
  end

  ysc=getuser(gcf,'ysc');
  if yflag || isempty(ysc) || ~isnumeric(ysc) || numel(ysc)~=3
     if yflag<=1
          ysc=get_scale('y'); 
     else ysc=get_scale_auto(1,yflag); 
     end
     setuser(gcf,'ysc',ysc);
  end

  inl(1); k=0; gg={[]};

  while 1
     fprintf(1,['   collect xy-data (%d) by clicking to the figure '... 
     '(stop by pressing enter)...\n'],k);
     g=ginput; if isempty(g), inl(1); break; end

     k=k+1; gg{k}=g;
  end

  if length(gg)==1, gg=gg{1}; end

end

% -------------------------------------------------------------------- %

function sc=get_scale(x) 

  while 1
     fprintf(1,...
     '   Determine %c-scale by clicking range [%c1 %c2] ...\n',x,x,x);
     figure(gcf); g=ginput;
     if size(g,1)==2, break; else
     fprintf(1,'   ERR exactly two points required!'); end
  end

  while 1
       i=input(sprintf(...
       '   Enter values [%c1 %c2] (including brackets): ',x,x));

     if isnumeric(i) && numel(i)==2, break; else
     fprintf(1,'   ERR two values required!'); end
  end
  inl 1

  if x=='x'
       sc=scale_x(g(:,1),i);
  else sc=scale_y(g(:,2),i);
  end

end

% -------------------------------------------------------------------- %

function sc=scale_x(x0,xnew)

  sc=[ x0(1), diff(xnew)/diff(x0), xnew(1) ];

  xl=get(gca,'XLim');
  xl= (xl-sc(1))*sc(2) + sc(3);

  for h=[ findall(gca,'type','line'); findall(gca,'type','image') ]'
     x=get(h,'XData'); x=(x-sc(1))*sc(2) + sc(3);
     set(h,'XData',x);
  end
  xlim(xl);

end

% -------------------------------------------------------------------- %

function sc=scale_y(y0,ynew)

  sc=[ y0(1), diff(ynew)/diff(y0), ynew(1) ];

  yl=get(gca,'YLim');
  yl= (yl-sc(1))*sc(2) + sc(3);

  for h=[ findall(gca,'type','line'); findall(gca,'type','image') ]'
     y=get(h,'YData'); y=(y-sc(1))*sc(2) + sc(3);
     set(h,'YData',y);
  end

  set(gca,'YLim',sort(yl),'YDir','normal');

end

% -------------------------------------------------------------------- %
% yflag: whether to operate on y-axis (default: x-axis)
% adapted from test routine figscale.m
% Wb,Jun28,18

function sc=get_scale_auto(yflag,flag) 

   ot={'tag','fig:autoscale'};
   h=findall(gca,'Type','Imag');
   if numel(h)~=1
      wbdie('failed to identify image data (e=%g)',numel(h)-1); end
   A=get(h,'CData');

   A=A/max(abs(A(:)));
   if ndims(A)==3, A=sqrt(sum(A.^2,3)); end

   [ny,nx]=size(A);
   nx2=floor(nx/2);
   ny2=floor(ny/2);

   xd=get(h,'XData'); xd=linspace(xd(1),xd(2),nx); dX=diff(xd(1:2));
   yd=get(h,'YData'); yd=linspace(yd(1),yd(2),ny); dY=diff(yd(1:2));

   xx=sum(A,1); dx=diff(xx);

      i=1:nx2;
      ipl=find(dx(i)==max(dx(i)));
      inl=find(dx(i)==min(dx(i))); il=[ipl,inl];

      i=nx2+1:nx-1;
      ipr=find(dx(i)==max(dx(i)));
      inr=find(dx(i)==min(dx(i))); ir=nx2+[ipr,inr];

   yy=sum(A,2); dy=diff(yy);

      j=1:ny2;
      jpu=find(dy(j)==max(dy(j)));
      jnu=find(dy(j)==min(dy(j))); jd=[jpu,jnu];

      j=ny2+1:ny-1;
      jpd=find(dy(j)==max(dy(j)));
      jnd=find(dy(j)==min(dy(j))); ju=ny2+[jpd,jnd];

   if numel(il)~=2 || numel(ir)~=2 || numel(ju)~=2 || numel(jd)~=2
      wbdie('failed to identify axis boundaries');
   end

   xr=[mean(xd(il)), mean(xd(ir))];
   yr=[mean(yd(ju)), mean(yd(jd))];

   h=plot3(xr([1 1 2 2 1]), yr([1 2 2 1 1]), ones(1,5),'c--',ot{:});

   if yflag
      if flag<3, q=il; else q=ir; end
      i=mean(q)+sort([-3 3]*diff(q));
      i=floor(i(1)):ceil(i(2));

      y1=sum(A(:,i),2);
      dy=diff(y1);

      q=0.3*max(abs(dy)); ip=find(dy>+q); in=find(dy<-q);
      if numel(ip)==numel(in)
         ym=mean([yd(ip); yd(in)],1)+0.5*dY;
         hm=plot3(xd(i([1 end])),[ym;ym],[1 1],'m',ot{:});
      end

      for i=1:min(10,numel(ym))
         fprintf(1,'%4g.  %8g\n',i,ym(i)); end
      s='y';
   else
      if flag<3, q=jd; else q=ju; end
      j=mean(q)+sort([-3 3]*diff(q));
      j=floor(j(1)):ceil(j(2));

      x1=sum(A(j,:),1);
      dx=diff(x1);

      q=0.3*max(abs(dx)); ip=find(dx>+q); in=find(dx<-q);
         ip(find(diff(ip)==1))=[];
         in(find(diff(in)==1))=[];

      if numel(ip)==numel(in)
         xm=mean([xd(ip); xd(in)],1)+0.5*dX;
         hm=plot3([xm;xm],yd(j([1 end])),[1 1],'m',ot{:});
      else
         wbdie('failed to determine x-ticks (try other top/bottom axis?)');
      end

      for i=1:min(10,numel(xm))
         fprintf(1,'%4g.  %8g\n',i,xm(i)); end
      s='x';
   end

   while 1
      q=input(sprintf(...
      '\n   Select pair of %s-ticks [-1 is last, k(eyboard)]: ',s),'s');
      if isequal(q,'k'), keyboard, continue; end
      eval(['i=[' q '];']);
      if isnumeric(i) && numel(i)==2, break; else
      fprintf(1,'   ERR two values required!'); end
   end

   m=numel(hm); j=find(i<0);
   if ~isempty(j)
     i(j)=m+i(j)+1;
   end

   set(hm(i),'LineW',3); grid on

   while 1
      q=input(sprintf(...
      '   Specify %s-range: ',s),'s');
      eval(['q=[' q '];']);
      if isnumeric(q) && numel(q)==2, break; else
      fprintf(1,'   ERR two values required!'); end
   end

   if yflag
        sc=scale_y(ym(i),q);
   else sc=scale_x(xm(i),q);
   end

   delete(findall(gca,ot{:}));

end

% -------------------------------------------------------------------- %

