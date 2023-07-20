function ah=breakaxes(ah_,varargin)
% function ah=breakaxes(ah[, opts])
%
%    specify (at least) two axes handles
%    to overlay with new axes set that shows wiggly separator
%
% Wb,Oct24,20

  if ~all(isaxis(ah_))
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end
  nah=numel(ah_);
  if ~nah
     wbdie('got %g axis handles (at least 2 required)',nah);
     return
  end

  setfig(get(ah_(1),'parent'));
  fh=gcf; ah=findall(gcf,'tag',mfilename);

  getopt('init',varargin);
     ph=getopt('ph',30);
     w =getopt('w',[]);
     nx=getopt('nx',32);
     tflag=getopt('-t');
     sfac =getopt('sfac',1);

     if getopt('--rm'), delete(ah); ah=[]; end
  rpos=getopt('get_last',0.5);

  pp=get(ah_,'Position'); if iscell(pp), pp=cat(1,pp{:}); end
  if isempty(w) && nah>1
     dp=max(pp(:,1:2),[],1) - min(pp(:,1:2),[],1);
     w=find(dp==max(dp),1);
  end

  p=[min(pp(:,1:2),[],1), max(pp(:,1:2)+pp(:,3:4),[],1)];
  p=[ p(1:2), p(3:4)-p(1:2) ];

  if nah==1
     q={isempty(w), isempty(rpos)};
     if ~xor(q{:}), wbdie(['invalid usage ' ... 
        '(need unique specification of location for single axis)']); end
     if q{1}, q=rpos; else q=w; end
     if ~iscell(q) || numel(q)~=2, q, wbdie('invalid usage'); end
     w=0; rpos=[];

     if     isequal(q{1},'S'), w=2; p(2)=p(2)-0.5*p(4);
     elseif isequal(q{1},'N'), w=2; p(2)=p(2)+0.5*p(4);
     end
     if w, rpos=q{2};
        if     rpos>.7, p(1)=p(1)+0.5*p(3); rpos=rpos-.5;
        elseif rpos<.3, p(1)=p(1)-0.5*p(3); rpos=rpos+.5; end

     elseif isequal(q{1},'E'), w=1; p(1)=p(1)+0.5*p(3);
     elseif isequal(q{1},'W'), w=1; p(1)=p(1)-0.5*p(3);
     else w, rpos, wbdie('invalid value for location'); end

     if isempty(rpos), rpos=q{2};
        if     rpos>.7, p(2)=p(2)+0.5*p(4); rpos=rpos-.5;
        elseif rpos<.3, p(2)=p(2)-0.5*p(4); rpos=rpos+.5; end
     end
     if ~w || isempty(rpos)
         w, rpos, wbdie('invalid value for location'); end

     if tflag, wblog('TST','using w=%g, rpos=%g',w,rpos); end
  elseif tflag, wblog('TST','using w=%g',w); 
  end

  if ~isempty(ah), q=[];
      p0=get(ah,'pos'); if ~iscell(p0), p0={p0}; end
      for i=1:numel(p0), q(i)=norm(p0{i}-p); end
      e=min(q);
      if e<0.1*norm(p(3:4)), ah=ah(find(q==e),1); cla
         if tflag, wblog('TST','using existing breakaxes'); end
      else ah=[]; end
  end

  if isempty(ah), ah=axes('Position',p); end
  if tflag
     set(ah,'XColor','r','YColor','r','Color','none');
     box(ah,'on'); axis on
  else axis off; end

  p=get(gca,'Pos'); r=p(4)/p(3); r2=r/3; dy=1.5*r2;

  setuser(ah,'ah',ah_);

  axis([[0 2],[0 2]]); hold on

  ph=ph*(pi/180); if w>1, ph=ph+pi/2; end

  sfac=sfac*0.025; if numel(sfac)==1, sfac(2)=sfac; end

  R=[cos(ph) -sin(ph); sin(ph) cos(ph)];
  R=diag([1 1/r].*sfac)*R*diag([1,r]);

  q=[1;1]; if w==1, q(2)=2*rpos; else q(1)=2*rpos; end
  q=repmat(q,1,nx);

  x=linspace(-1,1,nx);
  y=-r2*sin(x*pi);
  xy={ q + R*[x; y+dy]; q + R*[x; y-dy] };
  XY=[xy{1},fliplr(xy{2})]; n=size(XY,2);

  h=patch(XY(1,:),XY(2,:),'b','EdgeColor','none','FaceColor',[1 1 1]-1E-4);

  plot(xy{1}(1,:),xy{1}(2,:),'k');
  plot(xy{2}(1,:),xy{2}(2,:),'k');

  set(ah,'tag',mfilename);

  if ~nargout, clear ah; end
end

