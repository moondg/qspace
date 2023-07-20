function Iout=mp(M,varargin)
% Function I=mp(M [,opts]))
%
%    matrix plot:
%    simple wrapper to generate surface plot of a real-valued 2D matrix
%
% Options
%
%     '-q'              quiete mode
%     'xd',...          x-axis data
%     'yd',...          y-axis data
%
%     'diag', x         set diagonal to value x (eg NAN or 0, ...)
%     'trim', k         trim k of cols/rows all around
%     'tstr', ...       intro to title string ('z = ')
%
%     'zlim', [z1 z2]   ZLim for surface plot (default: 'auto')
%     'clim', [c1 c2]   CLim for surface plot (default: 'auto')
%     'view', [al be]   view for surface plot (default: 2D view)
%     '-3D'             3D view
%
%     '-smooth'         interpolate color plot (using shading interp)
%     '-cm'             use colormap cmap0w (exclusive to -lh)
%     '-lh'             use lighting (light2) (exclusive to -cm)
%     '-xy'             use axis xy instead of axis ij
%     '~res'            do not down-sample matrix (done by default)
%     '-img'            convert matrix to image and plot as image
%     '-log[xyz]'       log-scale for x-, y-, or z-axis
%
% Wb,Nov26,02 ; Wb,Jun24,08 ; Wb,Mar14,19

% depreated
%     '-gca'            plot into current axes set (same as mp0) // Wb,Mar21,19
% for last treatment of 3D data, see Archive/mp_190314.m

  if nargin<1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  [n1,n2,n3]=size(M); olg={};

     if n3>1,        wblog('ERR','got 3D array');     return; end
     if isvector(M), wblog('ERR','got vector');       return; end
     if ~isreal(M)
        wblog('WRN','got complex data (skipping imag)');
        M=real(M);
     end

  lxyz=zeros(1,3);

  getopt('init',varargin);
     xd=getopt('xd',1:n2); xd=reshape(xd,1,[]);
     yd=getopt('yd',1:n1); yd=reshape(yd,1,[]);

     if getopt('-logx'), lxyz(1)=1; olg{end+1}='XScale'; end
     if getopt('-logy'), lxyz(2)=1; olg{end+1}='YScale'; end
     if getopt('-logz'), lxyz(3)=1; olg{end+1}='ZScale'; end

     v3=getopt('-3D'); if ~v3
     v3=getopt('view',[]); end

     qflag=getopt('-q');

     x=getopt('diag',[]); if ~isempty(x), M=M-diag((diag(M)-x)); end
     trm=getopt('trim',0);

     img=getopt('-img'); if ~img
       smofig=getopt('-smooth');
       nores=getopt('~res');
       zl=getopt('zlim',[]);
     end
     xyflag=getopt('-xy');

     cmap=getopt('-cm'); ocm={};
     if ~cmap
        ocm=getopt('cm',{});
        if ~isempty(ocm), cmap=2; end
     end

     if ~cmap
        lght=getopt('-lh');
        cl  =getopt('clim',[]);
     else lght=0; cl=[]; end

     if ~isempty(v3) && v3, addrc=0;
     else addrc=~getopt('~addrc');
     end

     tstr=getopt('tstr','z =');

     if getopt('-gca')
        wblog('WRN','option mp(..,''-gca'') is deprecated [03/2019]');
        dispstack(dbstack);
     end

  getopt('check_error');

  if addrc
     M(end+1,:)=2*M(end,:)-M(end-1,:);
     M(:,end+1)=2*M(:,end)-M(:,end-1);

     if ~isempty(xd), xd(end+1)=2*xd(end)-xd(end-1); end
     if ~isempty(yd), yd(end+1)=2*yd(end)-yd(end-1); end
  end
  if trm, i=trm;
     M=M(i+1:end-i, i+1:end-i);
     if ~isempty(xd), xd=xd(i+1:end-i); end
     if ~isempty(yd), yd=yd(i+1:end-i); end
  end

  ah=gca; zr=full([ min(M(:)), max(M(:)) ]);

  if ~img
     if ~nores, s=get(0,'ScreenSize'); if s(3:end)>1
        f=get(gcf,'Position');
        a=get(ah,'Position'); 

        if all(a<=1) && all(f(3:end)>1), m=size(M); m_=m;
           f=fliplr(a(3:4).*f(3:4));
           p=round(m./f - 0.25);
           for i=find(p>1)
              if i==1, q=yd; else q=xd; M=M'; m=size(M); end

              l=mod(-m(1),p(i)); if l
                 M(end+1:end+l,:)=repmat(M(end,:),l,1); m=size(M);
              end

              s=[p(i), m(1)/p(i), m(2)];
              M=permute(mean(reshape(M,s),1),[2 3 1]);

              if ~isempty(q)
                 if lxyz(i), q=log(q); end
                 if l, q(end+1:end+l)=q(end); end
                 q=permute(mean(reshape(q,s(1:2)),1),[2 1]);
                 if lxyz(i), q=exp(q); end
              end

              if i==1, yd=q; else xd=q'; M=M'; end
           end
           if ~qflag, m=size(M);
              if ~isequal(m,m_), wblog(' * ',...
                 'downsampled %gx%g => %gx%g [ah @ %gx%g]',m_,m,round(f));
              end
           end
        end
     end, end

     h=surface(xd,yd,M);

     shading flat  ; if smofig
     shading interp; end

     if ~isempty(zl), set(ah,'ZLim',zl); end

     if numel(v3)==1, if v3, view(3); end
     elseif numel(v3)==2 && isnumeric(v3), view(v3);
     elseif ~isempty(v3), wblog('ERR','invalid 3D view specs'); end

  else
     h=image(xd,yd,M,'CDataMapping','scaled');
  end

  if ~isempty(cl), set(ah,'CLim',cl); end

  if xyflag, axis xy, else axis ij; end
  xytight;  box on

  if lght, light2; elseif cmap, cmap0w('ah',gca',ocm{:})
     if isempty(cl) && img
        q=max(abs(zr)); set(ah,'CLim',[-q q]);
     end
  end

  if ~isempty(olg), olg(2,:)={'log'}; set(ah,olg{:}); end

  q=[zr, diff(zr)]; i=find(q);
  qfac=max(floor(log10(abs(q(i)))));

  if abs(qfac)>3, q=q/(10^qfac);
       tstr=sprintf('%s [%.2g .. %.2g @ %.2g]*10^{%d}',tstr,q,qfac);
  else tstr=sprintf('%s [%.2g .. %.2g @ %.2g]',tstr,q);
  end

  label3('x','y','',tstr);

  if nargout,
     Iout=add2struct('-',h,q,v3,xyflag);
  end

end

