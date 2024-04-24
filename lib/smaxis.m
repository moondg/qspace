function ah = smaxis (dim1, dim2, varargin)
% FUNCTION: smaxis - setup multiple axis system
% USAGE:
%
%   ah = smaxis (dim1, dim2 [OPTIONS...]))
%   where dim1, dim2 - specifies number of axis to set
%
% OPTIONS:
%
%  'tag',..    tag for axis set (also indicateds to reuse axis set!)
%  'dx',...    internal margin in x
%  'dy',...    internal margin in y
%  'mrg',..    margin (outer and default for interior)
%  'X0',...    shift in x (1.5*mar)
%  'Y0',...    shift in y (1.5*mar)
%  'DX',...    plain linear shift (without rescaling)
%  'DY',...    plain linear shift (without rescaling)
%  'pos',..    explicit position array
%  'fpos',..   explicit figure position
%  '-fpos2'    square figure position
%  'zoom',.    global zoom *inside* figure window
%  'ZM',..     zoom of entire figure window
%  'tight'     axis sets placed tight
%  'landscape' orientation
%  'current'   take current figure
%
% set axis to active by:  axes(ah(i))
%
% (C) Wb,Feb12,01
% ------------------------------------------------------------------

% tags: title, window name
% => set by 'Name' field (see s=get(gcf,'Name') below).
% Wb,Apr15,13

  if nargin<2
     helpthis
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  getopt ('INIT', varargin);
     tight  = getopt ('tight');
     lands  = getopt('landscape');
     fpos   = getopt ('fpos',[]); if isempty(fpos)
        if getopt('-fpos2'), fpos=[900 410 710 720]; end
     end
     pos2   = getopt('pos2');
     m=getopt ('mrg',0.08);

     if lands
          my=1.5*m; x0=1.5*m; y0=1.2*m; DY=0.5*m;
     else my=    m; x0=1.5*m; y0=    m; DY=0;
     end

     if pos2
        fpos=[ 956 554 640 570 ];
     end

     cfflag = getopt ('-cf');
     DX     = getopt ('DX',  0);
     DY     = getopt ('DY', DY);
     mx     = getopt ('dx',  m);
     my     = getopt ('dy', my);
     tag    = getopt ('tag','');
     pos    = getopt ('pos',[]);
     zm     = getopt ('zoom',[]);
     ZM     = getopt ('ZM', 1);

     X0=getopt('X0',0); X1=getopt('X1',1-DX);
     Y0=getopt('Y0',0); Y1=getopt('Y1',1-DY);

     if getopt ('-tag'), ftag=1;
     elseif getopt ('-ftag'), ftag=2; % 'full' figure tage (including path)
     else ftag=0; end

  getopt('check_error');

  ss=get(0,'ScreenSize');

  sx=max(ss(3),1024);
  sy=max(ss(4),648);

  r=0.80;

  if tight, mx=0; my=0; end

  if lands
     fx=ZM*sx*0.50;
     fy=   fx*0.72;
  else
     fy=ZM*sy*0.736;
     fx=   fy*r;
  end

  fx=round(fx); fy=round(fy);

  sizex = ((X1-X0) - (0.8*m + x0 + mx*(dim2-1)) ) / dim2;
  sizey = ((Y1-Y0) - (1.0*m + y0 + my*(dim1-1)) ) / dim1;

% [ X1-X0, Y1-Y0; sizex, sizey ]
% wblog('TST','X: %g :: (%.3g+%.3g)*%.3g = %.3g [%+.3g]',...
%     X1-X0,sizex,mx,dim2,(sizex+mx)*dim2, DX+x0+X0 -(sizex+mx));
% wblog('TST','Y: %g :: (%.3g+%.3g)*%.3g = %.3g [%+.3g]',...
%     Y1-Y0,sizey,my,dim1,(sizey+my)*dim1, DY+1-y0+Y0 + my);

  if sizex<=0 || sizey<=0, error('Wb:ERR',...
    'neg. dimensions (invalid size specs, sizex=%g, sizey=%g',...
     sizex, sizey);
  end

% -------------------------------------------------------------------- %
  if isempty(fpos)
     if dim1>1 || dim2>1 || size(pos,1)>1
          fpos = [sx-fx-4, sy-fy-72, fx, fy];
     else fpos = get(0, 'DefaultFigurePosition');
     end
  end

  if ~isempty(tag)
      ff = findall(groot,'type','figure','tag',tag);
      if ~isempty(ff)
         set(groot,'CurrentFigure',ff(1));
         clf
      else
         figure('Pos', fpos);
      end

      set(gcf,'tag',tag);

    % [i,s]=system('hostid.pl'); // Wb,Apr24,24
      set(gcf,'NumberTitle','off');
      s=sprintf('Fig.%g/%s',get(gcf,'Number'),hostid);

      if isempty(findstr(s,tag))
         p=lower(dec2base(getpid,36));
         set(gcf,'Name', [s sprintf('[%s]',p) '/' tag]);
      end
  elseif ~cfflag
      figure('Pos',fpos);
  end

  cla; clf;

  if lands
  set(gcf,'PaperOrientation','landscape'); end

  ps=get(gcf,'PaperSize');

  mppr=1;
  set(gcf,'PaperPosition', [mppr, mppr, ps(1)-2*mppr, ps(2)-2*mppr])

% -------------------------------------------------------------------- %
  if ~isempty(pos)
      if (dim1>0 || dim2>0) && size(pos,1)~=dim1*dim2, error('Wb:ERR',...
        'dimension mismatch %d != %d*%d\n',size(pos,1),dim1,dim2);
      end

      pos=zoom_axpos(pos,zm);

      n=size(pos,1); ah=zeros(1,n);
      for i=1:n
         ah(i)=axes('Units','normalized','Position',pos(i,:));
      end

      figure(gcf);
      setax(ah(1,1));

      return
  end

  ah=zeros(dim1, dim2);
  if dim1==1 && dim2==1
     pos=zoom_axpos([1.8*m, 1.8*m, 1-2.6*m, 1-3.2*m],zm);
     ah=axes('Position',pos);
     setahs(ah); return
  end

  for i=1:dim1
      for j=1:dim2
          x  = DX   +x0+X0 + (j-1)*(sizex+mx);
          y  = DY+Y1-y0    - (i  )*(sizey+my) + my;
          dx = sizex;
          dy = sizey;

          pos=zoom_axpos([x,y,dx,dy],zm);

          ah(i,j) = axes('Units','norm','Position',pos,'Box','on');

          if (tight && i~=dim1); set(gca,'xtick',[],'ZTickMode','manual'); end
          if (tight && j~=1);    set(gca,'ytick',[],'ZTickMode','manual'); end

          if 0
             set(ah(i,j), 'FontSize', 12)
             set(get(ah(i,j),'xlabel'), 'FontSize', 12)
             set(get(ah(i,j),'ylabel'), 'FontSize', 12)
          end
      end
  end

  for i=1:dim1
  for j=1:dim2
      set(ah(i,j),'Units','points');
      set(ah(i,j),'Position',round(get(ah(i,j),'Position')));
      set(ah(i,j),'Units','normalized');
  end
  end

  setahs(ah);

  addt2fig Wb
  if isempty(ah), return; end

  setax(ah(1));

  if ftag
     l=dbstack;
     l=l(2).name;
     if ftag>1
        w=which(l,'-all');
        if numel(w)>1
           wblog('WRN','got multiple scripts %s (%g)',l,numel(w)); 
        end
        w=w{1};  i=find(w=='/'); if i(1)~=1, w
        error('Wb:ERR','\n   ERR got ususable full script id'); end
        mat=sprintf('%s%s_%s',w(1:i(end)),wbstamp,w(i(end)+1:end));
     else
        mat=sprintf('%s_%s',wbstamp,l)
     end
     setuser(gcf,'mat',mat);
  end

end

% -------------------------------------------------------------------- %
function pos=zoom_axpos(pos,z)

   if isempty(z) || all(z==1), return; end
   if numel(z)==1, z(2)=z; end

   s=size(pos);
   if ~isnumeric(pos) || all(s~=4) || numel(s)>2
      error('Wb:ERR','invalid axes position settings'); end
   if s(1)==4 && s(2)~=4, pos=pos.'; end

   z=repmat(z,size(pos,1),1); pos_=pos;
   pos(:,[1 2])= pos(:,[1 2])+ pos(:,[3 4]).*((1-z)/2);
   pos(:,[3 4])=               pos(:,[3 4]).*(z  );

end

% -------------------------------------------------------------------- %

