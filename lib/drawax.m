function [ha,z]=drawax(varargin)
% function h=drawax([ah,lineopts])
%
%    draw box across axis boundary to have axis on top
%    of data drawn (default: data appears on top of axis line!
%
% Options
%
%    'z',..   z-value (e.g. to move to foreground; [])
%    '--top'  draw top line on axes only
%    '--topr' draw top+right line on axes only
%
% Wb,Apr09,09

% rather use: set(gca,'Layer','top')  // Wb,Nov01,20

  if nargin && isaxis(varargin{1})
       ah=varargin{1}; varargin(1)=[];
  else ah=gca; end

  getopt('init',varargin);
     z=getopt('z',[]);
     if     getopt('--top'), ix=[1 2]; iy=[2 2];
     elseif getopt('--topr'),  ix=[1 2 2]; iy=[2 2 1];
     else ix=[1 1 2 2 1]; iy=[1 2 2 1 1];
     end
     qflag=getopt('-q');
  args=getopt('get_remaining');

  if ~qflag, dispstack(dbstack);
     wblog('WRN','consider set(gca,''Layer'',''top'')'); 
  end

  tag=mfilename;
  delete(findall(ah,'tag',tag));

  args={'tag',tag,args{:}};

  n=numel(ah); h=gobjects(1,n); r=[0.0001 0.9999];

  if numel(z)==1 && z
     z=repmat(z,1,5);
  end

  ish=zeros(1,n);

  for i=1:n, hi=ah(i); setax(hi);
     ish(i)=~isempty(regexpi(get(hi,'NextPlot'),'add'));
     if ~ish(i), hold on; end

     logx=isequal(get(hi,'XScale'),'log');
     x=xlim;           if logx, x=log(x); end
     x=x(1)+diff(x)*r; if logx, x=exp(x); end
     x=x(ix);

     logy=isequal(get(hi,'YScale'),'log');
     y=ylim;           if logy, y=log(y); end
     y=y(1)+diff(y)*r; if logy, y=exp(y); end
     y=y(iy);

     ha(i)=plot (x, y,'k','LineW',get(hi,'LineW'),args{:});

     if ~ish(i), hold off; end

     if ~isempty(z)
        if any(z)
           set(ha(i),'ZData',z);
        end
     else
        zmax=0;
        for h=get(gca,'Children')'
            if ~isequal(get(h,'Type'),'text')
               q=get(h,'ZData'); if ~isempty(q)
                  zmax=max(zmax,max(reshape(q,[],1)));
               end
            end
        end
        set(ha(i),'ZData',repmat(2*zmax,size(get(ha(i),'XData'))));
        za(i)=2*zmax;
     end
  end

  if nargout
     if isempty(z), z=zmax; end
  else clear ha, end

end

