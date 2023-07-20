function map_xdata(varargin)
% function map_xdata([ah,][type],args,...])
% Examples
%
%    map_xdata(gca,['lin',] xrange_from, xrange_to);
%    map_xdata(gca,xd); % assuming current x-data is simple index into xd
%
% Options
%
%   '--log'   switch to log-x scale for (e.g., because xd is log-spaced)
%
% Wb,Sep13,15

% examples
%   map_xdata('lin',[1 L],[0 J2])

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     lflag=getopt('--log');
  args=getopt('get_remaining');

  if isaxis(args{1})
       ah=args{1}; args(1)=[];
  else ah=gca; end

  lin_flag=0;
  if isequal(args{1},'lin') % make 'lin' optional // Wb,May26,17
     args(1)=[];
     lin_flag=1;
  end
  nargs=numel(args); 

  if nargs>=2 && ...
     isnumeric(args{1}) && numel(args{1})==2 && ...
     isnumeric(args{2}) && numel(args{2})==2, lin_flag=lin_flag+10;
  elseif lin_flag
     wbdie('invalid usage (''lin'')');
  end

  if lin_flag
     x0=args{1}; x1=args{2};
     FCT=@(x) (x-x0(1))*(diff(x1)/diff(x0)) + x1(1);
     args(1:2)=[];
  else
     if ~isnumeric(args{1}) || nargs>1 && ~ischar(args{2})
        args, wbdie('invalid usage (map)');
     end
     xd=args{1}; nx=numel(xd);
     args(1)=[];
  end

  getopt('init',args);
     xlb=getopt('xlb','');

  if isempty(xlb)
  xlb=getopt('get_last',''); else
  getopt('check_error'); end

  th=intersect(get(gca,'Children'),findall(ah,'type','text','Units','data'))';
  if ~isempty(th)
     if lin_flag
        for h=th
           p=get(h,'Position'); p(1)=FCT(p(1));
           set(h,'Position',p);
        end
     else
        set(th,'Units','normalized');
     end
  end

  na=numel(ah);

  lh=findall(ah,'type','line')';

  ph=findall(ah,'type','patch')';
  if ~isempty(ph)
     if lin_flag, lh=[lh, ph];
     else wbdie('invalid usage (got %g patches)',numel(ph)); end
  end

  if lin_flag
     xl=cell(1,na); for i=1:na, xl{i}=xlim(ah(i)); end

     for h=lh, set(h,'XData',FCT(get(h,'XData'))); end
     for i=1:na, set(ah(i),'XLim',FCT(xl{i})); end
  else
     for h=lh
        i=get(h,'XData');
        if any(i~=round(i)), wbdie(...
           'invalid existing xdata (expecting integers)');
        elseif any(i<1) || any(i>nx), wbdie(...
           'invalid existing xdata (index out of bounds: [%g,%g]/%g)',...
           min(i),max(i),nx);
        end
        set(h,'XData',xd(i));
     end
     for i=1:na
        setax(ah(i)); xtight
     end
  end

  if lflag, set(gca,'XScale','log'); end

  if ~isempty(xlb)
     for i=1:na
        h=get(ah(i),'XLabel'); s=get(h,'string');
        if ~isempty(s)
             set(h,'string',[s ' \rightarrow ' xlb]);
        else set(h,'string',xlb);
        end
     end
 end

end

