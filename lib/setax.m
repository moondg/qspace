function setax(varargin)
% Function: setax ([ax or fig handle, OPTS])
%
%    set line widths and font sizes of current figure
%
% Options
%
%   '-gcf'   apply to all axes sets in current figure
%   'lw',... line width for all 'line' objects
%   'lwf',.. (same as lw, but scale factor)
%   'aw',... line width for all 'axes' objects
%   'awf',.. (same as aw, but scale factor)
%   'fsa',.. font size for all 'axes' objects
%   'fst',.. font size for all 'text' objects
%   'lbl',.. font size for xyz labels and title
%
% Wb,Dec11,03  Wb,Sep10,07
% ------------------------------------------------------------------

  if nargin==1, a=varargin{1};
     if ischar(a)
        ah=findall(groot,'type','axes','tag',a);
        if isempty(ah)
           wbdie('axes tag ''%s'' not found');
        elseif numel(ah)>1
           for i=1:numel(ah)
              s=get(get(ah(i),'parent'),'Name'); %% 'untagged'
              if s(1)=='[', ah(i)=-ah(i); end
           end
           i=find(ah>0); if numel(i)==1, ah=ah(i);
           else
              wblog('WRN','found %g axes sets with tag ''%s'' (%g)',...
              numel(ah),a,numel(i)); ah=abs(ah(1));
           end
        end
     elseif isaxis(a)
        ah=a;
     elseif isnumeric(a) && isequal(round(a),a)
        ah=getuser(gcf,'ah'); n=numel(ah);
        if ~n || a>n
           wbdie('invalid setuser axes set or index (%g/%g) !?',a,n); end
        ah=ah(a);
     else wbdie('invalid usage'); end

     n=numel(ah);
     for i=1:n
        f=get(ah(i),'Parent');
        set(f,'CurrentAxes',ah(i)); set(0,'CurrentFigure',f);
     end
     return
  end

  if nargin<2, helpthis
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  a=varargin{1};
  if nargin && isnumeric(a) && (isaxis(a) || isfig(a))
       par=a; varargin(1)=[];
  else par=[]; end

  getopt('init',varargin);
     fflag=getopt('-gcf'  );
     lw   =getopt('lw', []);
     lwf  =getopt('lwf',[]);
     aw   =getopt('aw', []);
     awf  =getopt('awf',[]);
     fsa  =getopt('fsa',[]);
     fst  =getopt('fst',[]);
     fs   =getopt('fs', []); if ~isempty(fs), fsa=fs; fst=fs; end
     lbl  =getopt('lbl',[]);
  getopt('check_error');

  if ~isempty(lw) && ~isempty(lwf), eval(['help ' mfilename])
     wblog('ERR - <line> and <linef> are exclusive parameters!');
     return
  end

  if ~isempty(aw) && ~isempty(awf), eval(['help ' mfilename])
     wblog('ERR - <axes> and <axesf> are exclusive parameters!');
     return
  end

% ------------------------------------------------------------------

  if isempty(par)
     if fflag, par=gcf; else par=gca; end
  end

  ah=findall(par,'type','axes');
  lh=findall(par,'type','line');
  th=findall(par,'type','text');

  if ~isempty(lw ); set(lh,'LineWidth', lw ); end
  if ~isempty(aw ); set(ah,'LineWidth', aw ); end
  if ~isempty(fsa); set(ah,'FontSize',  fsa); end
  if ~isempty(fst); set(th,'FontSize',  fst); end

  if ~isempty(lbl)
     for i=1:numel(ah)
        set(get(ah(i),'XLabel'),'FontSize',lbl);
        set(get(ah(i),'YLabel'),'FontSize',lbl);
        set(get(ah(i),'ZLabel'),'FontSize',lbl);
        set(get(ah(i),'Title' ),'FontSize',lbl);
     end
  end

  if ~isempty(lwf)
      for h=lh; set(h, 'LineWidth', fact * get(h,'LineWidth')); end
  end

  if ~isempty(awf)
      for h=ah; set(h, 'LineWidth', fact * get(h,'LineWidth')); end
  end

end

