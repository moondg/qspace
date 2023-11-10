function h=poslabel(varargin)
% Function h=poslabel([ah,] n [,pos,opts])
%
%    position figure / panel label (a), ...
%
% Options
%
%    'istr',. append info string to label
%    'dx',..  (in normalized units in current axis set)
%    'fs',..  font size of label (default: 16)
%    '-wb'    white, i.e., erase  background
%    '-wt'    white text, e.g., on top of dark surface plot
%    '-rm'    remove existing labels
%    '-fix'   adapt existing label (e.g. in terms of options above)
%
%     remaining options are applied to text object.
%
% Wb,Apr09,11

% example usage
% poslabel(ah(1),1,{'NW',[-0.3 0]})

  ah0=gca; 
  if nargin && isaxis(varargin{1})
       ah=reshape(varargin{1},1,[]); setax(ah(1)); varargin(1)=[];
  else ah=ah0; end

  if numel(varargin)<1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  getopt('init',varargin);
     rmflag=getopt('-rm');
     fix   =getopt('-fix'); o={};

     if getopt('-wb'), o={'BackgroundC','w','Margin',0.1};
     elseif getopt('-wt'), o={'Color','w'};
     else o={}; end

     fs    =getopt('fs',16);
     dx    =getopt('dx',[]);
     istr  =getopt('istr','');
  args=getopt('get_remaining');
  nargs=numel(args); pos={}; 

  tag=mfilename; % 'poslabel'

  h=findall(ah,'type','text','tag',tag);
  if ~fix, delete(h); end

  if nargs && isempty(dx), q=args{1};
     if iscell(q) && numel(q)==2 && ischar(q{1})
        pos=q; args(1)=[]; nargs=numel(args);
     end
  end

  if ~nargs
     if fix
        if numel(h)
           if numel(h)>1, wbdie('got multiple labels !?'); end
           n=get(h,'string');
        else
           wblog('WRN','%s: no label to fix',mfilename);
           clear h, return
        end
     else
        if rmflag && ~nargout, clear h; return, end
        wbdie('invalid usage');
     end
  else
     n=args{1}; args(1)=[]; nargs=numel(args);
  end

  if nargs && isempty(pos), q=args{1};
     if ischar(q) && numel(q) && numel(q)<=2 && isequal(q,upper(q))
        pos=q; args(1)=[];
     elseif iscell(q) && numel(q)==2 && ischar(q{1}) && isempty(dx)
        pos=q; args(1)=[];
     end
  end

  if  isempty(pos), pos='NW'; else pos_=pos; end
  pos=postrans(pos);
  if ~isempty(dx)
     if numel(dx)~=2, wbdie('invalid dx'); end
     pos=pos+reshape(dx,1,[]);
  end

  if ischar(n)
     if ~isempty(n) && (n(1)=='(' || n(1)=='\' || n(1)==' ')
     s=n; else s=sprintf('(%s)',n); end
  else
     if ~isnumber(n), n, wbdie('invalid usage'); end
     s=sprintf('(%s)',char('a'+n-1));
  end
  h=postext(pos,s);

  if ~isempty(istr)
     set(h,'String',[get(h,'String') ' ' istr]);
  end

  if isempty(dx)
     if pos(2)>0.8, o={o{:},'VerticalAl','top'};
     elseif pos(2)<0.2, o={o{:},'VerticalAl','bottom'}; end
  end

  set(h,'tag',tag,o{:},'FontSize',fs,args{:});
  if ~nargout, clear h; end

  if ~isequal(ah,ah0), setax(ah0); end

end

