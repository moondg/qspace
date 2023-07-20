function  [th,ah] = header(varargin)
% function  [th,ah] = header ([pos,][dpos,][opts,], fmt, ... [,{text options}])
% add header/footer text to current figure
%
%   the string is given by {fmt,...} 
%   if a single argument is given, str=fmt, otherwise arguments
%   are run through sprintf
%
%   the last argument, if it is a cell object, is used as
%   set of options to the text object
%
%   the position can be specified explicitly (normalized units)
%   or through the following indicators
%
%      pos (equivalent)       description
%                           
%      -01 'hleft'   'NW'     header left (default)
%       00 'hcenter' 'N'      header center
%      +01 'hright'  'NE'     header right
%      -11 'fleft'   'SW'     footer left
%       10 'fcenter' 'S'      footer center
%      +11 'fright'  'SE'     footer right
%
%   dpos allows a relative shift to specified position (relative units).
%   if single number is specified, this is considered a vertical shift.
%
% Options
%
%   '-append'   append string to existing string, if any, separated by ';  '
%
% Wb,Nov19,01 Nov20,06

  if ~nargin && nargout==1
      th=findall(gcf,'tag','frame');
      return
  end

  if nargin<1, eval(['help ' mfilename]), return, end

% -------------------------------------------------------------------- %
  ah0=gca;

  ah = findall(gcf,'tag','frame');
  if isempty(ah) 
  ah=axes('Units', 'norm', 'Position', [0.002 0.002 0.996 0.996], ...
          'Color', 'none', ...
          'Box', 'on', ...
          'XTick', [], 'YTick', [], ...
          'Tag', 'frame');
  end
  set(ah, 'HitTest', 'off', ...
          'Selected','off', ...
          'Units', 'norm');

  if length(varargin)==1
     if strcmp(varargin{1},'off')
        set(ah,'Visible','off'); set(get(ah,'Children'),'Visible','off');
        return
     elseif strcmp(varargin{1},'on')
        set(ah,'Visible','on');  set(get(ah,'Children'),'Visible','on');
        return
     end
  end

  set(gcf,'CurrentAxes',ah);

% -------------------------------------------------------------------- %
  pos=[]; geth=0;
  if iscell(varargin{1})
     n=numel(varargin{1});
     if n
        tag=lower(varargin{1}{1});
        if n==2, dpos=varargin{1}{2};
        else wbdie('invalid usage (pos)'); end
     end
  else
     tag=lower(varargin{1});
     dpos=[];
  end

  if ischar(tag)
     switch tag
       case { -01, 'hleft',  'nw'}, pos=[ 0, 1]; geth=1; tag='hleft';
       case {  00, 'hcenter','n' }, pos=[.5, 1]; geth=1; tag='hcenter';
       case { +01, 'hright', 'ne'}, pos=[ 1, 1]; geth=1; tag='hright';
       case { -11, 'fleft',  'sw'}, pos=[ 0, 0]; geth=1; tag='fleft';
       case {  10, 'fcenter','s' }, pos=[.5, 0]; geth=1; tag='fcenter';
       case { +11, 'fright', 'se'}, pos=[ 1, 0]; geth=1; tag='fright';
     end
     if nargin==1 && geth
        th=findall(ah,'tag',tag);
        return
     end
  elseif numel(tag)==2 && isnumeric(tag)
     pos=varargin{1}; tag='user';
  end

  if isempty(pos), tag='hleft'; pos=[0 1];
  else varargin(1)=[]; end

  if isempty(dpos)
     if nargin>1 && isnumeric(varargin{1}) && numel(varargin{1})<=2
        dpos=varargin{1}; varargin(1)=[];
        if isscalar(dpos), dpos=[0 dpos]; end
        pos=pos+dpos;
     end
  end

  if numel(varargin) && isequal(varargin{1},'-append')
       add_flag=1; varargin(1)=[];
  else add_flag=0; end

  m=dbstack('-completenames');
  if numel(m)>1, i=regexp(m(2).file,'@');
     if isempty(i), m=m(2).name;
     else
        m=regexprep(m(2).file(i+1:end),'\.m$','');
        m=regexprep(m,'\/.*\/','\/');
     end
  else m='base'; end

  topts={}; if length(varargin)>0
     if iscell(varargin{1}) 
        topts=varargin{1}; varargin(1)=[];
     elseif iscell(varargin{end}) 
        topts=varargin{end}; varargin(end)=[];
     end
  end

  s=varargin{1}; varargin(1)=[]; narg=numel(varargin);

  i=findstr(s,'%MM'); if ~isempty(i)
    q=lineno('all','-s','nx',1); q=strrep(strrep(q,'.m',''),'_','\_');
    q=regexprep(q,' *-> *','::');
    s=strrep(s,'%MM',q);
  end

  s=strrep(s,'%H',hostname);
  s=strrep(s,'%M', strrep(m,'_','\_'));

  s=regexprep(s,'\\+n',char(10));

  if narg
     s=regexprep(s,'\\+','\\\\');
     s=sprintf(s,varargin{:});
  elseif ~isempty(findstr(s,'%'))
     wbdie('invalid remaining format specifiers');
  end
  str=s;

  tiny=[0.012 0.005];

  i=(pos<=0)-(pos>=1);
  pos=pos + i.*tiny;

  if     pos(1)<0.25, halign='left';
  elseif pos(1)>0.75, halign='right'; else halign='center'; end
  if     pos(2)>0.50, valign='top';   else valign='bottom'; end

  th=findall(ah,'tag',tag);

  if isempty(th) || isequal(tag,'user')
     p = get (gca, 'Position');
     pos = ([1/p(3)  0;  0  1/p(4)] * (pos - [p(1) p(2)])')';

     if ~isempty(dpos)
        pos=pos+dpos(1:2);
     end

     th=text(pos(1),pos(2),str,'Tag',tag, ...
       'Units', 'normalized', ...
       'HorizontalAlignment',halign,'VerticalAlignment',valign ...
     );
  else
     if ~isempty(dpos)
        dpos(end+1:3)=0;
        set(th,'Position', dpos+get(th,'Position'));
     end
     if add_flag
        s=get(th,'string');
        if ~isempty(s)
           if size(s,1)>1, str=strvcat(s,str);
           else str=[s ';  ' str]; end
        end
     end
     set(th,'string',str);
  end

  if ~isempty(topts), set(th,topts{:}); end

  mv2back(ah);
  setax(ah0);

  if ~nargout, clear th ah; end

end

