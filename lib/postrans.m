function [pos,opts]=postrans(varargin)
% Function [pos,opts]=postrans([opts,]pos)
%
%   Translates position such as 'NorthEast', ... (similar to legend)
%
% Output
%
%   pos   position in normalized units
%   opts    text alignment (e.g. as used in postext)
%
% Wb,Aug15,08 ; Wb,Jul27,17

  getopt('init',varargin);
     lflag=getopt('--leg');
  pos=getopt('get_last',[]);

  if nargin<1 || isempty(pos), eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if isnumeric(pos)
     if lflag || numel(pos)~=2
        wbdie('invalid position specs');
     end

     if nargout>1
        x=pos(1); y=pos(2);
        if x<0.25, ha='left';   elseif x>0.75, ha='right'; else ha='center'; end
        if y<0.25, va='bottom'; elseif y>0.75, va='top';   else va='middle'; end
        opts={'HorizontalAl',ha,'VerticalAl',va};
     end

     return
  end

  if iscell(pos), dx=pos{2}; pos=pos{1};
  else dx=[]; end

  if ~ischar(pos)
     wbdie('invalid position specs');
  end

  switch lower(pos)
     case {'north',    'n' }, pos='North';     i=[2,3];
     case {'northeast','ne'}, pos='NorthEast'; i=[3,3];
     case {'east',      'e'}, pos='East';      i=[3,2];
     case {'southeast','se'}, pos='SouthEast'; i=[3,1];
     case {'south',    's' }, pos='South';     i=[2,1];
     case {'southwest','sw'}, pos='SouthWest'; i=[1,1];
     case {'west',      'w'}, pos='West';      i=[1,2];
     case {'northwest','nw'}, pos='NorthWest'; i=[1,3];
     case {'center',   'cc'}, pos='Center';    i=[2,2];
     otherwise wbdie('invalid position `%s''',pos);
  end

  if nargout>1
     if i(1)<2, ha='left'; elseif i(1)>2, ha='right'; else ha='center'; end
     if i(2)<2, va='bottom'; elseif i(2)>2, va='top'; else va='middle'; end
     opts={'HorizontalAl',ha,'VerticalAl',va};
  end

if lflag, return; end

  X=[0.05, 0.5, 0.95];
  Y=[0.05, 0.5, 0.95];

  pos=[X(i(1)),Y(i(2))];

if isempty(dx), return; else n=numel(dx); end

  if n==2, pos=pos+dx;
  elseif n==1
     if i(1)==1 || i(1)==3
          pos(2)=pos(2)+dx;
     else pos(1)=pos(1)+dx; end
  else
     wbdie('invalid dx to position (%g)',n);
  end

end

% -------------------------------------------------------------------- %

