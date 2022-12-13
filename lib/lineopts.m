function opts=lineopts(varargin)
% Functions: opts=lineopts(lh [,flag])
% 
%    return lineopts from given line handle as cell array
%    if lh is an array, then only the line options of the
%    first element are returned.
% 
% Wb,Nov22,07 ; Wb,Aug04,16

  getopt('init',varargin);
     aflag=getopt('-a');
     cflag=getopt('-c');
  lh=getopt('get_last',[]);

  if isempty(lh) || ~ishandle(lh)
     helpthis, if nargin || nargout
     wberr('invalid usage'), end, return
  end
  lh=reshape(lh,1,[]);

  if ~all(isline(lh)), get(lh,'Type')
     wberr('got non-line handles');
  end

  lo={ 'Color','MarkerEdgeColor','MarkerFaceColor',  ...
	   'LineStyle','LineWidth','Marker','MarkerSize' };

  if     cflag, lo=lo(1:3);
  elseif aflag, lo=[lo, 'DisplayName','Tag'];
  end

  lo(2,:)=get(lh(1),lo);
  opts=reshape(lo,1,[]);

end

