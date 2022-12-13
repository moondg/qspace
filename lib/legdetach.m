function [ah,hh]=legdetach(l)
% Function [ah,hh]=legdetach(l)
%
%    detach legend from current axis
%    by creating a copy of the legend axis set.
%
% Wb,Apr17,08 ; Wb,Aug04,16

% by matlab>2013b legend is no longer of axis type,
% but is an object type of its own.
% for old version, see Archive/legdetach_160801.m

  h0=gca;
  ah=axes('Position',get(l,'Position'));

  if ~isnumeric(l)
     wberr('legdetach() depricated (merged with legdisp)'); 
  end

  h1=copyobj(findall(l,'type','text'),ah); set(h1,'tag','','UserD',[]);
  h2=copyobj(findall(l,'type','line'),ah); set(h2,'tag','','UserD',[]);
  hh=[h1;h2];

  set(ah, ...
    'XLim',get(l,'XLim'), 'XColor',get(l,'XColor'), ...
    'YLim',get(l,'YLim'), 'YColor',get(l,'YColor'), ...
    'Color',get(l,'Color'), 'visible',get(l,'visible'),...
    'tag','legdetached' ...
  );
  axis(ah,'off');
  set(gcf,'CurrentAxes',h0); legend off

end

