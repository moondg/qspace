function o=getlopts(k,co)
% function o=getlopts(k [,co])
%
%    get (unique) set of line options
%
% Wb,May30,11

  if nargin<2, co=get(gca,'ColorOrder'); end
  n=size(co,1); j=floor((k-1)/n);

  c=getcolor(k,'cm',co);

  switch mod(j,3)
    case 0, o={'Color',c};
    case 1, o={'Color',c,'LineSt','--'};
    case 2, o={'Color',1-0.3*(1-c)}; % ,'LineW',2
  end

end

