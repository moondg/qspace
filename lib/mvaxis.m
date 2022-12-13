function mvaxis(ah,dxy)
% Function: mvaxis - moves current axis by [dx dy] (units normalized)
% 
% Usage: mvaxis(ah, [dx dy])
% 
%     (axis handle is optional)
% 
% Wb, Aug01,05

  if nargin==0 | nargin>2, eval(['help ' mfilename]), return, end
  if nargin==1, dxy=ah; ah=gca; end

  dxy(end+1:4)=0;

  for h=reshape(ah,1,[])
     set(h,'Units','Normalized');
     set(h,'Position',get(h,'Position')+dxy);
  end

end

