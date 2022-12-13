function ah = axonax (h0, pos, flag)
% ah = axonax (h0, pos, flag)
%    ah   - new axis handle
%    h0   - existing axis handle
%    pos  - 'R'(ight), 'T'(op), else (RT)
%    flag - 1: inherit xy limits; 2: inherit y-labels (3: and shift them)
%
% creates new axis set at position of given axis, but transparent;
% this is useful for creating double left/right tick sets!
%
% NB! See also MatLab routine plotyy()
%
% Wb,Feb05,02

  if nargin==1, pos='R'; flag=0;
  elseif nargin~=3
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  pos=upper(pos);

  o={ 'Color','none', 'XAxisLocation', 'top', 'YAxisLocation', 'right' };
  if flag>0, o={o{:},'XLim',get(h0,'XLim'), 'YLim',get(h0,'YLim')}; end

  ah=axes ('Position', get(h0,'Position'),o{:},'tag',mfilename);
  setuser(h0,'ih',ah);

  if flag==2
     ytick = get(h0,'YTick');
     yy    = abs(ytick - mean(ytick));
     yy    = find (yy == min(yy));
     set(ah, 'YLim', get(ah,'YLim') - ytick(yy(1)));
  elseif flag==3
     ytick = get(h0,'YTick');
     set(ah, 'YLim', get(ah,'YLim') - ytick(1));
  end

  if any(pos=='R'), set(ah,'XTickLabel',[]); end
  if any(pos=='T'), set(ah,'YTickLabel',[]); end

  hold on

end

