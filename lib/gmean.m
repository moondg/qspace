function dm=gmean(dd,varargin)
% function dm=gmean(data,...)
%
%    Same usage as mean (all arguments are forwarded)
%    except that geometric mean is taken.
%
% Wb,Nov28,24

  i=find(dd<=0);
  if ~isempty(i), wblog('WRN',...
     'data <= 0 encountered (%d/%d)',numel(i),numel(dd)); end

  dm=exp(mean(log(dd),varargin{:}));

end

