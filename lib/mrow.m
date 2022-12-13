function dd = mrow(d0)
% Function: mrow
% Usage   : dd = mrow(d0)
%
%    a non-square d0 is transposed into a row like object 
%    eg. size(d0,2) >= size(d0,1)
%
% Wb,Jun18,03

  dd = d0;

  if size(dd,2)<size(dd,1), dd = dd.'; end

return

