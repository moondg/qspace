function [p,f,x]=fileparts2(f)
% function [p,f,x]=fileparts2(fname)
%
%    sligthly adapted version of MatLabs fileparts:
%
%    path name always contains trailing /
%    (if no path is specified, ./ is returned)
%
%    extension starts from the very first period.
%
% Wb,May21,09

  i=find(f=='/');
  if isempty(i), p='./'; else
     p=f(1:i(end));
     f=f(i(end)+1:end);
  end

  i=find(f=='.');
  if isempty(i), x=''; else
     x=f(i(1):end);
     f=f(1:i(1)-1);
  end

end

