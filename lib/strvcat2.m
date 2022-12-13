function s=strvcat2(varargin)
% Function: s=strvcat2(varargin)
%
%    same as strvcat, except that it aligns text
%    to the right (leading blanks instead of trailing ones)
%
% Wb,May02,07

  s=char(varargin); % that's what strvcat does anyway
  [m,n]=size(s);

  for i=1:m
     for j=n:-1:1, if s(i,j)~=' ', j=j+1; break, end; end
     s(i,:)=[s(i,j:end), s(i,1:j-1)];
  end

end

