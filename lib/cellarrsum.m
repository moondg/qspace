function A=cellarrsum(C)
% Function: cellarrsum
% Usage: A=cellarrsum(C);
%
%    A = sum_i C{i}
%
% Wb,Nov30,05

  A=C{1}; n=numel(C);
  for i=2:n, A=A+C{i}; end

end

