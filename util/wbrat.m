%  Usage: [q,I] = wbrat(q [,opts]);
%
%     get rational approximation of given number
%     based on continued fraction approach.
%
%  Options
%
%     '-h'      display this usage
%     '-v'      verbose mode
%     '-q'      quiet mode
%     '~r'      do not include square root check [i.e., using rat() on q^2]
%     '~s'      do not return string, but rath numbers q=P/Q
%     '-a'      also return continued fraction sequence in I
%     'nmax',.. maximimum number of iterations (6)
%     'eps1',.. tolerance on continued fraction (1E-6)
%     'eps2',.. rel. overall tolerance wrt. to in put number q (1E-12)
%
%  Output
%
%      q        same as input, but 'fixed' to rational approximation
%      I        info structure
%               .P  integer enumerator
%               .Q  integer denominator
%               .relerr relative error of rational representation
%               .n  total number of continued fraction iterations used
%               .param  (input) parameters nmax, eps1, eps2 used
%
%  See also MatLab's rat.m
%  AWb (C) July 2010
