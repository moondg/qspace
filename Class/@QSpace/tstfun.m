function tstfun(A)
% function tstfun(A)
% Wb,Feb23,11

  wblog('1'); whos A

% BEHAVIOUR WITHIN CLASS METHODS
% whenever a class method requires the functionality of the overloaded
% subsref or subsassign, it must call the overloaded methods with
% FUNCTION CALLS rather than using operators like '()', '{}' or '.'

  A(1)
  subsref(A, struct('type', '()', 'subs', {{1}}) )

  wblog('2');

  helpthis
  tstfun_aux(A)

end

function tstfun_aux(A)

  s=dbstack; s=s(1);
  helpthis

end

