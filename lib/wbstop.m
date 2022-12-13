
% NB! do not make *this a function, so it can see callers environment as is
  my_stack__=dbstack; dispstack(my_stack__(2:end)); clear my_stack__

  if isbatch
       wblog('WRN','ignoring wbstop / dbstop in batch mode)');
  else keyboard; end

