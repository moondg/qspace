function cd_start(p)
% function cd_start(p)
% auxilliary routine for matlab startup
% see ml, mbatch, etc.; former startup_aux.m
% Wb,Jun08,19

  builtin('cd',p);
  fprintf(1,'>> cd %s\n\n',repHome(pwd));

end

