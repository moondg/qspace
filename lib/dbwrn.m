
% Wb,Jun29,11
% dbclear all
  dbstop if error
  dbstop if warning

  setuser(0,'dbstop_if_WRN',1);
  setuser(0,'dbstop_if_ERR',1);

