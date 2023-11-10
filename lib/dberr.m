% Stop with debugger in case of error
% See also dbwrn()
% Wb,Jun29,11

% dbclear all
  dbstop if error

% setuser(0,'dbstop_if_WRN',1);
  setuser(0,'dbstop_if_ERR',1);

