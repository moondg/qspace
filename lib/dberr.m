function dberr(varargin)
% Wb,Jun29,11

% cmd='dbclear all ; dbstop if error';
% if nargin, cmd=[cmd, ' ; dbstop ' sprintf(' %s',varargin{:})]; end
% evalin('caller',cmd);

% dbclear all
  dbstop if error

% setuser(0,'dbstop_if_WRN',1);
  setuser(0,'dbstop_if_ERR',1);

end

