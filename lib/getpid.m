function pid=getpid(varargin)
% Function pid=getpid()
% Wb,Jan23,18 

% this returns root process of current matlab session (cf. 'ps -ax')
% e.g., including the command line options used when starting matlab
% same as when using mex-file -> getpid()
% For comparison: 'echo $PID' returns child process
% i.e. referring to the current prompt
  pid=feature('getpid');

end

