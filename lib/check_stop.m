function q=check_stop()
% function q=check_stop()
%
%    Abuse matlab's capture of output to effectively 
%    permit to interact with running job in an interactive sesseion.
%    This is rather fragile though, since not other system capture
%    command must be issued elsewhere in the meantime!
%
% Wb,Jun25,21

% sleep(5); check_stop
% => then type `stop' while sleeping, this actually works
%    but gives an error, because once returning to the command line
%    prompt, this tries to execute `stop' as new command
% Wb,Jun25,21

% check for `stop' string entered in interactive session! // Wb,Jun25,21
  [e,s]=system('echo'); s=regexprep(s,'\s','');
  if ~isempty(regexpi(s,'stop'))
     if ~nargout
         wberr('got interrupt string ''%s''',s);
     else q=1; end
  else q=0; end

end

