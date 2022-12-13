
% script to evaluate input arguments varargin{:}
% expecting L=12 someflag=2 etc.
% Wb,Jul11,16

% outsourced from run_Hamilton1D.m

  if ~isempty(varargin)
     cmd=sprintf('%s ; ',varargin{:});
     fprintf(1,'\n>> %s\n\n',cmd);

     try, eval(cmd); clear cmd
     catch l
         wblog('ERR','%s',l.message);
         dispstack(l);
     end
  end

  jid=str2num(getenv('JOB_ID'));
  tid=str2num(getenv('SGE_TASK_ID'));

