function [jid,tid]=getJobID(varargin)
% function [jid,tid]=getJobID([opts])
% Options
%
%   '-s'  return strings (rather than numbers).
%
% Wb,Jan09,17

  if ~isempty(getenv('PBS_O_HOME'))
     jid=get_job_id('PBS_JOBID','JOB_ID');
     tid=get_job_id('PBS_ARRAYID','PBS_TASKNUM','TASK_ID');
  else
     jid=get_job_id('SLURM_ARRAY_JOB_ID','SLURM_JOB_ID','JOB_ID');
     tid=get_job_id('SLURM_ARRAY_TASK_ID','TASK_ID');
  end

  if ~nargin
     jid=str2num(jid);
     tid=str2num(tid);
  elseif isequal(varargin,{'-s'})
     if nargout<2 && ~isempty(tid)
        jid=sprintf('%s.%02g',jid,str2num(tid));
     end
  else
     helpthis, wberr('invalid usage')
  end
end

% -------------------------------------------------------------------- %

function jid=get_job_id(varargin)
   jid='';
   for i=1:nargin
      jid=getenv(varargin{i}); if ~isempty(jid), break; end
   end
   jid=regexprep(jid,'\.[a-zA-Z].*','');
   jid=regexprep(jid,'\[.*\]','');
end

% -------------------------------------------------------------------- %

