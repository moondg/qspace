function [dkt_,b_]=set_global(iflag)
% Function: [dkt,batch]=set_global([opts])
%   
%    typically called within startup.m
%
% Options
%
%    --init     (re)intialize global variables from scratch
%    --test-dkt return determined values for ml_info.{desktop,batch} 
%               without touching ml_info.
%
% Wb,Aug07,07 ; Wb,Aug18,22

% NB! do not use variables
% as these are easily cleared, e.g., using `clear all'

% Wb,Aug18,22: introduced structure ml_env // former field `batch'
% with the additional field `desktop', while also moving mroot into ml_env
% => wblog color codes via escape sequences do not work with desktop
% => apply colors only for ml_envdesktop==1

  if ~nargin, iflag=0;
  elseif isequal(iflag,'--init'), iflag=1;
  elseif isequal(iflag,'--test-dkt')
     [dkt_,b_]=check_ml_mode(); return
  else wbdie('invalid usage'); end

  s=get(0,'UserData');
  if ~isstruct(s) && ~isempty(s)
     wbdie('invalid getuser(0) - empty or struct expected'); 
  end

  if iflag || ~isfield(s,'err_count')
     s.err_count=struct('err',0,'wrn',0,'tst',0);
  end

  if iflag || ~isfield(s,'ml_env')
     [dkt,b]=check_ml_mode();
     s.ml_env=struct('mroot',getenv('MYMATLAB'),'desktop',dkt,'batch',b);
  end

  set(0,'UserData',s);

end

function [dkt,b]=check_ml_mode()
   dkt=0; b=0;
   if usejava('Desktop'), dkt=1;
   elseif isdeployed(), b=1;
   elseif ~isempty(getenv('SGE_O_HOST')) && ...
          ~isempty(getenv('SGE_O_HOME')), b=2;
   elseif ~isempty(getenv('PBS_JOBID')) && ...
          ~isempty(getenv('PBS_O_WORKDIR')), b=3;
   else dkt=2;
   end
end

