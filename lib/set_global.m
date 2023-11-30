function [dkt,dbg,isb]=set_global(iflag)
% Function: [dkt,debug,batch]=set_global([opts])
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

% Wb,Nov22,23: added debug flag

% NB! do not use variables
% as these are easily cleared, e.g., using `clear all'

% Wb,Aug18,22: introduced structure ml_env // former field `batch'
% with the additional field `desktop', while also moving mroot into ml_env
% => wblog color codes via escape sequences do not work with desktop
% => apply colors only for ml_envdesktop==1

  if ~nargin, iflag=0;
  elseif isequal(iflag,'--init'), iflag=1;
  elseif isequal(iflag,'--test-dkt')
     [dkt,dbg,isb]=check_ml_mode(); return
  else wbdie('invalid usage'); end

  s=get(0,'UserData');
  if ~isstruct(s) && ~isempty(s)
     wbdie('invalid getuser(0) - empty or struct expected'); 
  end

  if iflag || ~isfield(s,'err_count')
     s.err_count=struct('err',0,'wrn',0,'tst',0);
  end

  if iflag || ~isfield(s,'ml_env')
     q={'mroot','desktop','debug','batch'}; q{2,1}=getenv('MYMATLAB');
     [q{2,2:4}]=check_ml_mode();
     s.ml_env=struct(q{:});
  end
  set(0,'UserData',s);

end

% -------------------------------------------------------------------- %
function [dkt,dbg,isb]=check_ml_mode()
   dkt=0; dbg=0; isb=0;
   if usejava('Desktop'), dkt=1;
   elseif isdeployed(), isb=1;
   elseif batchStartupOptionUsed(), isb=1;
   elseif ~isempty(getenv('SGE_O_HOST')) && ...
          ~isempty(getenv('SGE_O_HOME')), isb=2;
   elseif ~isempty(getenv('PBS_JOBID')) && ...
          ~isempty(getenv('PBS_O_WORKDIR')), isb=3;
   else dkt=2;
   end

   if dkt
      if str2num(getenv('ML_DEBUG')), dbg=1;
      elseif ~isempty(getenv('MATLAB_DEBUG')), dbg=2;
      elseif ~isempty(getenv('DEBUG')), dbg=3;
      end
   end
end

% -------------------------------------------------------------------- %

