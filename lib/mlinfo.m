function I=mlinfo(varargin)
% function I=mlinfo()
% See also startup_aux.m
% Wb,Sep21,15

  getopt('init',varargin);
     dflag=getopt('--start');
  getopt('check_error');

  if dflag
     I.started=datestr(now);
     I.finished='';
  else
     I.date=datestr(now);
  end

  I.host={ hostname, hostid };
  I.pid=getpid();
  I.pwd=pwd;
  I.cputime=cputime;

  I.display=(usejava('awt') && usejava('jvm'));
  I.status=struct( ...
     'ismcc',     ismcc, ...
     'isbatch',   isbatch, ...
     'isdeployed',isdeployed, ...
     'isjava',  [ usejava('awt'), usejava('jvm') ],...
     'desktop',   usejava('desktop'), ...
     'debug',   [ get_env('DEBUG'), get_env('ML_DEBUG') ], ... % Wb,Feb06,17
     'Display', {{getenv('DISPLAY'), get(groot,'DefaultFigureRenderer')}} ...
  );

% usejava('jvm')    %  1 even if matlab -nodisp, 0 if matlab -noj
% usejava('awt')    %  0 if matlab -nodisp
% usejava('swing')  %  0 if matlab -nodisp
% usejava('desktop') : 0 if matlab is run in command-line mode(!)

  [i,~]=system('command -v omp_threads.pl'); if ~i % user routine
  [i,s]=system('omp_threads.pl -x'); s=str2num(s); q=s(end);
  else q='(?)'; end

  I.nthreads=struct(...
    'cores',   q, ...
    'slots',   get_env('NSLOTS'), ...
    'OMP',     get_env('OMP_NUM_THREADS'), ...
    'MKL',     get_env('MKL_NUM_THREADS'), ...
    'QSP',     get_env('QSP_NUM_THREADS')  ...
  );

% NB! maxNumCompThreads is deprecated with matlab >= R2009
% n=str2num(n); i=maxNumCompThreads(n);
% setuser(groot,'maxNumCompThreads',i);
% Wb,Feb10,11
  I.cgs=struct(... % CGS info
    'store',      getenv('RC_STORE'),   ...
    'nthreads',   I.nthreads.QSP,       ...
    'verbose',    get_env('CG_VERBOSE') ...
  );
  % 'num_threads',str2num(getenv('CG_NUM_THREADS')), // Wb,Feb06,17
  % 'qs_num_threads',get_env('QS_NUM_THREADS'), ...
  % 'sp_num_threads',get_env('SP_NUM_THREADS'), ...

% NB! matlab understands strings 'true/false' as 1/0 values
% => correctly treated in if-statements or also str2num()!
% Wb,Feb12,16
  I.omp=struct(...
    'nthreads',   I.nthreads.OMP, ...
    'max_active_levels', get_env('OMP_MAX_ACTIVE_LEVELS'), ...
    'dynamic',    getenv('OMP_DYNAMIC'), ...
    'nested',     getenv('OMP_NESTED'),...
    'schedule',   getenv('OMP_SCHEDULE'), ...
    'stack_size', get_env('OMP_THREAD_STACK_SIZE') ...
  );

% https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications
% Wb,Feb12,16
  I.mkl=struct(...
    'nthreads',   I.nthreads.MKL, ...
    'nthreadD',   get_env('MKL_DOMAIN_NUM_THREADS'), ...
    'dynamic',    getenv('MKL_DYNAMIC'), ...
    'nested',     getenv('MKL_NESTED'), ...
    'threading_layer',getenv('MKL_THREADING_LAYER') ...
  );

  [jid,tid]=getJobID('-s'); % Wb,Jan09,17
  if isempty(tid)
       if ~isempty(jid), I.jobid=jid; else I.jobid={}; end
  else I.jobid={jid,tid};
  end

end

% -------------------------------------------------------------------- %
% Wb,Feb06,17

function q=get_env(varargin)
   q=str2num(getenv(varargin{1}));
end

% -------------------------------------------------------------------- %
% function q=num_threads(varargin)
%    n=nargin; if n==1, q=str2num(getenv(varargin{1})); return; end
%    
%    q=cell(1,n); mark=zeros(1,n);
%    for i=1:n
%       q{i}=getenv(varargin{i});
%       if ~isempty(q{i}), mark(i)=1; q{i}=str2num(q{i}); end
%    end
% 
%    if any(mark)
%         q(find(~mark))={1}; q=[q{:}];
%    else q=[]; end
% end
% -------------------------------------------------------------------- %

