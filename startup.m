
% disp([10 'Run ~/Matlab/startup ...'])

% ------------------------------------------------------------- %
% irrelevant (matlab by default looks in current directory first!)
% path('./',path);
% mr=pwd; 

  mr=getenv('MYMATLAB');

  mq='/Class/@QSpace'; estr='';
  if isempty(mr) && exist([ pwd mq ],'dir'), mr=pwd;
   % accept pwd if suitable, e.g.,
   % in case matlab was called with command line option '-sd ...'
     fprintf(1,'   startup.m: QSpace using (MYMATLAB) = %s\n\n',...
     strrep(mr,getenv('HOME'),'~'));
  elseif isempty(mr)
     estr='matLab path MYMATLAB not set';
  elseif ~exist(mr,'dir')
     estr='invalid matLab path MYMATLAB=%s';
  elseif ~exist([mr mq],'dir')
     estr=sprintf('invalid path MYMATLAB=%s (QSpace not found)',mr);
  elseif ~exist([mr '/DMRG'],'dir')
     estr=sprintf('invalid path MYMATLAB=%s [(MYMATLAB)/DMRG not found]',mr);
  end

  if estr, mr='';
   % error('QSpace:ERR:startup', ...)
     fprintf(2, [char(27) '[31m' ... 
     '   ERR startup.m %s\n' ...
     '   ERR startup.m PATH not set\n' char(27) '[0m\n\n' ],estr);
  end

% adapte for working package running 'ml -mp2' 
% with path mp2 set to *this directory within cto.dat
% Wb,Oct15,14

if mr
% NB! mex-files before other files to ensure correct
% behavior of `help some-mex-file' // Wb,Jan13,19
  path(path,[mr '/bin'   ]);
  path(path,[mr '/util'  ]);

  path(path,[mr          ]);
  path(path,[mr '/lib'   ]);
  path(path,[mr '/NRG'   ]);
  path(path,[mr '/DMRG'  ]);
  path(path,[mr '/tensor']);
  path(path,[mr '/setup' ]);
  path(path,[mr '/Class' ]);

  set_global % global count (wrn,err,tst); requires (mr)/lib
end

% ------------------------------------------------------------- %

  gr=groot;
  ss=get(gr,'ScreenSize'); ss=ss(3:4);
% fp=get(gr, 'DefaultFigurePos'); fp=fp(3:4);
  fp=[590 520]; % fp=[420 400];
  fp=[ ss-fp-[4 74], fp ] ;

% starting MatLab without display sets ScreenSize = [1 1 1 1] (!)
  if ~isbatch && all(ss>1)
  set(gr,'DefaultFigurePosition',fp); end

  set(gr,'DefaultFigureName', getenv('HOST'))
  set(gr,'DefaultFigurePaperType','A4');

% works better for mjpg to avoid coarse resolution for written text
% while the rest of the figure seems good; Utopia; default: Helvetica
  set(gr,'DefaultTextFontName', 'Arial');
  set(gr,'DefaultAxesFontName', 'Arial');

  set(gr,'DefaultTextFontSize',  12); % default: 10
  set(gr,'DefaultAxesFontSize',  12); % default: 10
  set(gr,'DefaultAxesLineWidth', 1.); % default: .5
  set(gr,'DefaultLineLineWidth', 1.); % default: .5

% rand('state',sum(100*clock)); // deprecated // Wb,Mar13,20
  rng shuffle

  clear mr mq estr gr ss fp

% startup_aux(which(mfilename), 'opengl neverselect (see readme.txt)');
% opengl neverselect

  if exist('startup_loc.m')==2, startup_loc; end

