function save2(varargin)
% Function: save2(fname [,OPTS, usual arguments for matlab save])
%
%    Save mat file but with a few more features.
%
% Options
%
%  '-d',.. directory to save to
%  '-q'    quiet
%  '-f'    overwrite file without asking if it exists
%  '-k'    keep time stamp of existing file (if any)
%  '-x',.. list of variables to exclude
%          NB! regular expressions are supported
%          NB! *tol* would not work! use '.*tol.*', instead
%
% Examples
%
%    save2('tmp.mat','-struct','Istruct')
%
% All remaining arguments are handed over to MatLab's save()
% Wb,Aug11,06 / Wb,May20,07 / Wb,Jan02,17

  getopt('INIT',varargin);
     dstr =getopt('-d','');
     excl =getopt('-x',{});

         if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end

     force=getopt('-f');
     kflag=getopt('-k');

  varargin=getopt('get_remaining');

  if length(varargin)<1
     eval(['help ' mfilename]);
     if nargin, wberr('invalid usage'); else return; end
  end

  fname=varargin{1}; varargin(1)=[];

  if ~isempty(dstr)
     ldir=cd;
     if exist(dstr,'dir'), cd(dstr);
     else cto(dstr,'nolog'); end
  else ldir=''; end

% NB! in order to avoid data loss of existing files while terminating
% saving the new data here (e.g. kill of process, reboot, etc.)
% ALWAYS  write file to hidden temporary file .fname~.mat first!
% then move the file to its final destination, given that a system
% /bin/mv operation is much faster as compared to the matlab write
% for large files // Wb,Jan02,17

  [p,f,x]=fileparts(fname); sx=[numel(f), numel(x)];
     if isempty(p), p='.'; end
     if ~sx(1), wblog('WRN','got filename ''%s'' !?',fname); end
     if isempty(regexpi(x,'\.mat$')), f=[f x]; x='.mat';
        if ~sx(1), wblog('WRN','using f=''%s'' -> ''%s%s''',f,f,x); end
     end
  ftmp=[p '/.' f '-tmp' x];
  f=[p '/' f x];

  xflag=exist(f,'file');
  if xflag
     if ~force && fexist(f), return; end
  else
     ftmp=f;
  end

  if vflag, s=repHome(basename(f));
     if length(s)<40
          wblog(1,'I/O','saving data to %s',s);
     else wblog(1,'I/O','saving data to file\n%s',s);
     end
  end

  m=0; nx=0; s=''; sx='';

  if ~isempty(excl)
     evalin('caller','setuser(groot,''whos'',whos);');
     vars=getuser(groot,'whos','-rm'); ix=[];

     if ~iscell(excl), excl={excl}; end

     for i=1:length(vars)
     for q=1:length(excl)
        if ~isempty(regexp(vars(i).name,excl{q}))
           ix=[ix,i]; break;
        end
     end
     end

     if ~isempty(ix)
        nx=length(ix); n=length(vars);
        sx=sprintf(' %s',vars(ix).name);
        if vflag>1
           if length(sx)>50, sx=[char(10) '  ' sx]; end
           fprintf('\n   skipping %d/%d variable(s): %s\n\n',nx,n,sx);
        end

        varx=vars(ix); vars(ix)=[]; m=length(vars);
        s=[' ', sprintf('%s ', vars.name)];
     end
  end

  sv=sprintf(' %s', varargin{:});
  cmd=sprintf('save %s%s%s',ftmp,s,sv);

  if vflag>1, fprintf('\n');
     if ~isempty(varargin) || ~isempty(s)
        if length(s)>50, s=sprintf(' <%d VARS>',m); end
        fprintf(1,'   save: %s%s%s\n', f,s,sv);
     end
     fprintf(1,'   file: %s\n   path: %s\n\n', f, pwd);
  elseif vflag
     if nx>0
        if length(sx)>40, wblog(1,'I/O skipping %d vars',nx);
        else, wblog(1,'I/O skipping%s',sx); end
     end
  end

  if isbatch
     try
        evalin('caller','Imain.finished=datestr(now);');
        I=getuser(groot,'Imain'); if isfield(I,'cpuhrs')
        evalin('caller',sprintf('Imain.cpuhrs=%g;',cputime/3600-I.cpuhrs)); end
        if vflag>1, evalin('caller','structdisp(Imain);'); end
     catch l
        fprintf(1,[...
         '\n   ERR failed to set Imain.finished time stamp' ... 
         '\n   ERR %s\n'],l.message); dispstack(l);
     end
  end

  wsz='MATLAB:save:sizeTooBigForMATFile';
  tag=regexprep(wsz,'.*:','');

  [msg,id]=lastwarn;
  if ~isempty(findstr(id,tag)), lastwarn(''); end

  % ====================== %
    warning('off',wsz);
    evalin ('caller',cmd); 
    warning('on',wsz);
  [msg,id]=lastwarn;
  if ~isempty(findstr(id,tag))
     try
        s=dbstack(1);
        if ~isempty(s)
             s=sprintf('%s:%g',s(1).name,s(1).line);
        else s='workspace'; end

        fprintf(1,'\n> WRN %s %s; retry ...\n> WRN %s -v7.3\n',...
          s,regexprep(msg,'cannot be saved.*','cannot be saved'),cmd);

        v=regexprep(msg,'.*Variable ''(.*)''.*',' $1'); fprintf(1,'\n');
        if isempty(findstr(sv,'-struct'))
           evalin('caller', ['whos ' v]);
        else
           evalin('caller', ['whos2(''' varargin{2} ''')']);
        end; fprintf(1,'\n');

        evalin('caller', [cmd ' -v7.3']); 
     catch l
        wblog('ERR','check save2 ...'); cmd, sv, varargin
        disp(l.message); dispstack(l);
     end
  end

  if xflag
     if kflag
        cmd=sprintf('touch -r ''%s'' ''%s''',f,ftmp);
        system(cmd);
     end

     cmd='/bin/mv -f'; if vflag>1, cmd=[cmd 'v']; end
     cmd=[cmd sprintf(' ''%s'' ''%s''',ftmp,f)];
     system(cmd);
  end

  if ~isempty(ldir), cd(ldir); end

end

