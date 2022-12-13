function check_interrupt(fname,varargin)
% function check_interrupt(fname)
%
%    Dynamic parameter change of running (cluster) job by checking
%    whether specified script <fname> exists. If it does, it is
%    executed once, which is ensured by subsequently archiving it
%    by appending current time stamp in the format '-dd-HHMM'
%    to the filen name.
%
%  Possible uses
%
%  - change parameters on the fly for running (cluster) jobs
%  - even temporarily stop a job, e.g.  by putting them into a
%    temporary wait loop or to perform a particular data analysis
%    on the current data set that save to files etc.
%
% Wb,Mar10,21

% E.g. see tst_Hamilton1D which call *this at the end of a full sweep.

  getopt('init',varargin);
     tflag=getopt('-t');
  getopt('check_error');

  q=dir(fname);
  if isempty(q)
     if tflag, wblog(' * ','failed to find script %s',repHome(fname)); end
     return
  end

  fprintf(1,'\n>> running script ./%s\n',q.name);
  system(['perl -we $'' ' ... 
    'open(FH,\''<\'',\''' q.name '\'') || die "invalid file"; ' ...
    'foreach (<FH>) { if (!/^\s*$/) { printf " | %s", $_; }}; ' ... 
    'close(FH); '' '
  ]);

  F=[q.folder '/' q.name];
  fid=fopen(F,'r'); ll={};
  if fid
     while 1
        l=fgetl(fid); if ~ischar(l), break; end
        ll{end+1}=l;
     end
     fclose(fid);
  end

  if ~isempty(ll)
     l=sprintf('%s\n',ll{:});
     try evalin('caller',l);
     catch me
        wblog('ERR','script ''%s'' failed to run\nmsg: %s',...
           repHome(fname),me.message); 
        dispstack(me.stack);
        disp(me)
     end
  end

  s={'','', datestr(now,'dd-HHMM') };
  if tflag, s{2}='echo TST '; else s{1}='mv '; end

  fprintf(1,'>> %s',s{1});
  Fx=[ regexprep(F,'\.m$','') '-' s{3} '.m' ];
  system([s{2} 'mv -v ''' F ''' ''' Fx '''']);
  fprintf(1,'\n');

end

