function s=wbstamp(varargin)
% function s=wbstamp(varargin)
% Wb,Jun18,07

% see also wbstamp.pl

  getopt('INIT',varargin);
     lflag=getopt('-l'); if ~lflag && getopt('-L'), lflag=2; end
     fflag=getopt('-f');
     bflag=getopt('-b');
     sflag=getopt('-s');
     pflag=getopt('-p');
     iflag=getopt('-init');
     cflag=getopt('-clear');
  getopt('check_error');

  if cflag, getuser(groot,'wbstamp','-rm'); end
  if ~iflag
     s=getuser(groot,'wbstamp');
     if ~isempty(s) && ischar(s) && ~fflag
        wblog(' * ','using wbstamp=''%s''',s);
        return
     end
  end

  s=['Wb', datestr(now,'yymmdd') ];
  if bflag, return; end

  if pflag
     s=[s '-' sprintf('%03x',mod(getpid,4096))];

  elseif lflag
     if lflag==1,
          s=[s '_' datestr(now,'HHMM')];
     else s=[s '_' datestr(now,'HHMMSS')]; end

  elseif ~sflag, [jid,tid]=getJobID();
     if ~isempty(jid)
        if isempty(tid)
             s=[s, sprintf('_%d',jid)];
        else s=[s, sprintf('_%d.%02d',jid,tid)]; end
     else
        c=clock; s(end+1)='a'+c(4);
     end
  end

  if iflag
     wblog(' * ','setting persistent wbstamp=''%s''',s);
     setuser(groot,'wbstamp',s);
     if ~nargout, clear s, end
  end

end

