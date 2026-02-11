function [rval,wesc]=wblog(varargin)
% function wblog([stack-level,] fmt,...)
%
%    logging routine with output format:
%    filename:linenr (output as with sprintf(fmt,...))
%
%    If stack-level>0 is specified (default: 0),
%    this shows filename:linenr from respective higher dbstack entry.
%
% Format string
%
%    \n  leads to newline with usual file:line header
%    \N  leads to old fashioned newline, i.e. with no header
%
%    except for the very first appearances of \N or \n in fmt,
%    which are interpreted as global leading \n
%
% Wb,Nov18,05: allow to redirect output to file
%
%    wblog('setfid',fid)     - redirects the following wblog's to file fid
%    wblog('setfout',fname)  - redirects the following wblog's to file fname
%    wblog('unsetfid')
%    wblog('unsetfout')      - unset the output to file
%
% Based on lineno() from C.Denham which uses dbstack()
% for more info: see help to lineno()
%
% Wb,Oct08,03

  persistent fid newfile llday

  if ~nargin
     if ~isdeployed, eval(['help ' mfilename]); end
     return
  end
  hl_check=0;

  if nargin && isequal(varargin{1},'--hl-check')
     rval=wblog_hl_check();
     if nargin<=1 || ~ischar(varargin{2}), return; end
     hl_check=1;

  elseif nargin==1
     if isequal(varargin{1},'--ping')
        llday=check_nextday(llday,fid); return
     elseif ~isempty(regexp(varargin{1},'^--got-')), q=varargin{1}(7:end);
        u=getuser(0,'err_count'); rval=0;
        if     regexpi(q,'^ERR$'), if isfield(u,'err'), rval=u.err; end
        elseif regexpi(q,'^WRN$'), if isfield(u,'wrn'), rval=u.wrn; end
        elseif regexpi(q,'^TST$'), if isfield(u,'tst'), rval=u.tst; end
        else wbdie('invalid usage (%s)',q); end
        return
     elseif isequal(varargin{1},'unsetfid')
        fid=[]; return
     elseif isequal(varargin{1},'unsetfout')
        if ~newfile
           wblog('ERR','fid set by calling routine.');
        elseif ~isempty(fid)
           fclose(fid);
        else wblog('ERR','No file to close.'); end

        fid=[]; return
     end

  elseif nargin==2

     if isequal(varargin{1},'setfid')

        if ~isempty(fid)
        wblog('ERR','must close active fid first.'); return; end

        fid=varargin{2}; newfile=0;
        try, fprintf(fid,'');
        catch,fid=[]; wblog('ERR','Invalid fid.'); end
        return

     elseif isequal(varargin{1},'setfout')

        if ~isempty(fid)
        wblog('ERR','must close active fid first.'); return; end

        fname=varargin{2};
        fid=fopen(fname,'w'); newfile=1;
        try, fprintf(fid,'');
        catch,fid=[]; wblog('ERR','Cannot write to file %s.',fname); end
        return

     end
  end

  if isnumber(varargin{1})
       idx=varargin{1}; varargin(1)=[];
  else idx=0; end

  s=lineno_aux(2+idx);
  t=datestr(now);

  tag=''; TAG=''; nargs=numel(varargin);

  if hl_check, tag=varargin{2};
  elseif nargs>1 && ischar(varargin{2}) && isempty(find(varargin{1}=='%'))
     tag=varargin{1}; varargin(1)=[];
  elseif nargs
     [i,j,x]=regexp(varargin{1},'^\s*(.)\s+');
     if ~isempty(i) && x{1}(1)<=3 && j(1)<5
        TAG=sprintf('LL%d',x{1}(1)); tag=varargin{1}(1:3);
        varargin{1}=varargin{1}(5:end);
     else
        l=regexp(varargin{1},'\s');
        if ~isempty(l) && l(1)<8, l=l(1); q=varargin{1}(1:l-1);
           if isempty(strfind(q,'%'))
              [i,j,x]=regexp(q,'^(WRN|ERR|LL\d|NB!|CST|Error|Warning)');
              if ~isempty(i), x=x{1};
                 tag=q(x(1):x(2));
                 varargin{1}=varargin{1}(l+1:end);
              end
           end
        end
     end
  end

  cstflag=0; if isempty(TAG), TAG=tag; end

  if ~isempty(tag)
     switch tag
       case {'LL1'}, tag='*  ';
       case {'LL2'}, tag=' * ';
       case {'LL3'}, tag='  *';
       case {'CST'}, cstflag=1; tag='';
       case {'Error'}, tag='ERR';
       case {'Warning'}, tag='WRN';
     end
     if ~cstflag, tag=sprintf('%3s ',tag); end

     if ~isempty(tag)
        u=get(0,'UserData'); f=deblank(lower(tag));
        if isempty(u) || ~isfield(u,'err_count')
        set_global; u=get(0,'UserData'); end

        g=u.err_count;
        if isfield(g,f)
           g=setfield(g,f, getfield(g,f)+1);
           u.err_count=g;
           set(0,'UserData',u);
        end
     end
  end

  [use_col,iterm]=wblog_hl_check();
  wesc={0,'',''};

  if isempty(fid) && iterm
     if     ~isempty(regexp (tag,'WRN')), wesc={1, 9};
     elseif ~isempty(regexp (tag,'ERR')), wesc={2, 1};
     elseif ~isempty(regexp (tag,'NB!')), wesc={4,34};
     elseif ~isempty(regexpi(tag,'ok!')), wesc={4, 2};
     elseif ~isempty(regexpi(tag,'^(ok\.|==>|SUC|<i>)\s*$'))
          wesc={4,{150,200,150}};
     end
  end

  werr=bitand(wesc{1},3);

  if wesc{1} && ~use_col
     wesc={0,'',''};
  end

  if wesc{1}
     esc_=[char(27) '['];

     if numel(wesc{2})==1
        if wesc{2}<8
             wesc{2}=[esc_ sprintf('3%dm',wesc{2})];
        else wesc{2}=[esc_ sprintf('38;5;%dm',wesc{2})];
        end
     elseif numel(wesc{2})==3
          wesc{2}=[esc_ sprintf('38;2;%d;%d;%dm',wesc{2}{:})];
     else wesc
        wbdie('invalid wesc');
     end
     wesc{3}=[esc_ '0m'];
  end

  if hl_check
     wesc=wesc([2 3 1]);
     return
  end

  hstr=[ sprintf('%-20.20s %s  ', shortfstr(s,20), t(13:end)) wesc{2} tag ];
  fmt=varargin{1};

  n=0; k=0; i=1; l=length(fmt);
  while i<=l
     if fmt(i)==char(10), i=i+1; k=i; n=n+1;
     elseif i<l
        if (isequal(fmt(i:i+1),'\n') ...
         || isequal(fmt(i:i+1),'\N')), i=i+2; k=i; n=n+1;
        else break; end
     else break; end
  end
  if n, fprintf(1,repmat('\n',1,n)); end
  if k, fmt=fmt(k:end); end

  nl=char(10); NL=[ wesc{3} nl hstr ];

  fmt=strrep(fmt, nl, NL);
  fmt=strrep(fmt,'\n',NL);
  fmt=strrep(fmt,'\N',nl); q=length(fmt)>1;

  if     q && isequal(fmt(1:2),'\\'), fmt=[wesc{2} fmt(3:end)];
  elseif q && isequal(fmt(1:2),'\r'), fmt=['\r' hstr fmt(3:end)];
  else fmt=[ hstr fmt ]; end

  if isequal(fmt(max(1,end-1):end),'\\'), fmt=[ fmt(1:end-2) wesc{3} ];
  else fmt=[ fmt wesc{3} nl ]; end

  llday=check_nextday(llday,fid);

  if isempty(fid)
     fprintf(1,fmt,varargin{2:end});
     if wesc{1}, fprintf(1,wesc{3}); end
  else
     fprintf(fid,fmt,varargin{2:end});
  end

  if ~werr, return; end

% -------------------------------------------------------------------- %
% also respond to dbwrn and dberr // Wb,May31,19
  q=werr;
  if q
     if     (q&1) && ~isset(getuser(0,'dbstop_if_WRN')), q=0;
     elseif (q&2) && ~isset(getuser(0,'dbstop_if_ERR')), q=0; end
  end
  if ~q
     if ~isempty(who('global','DEBUG')), global DEBUG
     else
        DEBUG=getenv('DEBUG'); if isempty(DEBUG)
        DEBUG=getenv('ML_DEBUG'); end
     end
     if ~isempty(DEBUG), q=werr; end
  end

  if ~q, return; end

  S=dbstatus; if isempty(S), return; end
  q_=0; for i=1:numel(S)
     if     isequal(S(i).cond,'warning'), q_=bitor(q_,1);
     elseif isequal(S(i).cond,'error'  ), q_=bitor(q_,2); end
  end

  if bitand(q,q_)
   % use regular keyboard stop, so code can be resumed (e.g. after WRN)
   % then use e.g. dbquit to *not* resume program // Wb,May31,19
   % if bitand(q,2), inl(1); dbstack; inl(1); keyboard % beep; 
   % elseif q,       inl(2); dbstack; inl(1); keyboard % beep; 
   % end
     wbstop
  else
     if bitand(q,2) inl(2); error('Wb:ERR',''); end
  end
end

% -------------------------------------------------------------------- %
% NB! since wblog() is called within getopt avoid
% wblog() to use routines with varargin calling getopt()

function lstr = lineno_aux(id)

  [S,index]=dbstack;
  if 1+id>numel(S), lstr='(cmdline)'; return; end

  S=S(1+id);

  if ~isequal(S.file, [S.name '.m']);
       lstr = sprintf('%s>%s:%d', S.file, S.name, S.line);
  else lstr = sprintf('%s:%d', S.file, S.line); end

end

% -------------------------------------------------------------------- %
% check whether to use colored output
% outsourced from main // Wb,Jul13,23

function [q,iterm]=wblog_hl_check()
   q=0;

   iterm=isdesktop()>1; if ~iterm, return; end

   qs=getenv('QS_LOG_COLOR');
   if isempty(qs), q=1;
   else q=str2num(qs);
      if isempty(q) || numel(q)~=1 || q<0
         wbdie(-2,'invalid QS_LOG_COLOR = %s',qs); 
      elseif q<0, q=0; end
   end
end

% -------------------------------------------------------------------- %
% check day of last log

function llday=check_nextday(llday,fid)

  d=datevec(now); d=d(3);
  if ~isequal(d,llday)
     if ~isempty(llday)
        if isempty(fid), f=1; else f=fid; end
        fprintf(f,'\n\n>> TODAY %s :: (%d->%d)',datestr(now),llday,d);
        fprintf(f,'\n\n');
     end
     llday=d;
  end

end

% -------------------------------------------------------------------- %
