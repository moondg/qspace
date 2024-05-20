function echo_c(varargin)
% function echo_c([opts,] infostring1, infostring2,...)
%
%    Echo (print) command in this or next non-empty line in
%    calling function or script.
%   
%    Auxilliary logging routine that is a more elaborate
%    alternative to matlab's echo on/off command;
%    uses auxiliary Perl script.
%
% Options
%
%    echo_c -H    % just print header line based in input info strings
%    echo_c -h    % same but without extra separator newline after header
%
% NB! The display of commands can be deactivated as follows
%
%    echo_c -q    % quiet, but keep info strings
%    echo_c -Q    % Quiet (also skips info strings)
%
% Wb,May10,24

  persistent vflag

  hflag=0;

  if isempty(vflag), vflag=2; end
  if nargin && ~isempty(varargin{1}) && isequal(varargin{1}(1),'-')
     done=1;
     if     isequal(varargin{1},'-q'), vflag=1;
     elseif isequal(varargin{1},'-Q'), vflag=0;
     elseif isequal(varargin{1},'-v'), vflag=2;
     elseif isequal(varargin{1},'-s')
        wblog(' * ','currently using vflag=%g',vflag); return
     elseif isequal(lower(varargin{1}),'-h'), done=0;
        if varargin{1}(2)=='H', hflag=1; else hflag=2; end
        varargin(1)=[]; 
     else done=0; end
     if done
        wblog(' * ','switching to mode %s (vflag=%g)',varargin{1},vflag); 
        if nargin>1, varargin(1)=[]; else return; end
     end
  end
  nargs=numel(varargin);
  if ~vflag || vflag<2 && ~nargs, return; end

  [S,I]=dbstack; if numel(S)<2, return; end

  f=S(2).file; F=which(f); if isempty(F), return; end
  l=S(2).line; L=num2str(l); sout=''; e=0;

% auxilliary Perl script to search around current source line caller
  if vflag>1 && ~hflag, this = mfilename; cmd = [
    'open(FH,''<'',''' F ''') or die "' f ': $!"; ' ... 
    'my $l=' L '; my $cstr=""; ' ...
    'while (<FH>) { if ($. == ' L ') { while ($_) { ' ...
       'if (!s/[,;\s]*if.*' this '.*\bend\b[,;\s]*//) { ' ...
       'if (!s/[,;\s]*' this '\(''[^'']*''\)[,;\s]*//) { ' ...
       'if (!s/[,;\s]*' this '\(\[[^\]]*\]\)[,;\s]*//) { ' ...
            's/[,;\s]*' this '[,;\s]*//; }}} ' ...
       'if (/^\s*$/) { $_=<FH>; } ' ... 
       'elsif (/^\s*\%\s*/) { chomp($cstr=$''); $_=<FH>; }' ...
       'else { $l=$.; last; }' ... % skip empty lines
    '}; last; }}; close FH; ' ...
    'if ($_) { chomp; s/^\s*//; print "$l\n$cstr\n$_"; }' ]; ...
   % or { die; } // accept echo_c at end of function or script

     [e,sout]=system([ 'perl -we ''' regexprep(cmd,'''','''\\''''') '''']);
  end

  if e && ~nargs, return; end

  e1=''; e2=''; em='';
  if wblog('--hl-check') % check terminal mode
   % enable high-lighted colored output via escape sequences
   % cf. QSpace documentation Sec. B.8 on Logging
     e1=[char(27) '[34m']; % blue
     e2=[char(27) '[38;5;240m']; % gray
     em=[char(27) '[0m'];
  end

  if nargs, l=72;
     if nargs==0
        istr=varargin{1}; l_=repmat('━',1,max(1,l-4-length(istr)));
        fprintf(1,[e1 '━━ %s %s\n' em],istr,l_);
     else
        l_=repmat('━',1,l-2);
        fprintf(1,[e1 '┏' l_ '┓\n']);
        for i=1:nargs
             fprintf(1,['┃  %-'  num2str(l-4) 's┃\n'],varargin{i});
        end; fprintf(1,['┗' l_ '┛' em '\n']);
     end
     if hflag<2, fprintf(1,'\n'); end
  end

  if ~e && ~isempty(sout)
     i=find(sout==char(10)); cstr={'',''};
     if numel(i)==2 % should always be the case
        L=sout(1:i(1)-1);
        cstr{1}=sout(i(1)+1:i(2)-1);
     elseif isempty(i), i=0;
     end
     sout=sout(i(end)+1:end);

   % strip comment
     sout=regexprep(sout,' *%\s+(.*)(?@cstr{2}=$1;)','');
     for i=1:2 % skip empty comments
        if ~isempty(regexp(cstr{i},'^\s*$')), cstr{i}=''; end
     end
     if ~isempty(cstr{1})
      % print comment with echo_c *before* command
        fprintf(1,[e1 '' e2 '%%  %s' em '\n'],cstr{1});
     end
     if ~isempty(cstr{2})
      % print comment with command *after* command
        fprintf(1,[e1 '>> %-40s ' e2 '%% %s:%s\n' e2 '%%  %s' em '\n'],...
        sout,f,L,cstr{2});
     else
        fprintf(1,[e1 '>> %-40s ' e2 '%% %s:%s' em '\n'],sout,f,L);
     end
  end

end

