function banner(varargin)
% Function banner([wtype,] str ...)
%
%     print text (handed over to sprintf) as banner
%     with wtype setting the flashiness
%       1   most flashy
%       2   intermediate
%       3   less flashy
%       4   least flashy
%       5   color coded
%      '%'  show as matlab comment with leading % (otherwise similar to 4)
%    'box'  draw utf8 box around text shown across entire terminal with
%
% NB! info string is also written to window title

% wtype=2-isbatch;

  wtype=[];

  if nargin>1 && length(varargin{1})<4
     if isnumber(varargin{1}), wtype=varargin{1}; varargin(1)=[];
     elseif ischar(varargin{1})
        if ~isempty(regexpi(varargin{1},'^box$'))
           wtype=varargin{1}; varargin(1)=[];
        end
     end
  end
  if isempty(wtype), wtype=1; end

  persistent tlast

  nargs=numel(varargin);
  if nargs && ~ischar(varargin{1}), wbdie('invalid usage'); end

  if nargs>1,   sout=sprintf(varargin{:});
  elseif nargs, sout=varargin{1};
  else sout=''; end

  L=getcols(); if L<60 || L>160, L=85; end

  if isequal(wtype,5)
     nl=char(10); 
     if isdesktop>1, e=char(27);
          wesc_={ [e '[32m'], [e '[0m'] };
     else wesc_={'',''}; end

     fprintf(1,[ nl ...
        wesc_{1} '>> ' sout ' ' repmat('-',1,L-length(sout)-4) ...
        wesc_{2} nl nl]);
     return
  end

  lsep=''; istr='  '; 
  use_box=0;

  switch wtype
    case 1
      if ~isempty(tlast)
         dt=now-tlast; if dt<1
              s=['  [' datestr(dt,   'HH:MM:SS') ']'];
         else s=['  [' datestr(dt,'dd-HH:MM:SS') ']']; end
      else s=''; end
      lsep=repmat('*',1,L); tlast=now;
      fprintf(1,'\n\n%s\n\n  %s/%s\n  %s%s\n',...
         lsep, hostname, pwd, datestr(now),s);
      if ~isempty(sout)
         lsep=['\n' lsep '\n\n']; fprintf(1,'\n');
     end
    case 2
      lsep=repmat('*',1,L);
      fprintf(1,'\n\n%s\n\n  ',lsep);
      lsep=['\n' lsep '\n\n'];
    case 3
      lsep=repmat('=',1,L);
      fprintf(1,'\n%s\n  ',lsep);
      lsep=[lsep '\n\n'];
    case 4
      lsep=repmat('-',1,L);
      fprintf(1,'\n%s\n  ',lsep);
      lsep=[lsep '\n\n'];
    case '%'
      lsep=repmat('-',1,L-2);
      fprintf(1,'\n%% %s\n%% ',lsep);
      lsep=['%% ' lsep '\n\n'];
    case {'box','Box','BOX'}
      if     isequal(wtype,'BOX'), use_box=3; wtype='box';
      elseif isequal(wtype,'Box'), use_box=2; wtype='box';
      else                         use_box=1; end
      if use_box==1
           bb={'─','│','┌' '┐','└','┘'};
      else bb={'━','┃','┏','┓','┗','┛'}; end

      if use_box>2, S=dbstack(); 
         if numel(S)>1, S=S(2);
            sout=[sprintf('%s:%d\n', S.file,S.line), sout];
         end
      end

      lsep=repmat(bb{1},1,L-2);
      fprintf(1,['\n' bb{3} lsep bb{4} '\n']);
      istr={[bb{2} ' '],bb{2}};
      lsep=[bb{5} lsep bb{6} '\n\n'];

    otherwise,
  end

  if use_box
    if ~isempty(regexp(sout,'\x27|\\e\>'))
       use_box=0; istr=istr{1};
    end
  end

  if isequal(wtype,'%'), i=' % '; i(1)=10;
     fprintf(1,[istr '%s\n'],strrep(sout,i(1),i));
  elseif use_box
     sout=regexprep(sout,'\\n\>\s*','\n');
     ss=textscan(sout,'%s','whitespace','\n'); ss=ss{1};

     fmt0=[istr{1} sprintf('%%-%ds',L-3) istr{2} '\n'];

     for i=1:numel(ss), [~,dn]=strlen(ss{i});
        if dn
             fmt=[istr{1} sprintf('%%-%ds',L-3+dn) istr{2} '\n'];
        else fmt=fmt0; end
        fprintf(1,fmt,ss{i});
     end
  else
     sout=regexprep(sout,'(\n|\\n\>)',['\n' istr]);
     fprintf(1,[istr '%s\n'],sout);
  end

  if ~isempty(lsep), fprintf(1,lsep); end
  if wtype==1, inl 1; end

% write info also to window title
% if ~isbatch && wtype<=2
%  % NB! might fail if `sout' contains carriage return or backslashes!
%    sout=strrep(sout,'\','');
%  % sout=strrep(sout,char(10),'; ');
%    i=find(sout==char(10)); if ~isempty(i), sout=sout(1:i(1)-1); end
%
%    try 
%        eval(['! echo -ne "\033]0;matlab @  ${HOST}  ' ...
%        regexprep(repHome(pwd),'\$','\\$') ' : ' sout '\007"'])
%    catch
%    end
% end

end

