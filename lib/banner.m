function banner(varargin)
% Function banner([flag,] str ...)
%
%     print text (handed over to sprintf) as banner
%     with flag setting the flashiness:
%       1   most flashy
%       2   intermediate
%       3   less flashy
%       4   least flashy
%       5   color coded
%      '%'  show as matlab comment with leading % (otherwise similar to 4)
%
% NB! info string is also written to window title

% flag=2-isbatch;
  if nargin>1 && isscalar(varargin{1})
       flag=varargin{1}; varargin(1)=[];
  else flag=1; end

  if isempty(varargin) || ~ischar(varargin{1}), helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  persistent tlast

  if numel(varargin)>1
       s=sprintf(varargin{:});
  else s=varargin{1};
  end

  L=getcols(); if L<74, L=74; end

  if isequal(flag,5)
     nl=char(10); 
     if isdesktop>1, e=char(27);
          wesc_={ [e '[32m'], [e '[0m'] };
     else wesc_={'',''}; end

     fprintf(1,[ nl ...
        wesc_{1} '>> ' s ' ' repmat('-',1,L-length(s)-4) ...
        wesc_{2} nl nl]);
     return
  end

  switch flag
    case 1
      if ~isempty(tlast)
           s=now-tlast; if s<1
                s=sprintf('  [%s]',datestr(s,   'HH:MM:SS'));
           else s=sprintf('  [%s]',datestr(s,'dd-HH:MM:SS')); end
      else s=''; end
      lsep=repmat('*',1,L); tlast=now;
      fprintf(1,'\n\n%s\n\n  %s/%s\n  %s%s\n\n  ',...
      lsep, hostname, pwd, datestr(now),s);
    case 2
      lsep=repmat('*',1,L);
      fprintf(1,'\n\n%s\n\n  ',lsep);
    case 3
      lsep=repmat('=',1,L);
      fprintf(1,'\n%s\n  ',lsep);
    case 4
      lsep=repmat('-',1,L);
      fprintf(1,'\n%s\n  ',lsep);
    case '%'
      lsep=repmat('-',1,L-2);
      fprintf(1,'\n%% %s\n%% ',lsep);
    otherwise,
  end

  if isequal(flag,'%'), i=' % '; i(1)=10;
     fprintf(1,'%s\n',strrep(s,i(1),i));
  else
     s=regexprep(s,'\n([^#%])','\n  $1');
     fprintf(1,'%s\n',s);
  end

  switch flag
    case  1,  fprintf(1,'\n%s\n\n',lsep);
    case  2,  fprintf(1,'\n%s\n\n',lsep);
    case  3,  fprintf(1,'%s\n\n',lsep);
    case  4,  fprintf(1,'%s\n\n',lsep);
    case '%', fprintf(1,'%% %s\n\n',lsep);
  end

  if flag==1, inl 1; end

% write info also to window title
% if ~isbatch && flag<=2
%  % NB! might fail if `s' contains carriage return or backslashes!
%    s=strrep(s,'\','');
%  % s=strrep(s,char(10),'; ');
%    i=find(s==char(10)); if ~isempty(i), s=s(1:i(1)-1); end
%
%    try 
%        eval(['! echo -ne "\033]0;matlab @  ${HOST}  ' ...
%        regexprep(repHome(pwd),'\$','\\$') ' : ' s '\007"'])
%    catch
%    end
% end

end

