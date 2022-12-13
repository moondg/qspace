function clr(varargin)
% function clr([opts])
 %
%    clear entire workspace (and terminal window).
%
% Options
%
%   -k  keep figures, classes, debug stops
%
% Wb,Jan05,07

  try
  getopt ('init', varargin);
     keep = getopt ('-k');
     figf = getopt ('-f');
     eflag= getopt ('-e');  % call subsequent 'dbstop if error'
     rflag= getopt ('-r');
  varargin=getopt('check_error');
  catch
     return
  end

% ----------------------------------------------------------------------- %
  if ~keep
     if ~figf
     evalin('caller','close all'); else untagf -a -q; end
     evalin('caller','clear functions; dbclear all; clear classes');
  end

  evalin('caller','clear all; clear global');

  if eflag, dbstop if error; end
  if rflag, system('reset'); end

  if keep || figf, untagf -a -q; return; end
  if rflag, return; end

  L=str2num(getenv('LINES'));
  if ~isempty(L), if L<80, L=80; end
  else
     L=evalc('!tput cols 2>&1');
     if ~isempty(L)
        L=str2num(L); if L>90, L=90; end
     end
     if isempty(L), L=80; end
  end

  s='   C L E A R   A L L';

  if isdesktop>1, e=char(27);
       wesc_={ [e '[48;5;15;38;5;0m'], [e '[0m'], [e '[38;5;0m' '>>'] };
  else wesc_={'','',''}; end % print '>>' black, and hence invisible ^1)

  l=[0 length(s)]; l(3)=L-sum(l);

  s=[ '\n' wesc_{1} repmat(' ',1,l(1)), s, ...
      repmat(' ',1,l(3)) wesc_{2} '\n' ...
      wesc_{3} wesc_{2} '\n'
  ];

  fprintf(1,s);
  clc; return

end

