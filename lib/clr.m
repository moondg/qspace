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
% NB! use `clear classes' when constructors are changed // Wb,Ap21,06
% NB! `clear functions' also clears dbstop's

  if ~keep
     if ~figf
     evalin('caller','close all'); else untagf -a -q; end
     evalin('caller','clear functions; dbclear all; clear classes');
  end

  evalin('caller','clear all; clear global');

  u=get(0,'UserData');
  [u,i]=clear_fields(u,'err_count','tstEnv');
  if i, set(0,'UserData',u); end

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
  clc

end

% -------------------------------------------------------------------- %
% initialize fields (numbers to 0, cells to empty)
% where the return value i counts the number of changes
% Wb,May06,25

function [u0,i]=clear_fields(u0,varargin)

  i=0;
  for k=1:numel(varargin), f0=varargin{k};
     if isfield(u0,f0), u1=getfield(u0,f0);
        if isstruct(u1), j=0;
           for f1=fieldnames(u1)', f1=f1{1};   q =getfield(u1,f1);
              if isnumber(q) && q,             u1=setfield(u1,f1, 0); j=j+1;
              elseif iscell(q) && ~isempty(q), u1=setfield(u1,f1,{}); j=j+1;
              end
           end
           if j,                           u0=setfield(u0,f0,u1); i=i+j; end
        elseif isnumber(u1) && u1,         u0=setfield(u0,f0, 0); i=i+1;
        elseif iscell(u1) && ~isempty(u1), u0=setfield(u0,f0,{}); i=i+1;
        end
     end
  end

end

% -------------------------------------------------------------------- %

