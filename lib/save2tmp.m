function save2tmp(varargin)
% Function: save2tmp([usual save without filename])
%
%    Save variables to temporary file $ML/tmp.mat to exchange data
%    between MatLab sessions. Arguments are the same as for MatLab's
%    save except that no file name is specified.
%
% Options
%
%   '-<#>'    temporary file with id <#> (#=1..9)
%   '-tid' #  same as above
%   '-a'      append to file
%   '-db'     debug mode: create independent temfile in local directory.
%   '-x',..   exclude pattern
%
%   All remaining arguments are handed over to MatLabs save().
%
% See also loadtmp.m, save2.m
% Wb,Nov12,09

  getopt('init',varargin);
     aflag =getopt('-a');
     dbflag=getopt('-db');
     tid   =getopt('-tid','');
     xpat  =getopt('-x','');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg && ischar(varargin{1}) && ~isempty(regexp(varargin{1},'^-[0-9]+$'))
     tid=varargin{1}(2:end);
     varargin(1)=[]; narg=narg-1;

     q=str2num(tid); if q<0 || q>9
     wbdie('invalid tid=%s',tid); end
  end

  if dbflag || isbatch()
     s=dbstack; s=s(min([2, length(s)]));
     f=[ pwd '/' wbstamp('-l') sprintf('_%s_%04d_debug',s.name,getpid) '.mat' ];
  else
     f=[ getenv('HOME') '/Matlab/tmp' tid '.mat' ];
  end

  if ~isempty(xpat)
     if isempty(varargin), varargin=evalin('caller','who'); end
     i=find(cellfun(@numel,regexp(varargin,xpat)));
     if ~isempty(i)
        s=sprintf(', %s',varargin{i});
        if numel(i)>1, s(1:2)='s '; else s=s(2:end); end
        fprintf(1,'\n   skipping var%s\n',s);
        varargin(i)=[];
     end
  end

  cmd=['save ' f ' ' sprintf(' %s', varargin{:}) ];

  if aflag, cmd=[cmd ' -append']; end

  evalin('caller', cmd); 

  if ~isempty(varargin)
       s=sprintf(' %s', varargin{:}); s=['(' s(2:end) ')'];
  else s=''; end
  fprintf(1,'\n   data saved to %s%s\n\n',basename(f),s);

end

