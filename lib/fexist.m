function [i,fname]=fexist(fname,varargin)
% Function: [i,fname]=fexist(fname [,opts])
%
%    If file exists in current directory or specified with path.
%
% Options
%
%    -t   test flag
%
%    -W   same as -w, but in batch mode writes WRN overwrite.
%    -w   checking whether file exists prior to (over-)writing it
%         in non-batch mode, this returns i<0 if file exists.
%         within interactive desktop, this asks whether file
%         should be overwritten. If [y1] was answered,
%         this returns i<0, indicating that file shall be overwritten.
%         Instead, one may specify also an alternative file name here
%         which is returned as second argument then.
%
%         Hence with opion -w, a return value <0 indicates
%         that file name exists but shall be overwritten.
%
% See als o isfile2.m for plain check whether file exists (eg. for reading only)
% Wb,May20,07

% [02/24/2024] -w option added to suppress `WRN file is overwritten'
% for plain call fexist(file)

  if ~nargin
     if ~helpthis(nargout), wbdie('invalid usage'); end
     return
  end

  if ~ischar(fname), wbdie('invalid usage'); end
  if isempty(fname)
     wbdie('invalid usage (got empty string for file name)');
  end

  if isempty(regexp(fname,'\/')), fname=['./' fname]; end

  i=exist(fname,'file');
  if isempty(regexp(fname,'\.[a-z]{1,3}')), x='';
     if i
        q=which(fname,'-all'); n=numel(q);
        [~,~,x]=fileparts(q{1});
        fname=[fname x];
        if n>1
           wblog('WRN','%d matches found on path (assuming *%s)',n,x);
        end
     else
        fname=[fname '.mat'];
        i=exist(fname,'file');
     end
  end

  if ~i || nargin<2, return; end

  getopt('INIT',varargin);
     tflag=getopt('-t');
     if getopt('-w'), wflag=1;
     elseif getopt('-W'), wflag=2; else wflag=0; end

  if ~tflag
     tflag=getopt('get_last',0);
     if ~isnumeric(tflag), wbdie('invalid usage (tflag)'); end
  else getopt('check_error'); end

  if ~wflag, return; end

  if isbatch
     if wflag>1
        wblog(1,'WRN','overwriting %s', repHome(fname));
     end
     i=-i; return
  end

  [p,n,x]=fileparts(fname); if isempty(p), p='.'; end
  q=sprintf(['\n' ...
    '=> file %s exists.\n' ...
    '   overwrite? { [1y]|0n| new filename } '], repHome(fname));

  if tflag
     fprintf(1,'%s\n',q(1:strfind(q,'s.')+2));
     return
  end

  if wblog('--hl-check')
     dfmt=[char(27) '[32m$1' char(27) '[0m'];
     q=regexprep(q,'\[([^\]]*)\]',dfmt);
  end

  q=input(q,'s'); inl(1);
  if isempty(q); i=-i; return; end

  if isequal(lower(q),'keyboard') || isequal(lower(q),'dbstack')
     keyboard; [i,fname]=fexist(fname); return
  elseif isequal(lower(q),'y'), i=-i; return
  elseif isequal(lower(q),'n'), return
  else
     [x,ok]=str2num(q);
     if ok, if x, i=-i; end
     elseif nargout<2, wbdie('file name provided is ignored');
     else [i,fname]=fexist([p '/' q]);
     end
  end

end

