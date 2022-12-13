function vi(fname0)
% function vi(fname) -  opens given filename with vi (UNIX)
%
%    fname just must be given such that MatLab is able to find it
%    with the which command (this means file extension does not
%    need to be there necessarily
%
% Wb,Apr04,01 ; Wb,Sep27,12

  if nargin~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  i=regexp(fname0,':\d+$');
  if ~isempty(i)
       q=fname0(i:end); fname0=fname0(1:i-1);
  else q=''; end

  fname=which(fname0);
  if isempty(fname)
     if ~strcmpi(fname0(end-1:end), '.m'), fname0 = [fname0, '.m']; end
     fname=fname0;
  end

  if ~isempty(strfind(fname,'built-in'))
     fname=regexprep(fname,'.*\((.*)\)','$1.m');
  elseif ~isempty(q)
     fname=[fname q];
  end
  
  if isunix
     cmd = ['! unset LD_LIBRARY_PATH ; vim ' fname];
     eval(cmd)
  else
     cmd = ['edit ' fname];
     eval(cmd)
  end

end

