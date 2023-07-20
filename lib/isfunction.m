function [i,istr]=isfunction(cmd)
% function i=isfunction(cmd)
%
%    check whether input is a function
%    if it is a function handle, i=9 is returned
%    otherwise it is assumed that cmd is a string that refers
%    to a function or a matlab script (this checks the header
%    line for the keyword 'function').
%
% Return values for i
% this (mostly) uses the return codes of matlab's exist() function
%
%     0  not a function
%    -2  matlab scripts
%     2  m-file on the matlab search path, yet also when path is included
%        this value is set to -2 if the corresponding fiel does not have
%        the keyword `function' in the heaader of its source code
%     3  mex-file
%     5  matlab-builtin
%
%     9  input is a function handle
%
% Wb,Jul13,11

% Wb,Dec07,19: changed i=0 -> -2 for matlab scripts

  if nargin~=1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if ~ischar(cmd)
     if isa(cmd,'function_handle')
        i=9; return
     else wbdie('invalid usaage');
     end
  end

  i=exist(cmd,'file'); istr='';

  if i==3 || i==5, return; end
  if i~=2, i=0; return; end

  w=which(cmd,'-all'); if isempty(w), i=0; return; end
  if numel(w)>1
     disp(strvcat(w{:}));
     wblog('WRN','m-file ''%s'' exists %g times!',cmd,i);
  end

  fid=fopen(w{1},'r'); l=fgetl(fid); fclose(fid);
  if isempty(regexp(l,'^\s*function')), i=-2; end

end

