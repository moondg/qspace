function [i,fname]=fexist(fname,tflag)
% Function: [i,fname]=fexist(fname [,opts])
%
%    If file exists and shall not be overwritten, 1 is returned,
%    0 otherwise. In interactive mode, alternative file name can be
%    specified, which is returned as second argument.
%
% Options
%
%   '-t'  test flag
%
% See also isfile2.m for plain check whether file exists (eg. for reading only)
% Wb,May20,07

  i=0;

  if exist(fname,'file')
     [p,n,x]=fileparts(fname); if isempty(p), p='.'; end
     if isbatch
        wblog(1,'WRN','overwriting %s', repHome(fname));
     else
        q=sprintf(['\n' ...
          '=> file %s exists.\n' ...
          '   overwrite? { [1y]|0n| new filename } '], repHome(fname));
        if nargin<2 || ~tflag, q=input(q,'s'); inl(1);
        else fprintf(1,'%s\n',q(1:findstr(q,'s.'))); q=0; return; end
        if isempty(q); return; end

        if isequal(lower(q),'keyboard') || isequal(lower(q),'dbstack')
           keyboard; [i,fname]=fexist(fname); return
        elseif isequal(lower(q),'y'), i=0; return
        elseif isequal(lower(q),'n'), i=1; return
        end

        [x,ok]=str2num(q);
        if ~ok
           if nargout<2, wberr('file name provided is ignored'); end
           [i,fname]=fexist([p '/' q]);
        else i=~x; end
     end
  end

end

