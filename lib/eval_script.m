function eval_script(f)
% function eval_script(f)
%
%    Evaluate script by reading it as text and running eval on it (in caller)
%    Intended for mcc-batch jobs which cannot run additional un-compiled
%    parameter scripts on demand.
%
% Wb,Jul16,19

  f_=f;
  if isempty(regexp(f,'\.m')), f=[f '.m']; end
  if ~exist(f)
     if isempty(regexp(f,'^[\.\/]')), f=['./' f]; end
  end

  fid=fopen(f,'r'); ll={};
  if ~fid, error('Wb:ERR','\n   ERR invalid file ''%s''',f_); end

  while 1
     l=fgetl(fid); if ~ischar(l), break; end
     ll{end+1}=l;
  end
  fclose(fid)

  l=sprintf('%s\n',ll{:});
  fprintf(1,'>> %s (%d lines)\n',f_,numel(ll));

  evalin('caller',l);

end

