function sysinfo(varargin)
% Function: sysinfo([whos])
%
%    Example:
%    sysinfo([whos; whos('global')]);
%
% See also getmem.m
% Wb,Aug01,07

  getopt('init',varargin);
     wflag=getopt('-w');
  varargin=getopt('get_remaining'); narg=length(varargin);

  l=repmat('=',1,82);
  fprintf(1,'\n%s\n',l);

  system('t2 -du -nl 3'); gotw=1;

  if ~narg
     if wflag
        evalin('caller','setuser(groot,''whos'',[whos; whos(''global'')]);');
        w=getuser(groot,'whos','-rm');
     else w=[]; end
  else w=varargin{1};
     if ~isfield(w,'bytes') || narg>1
     w=[]; wblog('ERR','invalid usage'); end
  end

  if ~isempty(w)
     mark=zeros(size(w));
     for i=1:numel(w), if w(i).bytes>1E6, mark(i)=1; end, end
     i=find(mark); inl(1);
     if ~isempty(i)
        evalin('caller',['whos ' strhcat({w(i).name})]);
     end

     fprintf(1,'  Total of %s (%g variables).\n',...
     num2str2(sum(cat(1,w.bytes)),'--bytes'), numel(w));
  end

  eval(sprintf('! pstime.pl %g',getpid));

  fprintf(1,'%s\n\n',l);

end

