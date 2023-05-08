function dispstack(S,k)
% Function dispstack(S [,k])
% Wb,Aug24,07

  if ~nargin
     S=dbstack; S=S(2:end);
  elseif nargin>2
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  if isempty(S)
     fprintf(1,'\n   (empty dbstack)\n\n');
     return
  end

  if nargin<2, k=-1; end

  if isobject(S), S=obj2struct(S); end

  if ~isfield(S,'file')
     if isfield(S,'stack')
        if isfield(S,'message'),    fprintf(1,'\n  Message: %s\n', S.message); end
        if isfield(S,'identifier'), fprintf(1,'       id: %s\n',S.identifier); end
        S=S.stack;
     else wberr('invalid usage'); end
  end

% indicator for current frame // Wb,Apr24,13
% (see bt.m)
% if ~isfield(S,'curr'), S(1).curr=0; end
% for i=1:numel(S)
%     if isempty(S(i).curr), S(i).curr=0; end
% end

  n=numel(S); vflag=0; inl 1

  if vflag
     fprintf(1,'\nDBSTACK\n\n');
     for i=n:-1:1
         s=sprintf('%s:%d', S(i).file, S(i).line);
         fprintf(1,'%4d  %-40s %s\n', n-i+1,s,S(i).name);
     end
  else
     for i=n:-1:1
         f=S(i).file; s={' ','',''};
         if ~isempty(find(f=='@',1))
            s{2}=regexprep(f,'.*/(@[^/]*/).*','$1');
         else q=''; end

         f=regexprep(f,'^.*/','');

         if i==k, s{1}='*'; end
         f2=S(i).name;
            if isequal(f2,f(1:end-2)), f2=''; else f2=['(' f2 ')']; end
            s{3}=sprintf('%s:%d',f,S(i).line);
         fprintf(1,'%4d%s  %s%-30s %s\n',i,s{:},f2);
     end
  end

  inl 1

end

