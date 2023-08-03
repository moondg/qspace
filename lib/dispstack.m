function dispstack(S,k)
% Function dispstack(S [,k])
% Wb,Aug24,07

  if ~nargin
     [S,k]=dbstack('-completenames'); S=S(2:end);
  elseif nargin==1 && isnumber(S), i=S;
     [S,k]=dbstack('-completenames'); S=S(2+i:end);
     wblog('TST','k=%d, i=%d',k,i); 
  elseif nargin<2, k=-1;
  elseif nargin>2
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if isempty(S)
     fprintf(1,'\n   (empty dbstack)\n\n');
     return
  end

  if isobject(S), S=obj2struct(S); end

  if ~isfield(S,'file')
     if isfield(S,'stack')
        if isfield(S,'message'),    fprintf(1,'\n  Message: %s\n', S.message); end
        if isfield(S,'identifier'), fprintf(1,'       id: %s\n',S.identifier); end
        S=S.stack;
     else wbdie('invalid usage'); end
  end

% indicator for current frame // Wb,Apr24,13
% (see bt.m)
% if ~isfield(S,'curr'), S(1).curr=0; end
% for i=1:numel(S)
%     if isempty(S(i).curr), S(i).curr=0; end
% end

  n=numel(S); vflag=0; inl 1

  use_col=wblog('--hl-check');
  if use_col
       e1=[ char(27) '[38;5;8m' ]; em=[ char(27) '[0m'];
  else e1=''; em='';
  end

  if vflag
     fprintf(1,'\nDBSTACK\n\n');
     fmt=[ '%s' e1 ':%d' em];

     for i=n:-1:1
         s=sprintf(fmt, S(i).file, S(i).line);
         fprintf(1,'%4d  %-40s %s\n', n-i+1,s,S(i).name);
     end
  else
     fmt='%s:%d';
     for i=n:-1:1
         f=S(i).file; s={' ','',''};
         if any(f=='@')
            s{2}=regexprep(f,'.*/(@.*\/).*','$1');
         end

         if i==k, s{1}='*'; end

         f=regexprep(f,'^.*/','');
         s{3}=sprintf(fmt,f,S(i).line);

         f2=S(i).name;
         if isequal(f2,f(1:end-2)), f2='';
         else f2=[e1 '(' f2 ')' em];
            l=length(s{2})+length(s{3});
            if l<40, s{3}=[s{3}, repmat(' ',1,40-l)]; end
         end
         if use_col
            if ~isempty(s{2}), s{2}=[e1 s{2} em]; end
            s{3}=regexprep(s{3},'(\.[a-z]:\d+)\>',[e1 '$1' em]);
         end
         fprintf(1,'%4d%s  %s%s %s\n',i,s{:},f2);
     end
  end
  fprintf(1,'\n');

end

