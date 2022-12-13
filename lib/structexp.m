function structexp(S,varargin)
% Function: structexp(S)
%
%    Expand given structure by setting every field
%    as variable in the calling workspace
%
% Wb,May11,07

  getopt('init',varargin);
     xflag=getopt('eval'         );
     xpatv=getopt('exceptions',{});
     xpatf=getopt('-x',{});
  varargin=getopt('get_remaining');

  if isempty(S), return; end

  if ~iscell(xpatv), xpatv={xpatv}; end

  if ~isempty(varargin)
     ff=varargin;
     for i=1:numel(ff)
        if ~ischar(ff{i}), ff{i}, wberr('invalid field');
        elseif ~isfield(S,ff{i}), wberr('invalid field ''%s''',ff{i}); end
     end
  else ff=fieldnames(S); end

  if ~isempty(xpatf)
     n=numel(ff); mark=zeros(1,n);
     for i=1:n
        if ~isempty(regexp(ff{i},xpatf))
        mark(i)=1; end
     end
     i=find(mark); if ~isempty(i), ff(i)=[]; end
  end

  if ~xflag
     for k=1:length(ff)
     assignin('caller',ff{k},getfield(S,ff{k})); end
     return
  else
     for k=1:length(ff)
        x=getfield(S,ff{k});
        if ~ischar(x), i=1; else i=0;
           for l=1:length(xpatv)
              if ~isempty(findstr(x,xpatv{l}))
              i=1; break; end
           end
        end

        if i
           assignin('caller',ff{k},getfield(S,ff{k}));
        else
           try, evalin('caller',[ ff{k} '=' x ';']);
           catch i
              fprintf(2,...
                '\n%s\n%s - failed to evaluate expression\n\n  `%s=%s;''\n\n',...
                i.message,lineno,ff{k},x ...
              );
              assignin('caller',ff{k},x);
           end
        end
     end
  end

end
