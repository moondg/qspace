function s=sizeof(varargin)
% function s=sizeof(var1,var2,... [,'-s'])
%
%    Return size in memory of specified variables
%    or objects. Options '-s' returns size as string.
%
% Wb,Apr07,23

  vars={}; opts={}; l=0; ss=[]; nn={};
  nflag=0;

  for i=1:nargin, v=varargin{i};
     if ischar(v)
        if ~isempty(regexp(v,'^-')), opts{end+1}=v; continue; end
        l=l+1; nn{l}=inputname(i);
        if ~isempty(nn{l}), nflag=nflag+1;
             ss(l)=evalin('caller',['sizeof(' n ')']);
        else ss(l)=sizeof_1(v);
        end
     else
        l=l+1;
        ss(l)=sizeof_1(v);
     end
  end

  if ~isempty(opts)
     if ~isequal(opts,{'-s'})
        wbdie('invalid opts%s',sprintf(' %s',opts{:})); end
     sflag=1;
  else sflag=0; end

  if nargout && ~sflag, s=ss; return
  elseif isempty(ss)
     if nargout, if sflag, s=''; else s=[]; end, end
     return
  end

  s=cell(size(ss));
  for i=1:l, s{i}=num2str2(ss(i),'--bytes'); end

  if ~nargout
     nflag=~isempty(nn); if nflag, nn(end+1:l)={''}; end
     fprintf(1,'\n'); 
     for i=1:l
        if nflag, if isempty(nn{i}), q=''; else q=nn{i}; end
             fprintf(1,'   %-12s %s\n',q,s{i});
        else fprintf(1,'   %s\n',s{i}); end
     end
     fprintf(1,'\n'); clear s
  elseif numel(s)==1, s=s{1};
  end

end

function s=sizeof_1(v)
   i=whos('v'); s=i.bytes;
end

