function sz=str2bytes(varargin)
% function sz=str2bytes(sz1,sz2,..)
% Wb,Apr07,23

  sz=nan(1,nargin);
  for l=1:nargin
     v=varargin{l}; x=[]; i=regexpi(v,'[a-df-z]');
     if isempty(i), x=str2num(v); if ~isempty(x), sz(l)=x; end
     else
        i=i(1); n=[];
        x=str2num(v(1:i-1)); if isempty(x), continue; end
        switch upper(v(i))
           case 'K', if regexp(v(i:end),'kB?\s*$'), n=10; end
           case 'M', if regexp(v(i:end),'MB?\s*$'), n=20; end
           case 'G', if regexp(v(i:end),'GB?\s*$'), n=30; end
           case 'T', if regexp(v(i:end),'TB?\s*$'), n=40; end
           case 'B', if regexp(v(i:end),'(b|bytes?)$'), n=0; end
           otherwise wblog('TST','got unit ''%s'' !?',v(i:end));
        end
        if ~isempty(n), sz(l)=x*(2^n); end
     end
  end

end

