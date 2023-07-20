function printf(varargin)
% function printf(fmt,...)
%
%    Emulating regular printf()
%    while also accepting escape sequences, etc.
%
% Wb,Aug13,05 ; Wb,Aug02,22

  fid=1; k=1;
  if nargin && isnumber(varargin{1})
     fid=varargin{1}; k=2;
  end

  if nargin<k, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  dkt=isdesktop();

  for i=k:nargin, q=varargin{i};
     if any(q=='\'), q=regexprep(q,'\\n\>',char(10)); end
     if any(q=='\')
        if dkt>1
             q=regexprep(q,'\\e\>','\x1B');
        else q=regexprep(q,'(\\e|\x1B)\[[\d;]+m',''); end
     end
     varargin{i}=q;
  end

  fprintf(fid,varargin{k:end});

end

