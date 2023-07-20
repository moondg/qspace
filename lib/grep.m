function sx=grep(s_,varargin)
% function sx=grep(ss, pat [,opts])
%
%    find entries in cell string ss matching regexp pattern pat;
%    if ss is not a cellstring, it will be converted to one
%    separating entries by newlines, with newlines removed.
%
% Options
%
%    -i    ignore case, which uses regexpi() instead of regexp()
%    -c    always return as cell array
%    -s    return data as string
%
% Wb,Oct29,20

  getopt('init',varargin);
     iflag=getopt('-i');
     sflag=getopt('-s');
     cflag=getopt('-c');
  pat=getopt('get_last',[]);

  if isempty(pat)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  ss=s_;
  if ischar(ss), sflag=2;
     ss=textscan(ss,'%s','whitespace','\n'); ss=ss{1};
  end
  if iflag
       q=regexpi(ss,pat);
  else q=regexp (ss,pat); end

  sx=ss(find(cellfun(@isempty,q)==0));
  if sflag && ~cflag
     if numel(sx)==1, sx=sx{1};
     else
        sx=reshape(sx,1,[]);
        sx(2,:)={char(10)}; sx{end}=''; 
        sx=[sx{:}];
     end
  end

end

