function varargout=syntax_hl(wtype,varargin)
% function [q,e1,em]=syntax_hl(wtype [,opts])
%
%     Manage QSpace specific syntax highlighting (color coding)
%     via escape codes in string outputs such as QSpace/info.m etc.
%
%     Color output can be altered by the environmental variable
%     QS_LOG_COLOR (default: 1 in terminal mode, 0 with matlab desktop ^1)
%     which acts as a bit pattern. Setting QS_LOG_COLOR to ..
%
%         =0  does not use any color coding
%         >1  increases color coding
%
%     ^1) WRN! escape sequences are terminal specific;
%         these are not supported by the matlab desktop
%    
% See MLIB/wblog.m
% Wb,Jun30,23

% see also / adapted from MLIB/wblog.m
  if nargout>3 || nargin>1, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  q=str2num(getenv('QS_LOG_COLOR'));
  if isempty(q)
     if isdesktop>1, q=1; else q=0; end
  elseif q<0, q=0; end

  wesc={ q,'','' };
  if ~q, varargout=wesc(1:nargout); return; end

  switch wtype
    case { 'QS:info' }
     % NB! use non-gray colors for markers only
     % reserve gray for plain conj flag
     % NB! ensure similar contrast on white vs. black terminal background
     % getcolor(1) ~ [0 114 190] // blue
     % getcolor(2) ~ [217 83 25] // red/orange
     % getcolor(5) ~ [19 172 50] // green
     % color for marked itags (')
       wesc{2}=[22 200 58];
       if q>1
          wesc{2}={ wesc{2}, [128 128 128], round(0.50*wesc{2}) };
       end
    otherwise
    if nargout<3
       if isstr(wtype),  wbdie('invalid wtype=%s',wtype); 
       else disp(wtype); wbdie('invalid wtype'); end
    end
  end

  if wesc{1}
     q=str2num(getenv('QS_LOG_COLOR'));
     if ~isempty(q) && ~q, wesc={0,'',''}; end
  end

  if ~wesc{1}
     varargout=wesc(1:nargout);
     return
  end

  wc=wesc{2}; esc_=[char(27) '['];
  if ~iscell(wc), wc={wc}; end

  for i=1:numel(wc), k=wc{i};
     if numel(k)==1
        if k<0 || k>255 || k~=round(k), wbdie('invalid color code %g',k);
        elseif k<8, s='3'; else s='38;5;'; end
        wc{i}=[ esc_ s num2str(k) 'm' ];
     elseif numel(k)==3
        k=round(k); if any(k<0 | k>255), wbdie('invalid color code %g',k); end
        wc{i}=[esc_ '38;2;' sprintf('%d;%d;%dm',k)];
     else
        disp(wesc)
        wbdie('invalid wesc{2}{%d}',i);
     end
  end

  if nargout>2, wesc{3}=[esc_ '0m']; end

  if numel(wc)==1
       wesc{2}=wc{1};
  else wesc{2}=wc; end

  varargout=wesc(1:nargout);

end

