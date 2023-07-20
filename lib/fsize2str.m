function s=fsize2str(s)
% Function s=fsize2str(s)
% Wb,Apr08,08

   if nargin~=1
      eval(['help ' mfilename]);
      if nargin || nargout, wbdie('invalid usage'), end, return
   end

   if ischar(s)
      w=dir(s);
      if numel(w)~=1
         wblog('ERR','invalid usage (file not found; %g)',numel(w));
         s='???'; return
      end
      s=w.bytes;
   end

   s=num2str2(s,'--bytes');

end

