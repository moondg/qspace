function C=sumops(CC,varargin)
% Wb,Nov16,11

   getopt('init',varargin);
      fac=getopt('fac',[]);
   varargin=getopt('get_remaining'); narg=length(varargin);

   if narg==2, ii=sub2ind(size(CC),varargin{1},varargin{2});
   elseif nargin==1, ii=varargin{1};
   else ii=1:numel(CC); end

   if isempty(fac), fac=ones(size(ii));
   elseif numel(fac)~=numel(ii), error('Wb:ERR',...
   '\n   ERR invalid weight factors (dimension mismatch)'); end

   if isempty(ii), C=CC([]); return; end

   C=CC(ii(1));

   for i=1:numel(ii), Ci=CC(ii(i));
      Ci.op = fac(i)*Ci.op;
      Ci.hc = fac(i)*Ci.hc;

      if fac(i)==+1, s='+'; elseif fac(i)==-1, s='-';
      else s=sprintf('%+g',fac(i)); end

      if i==1
         if fac(i)==1, s=''; else s=[ s(1) ' ' s(2:end)]; end
         C=Ci; C.istr=[ s C.istr ];
      else
         C.istr = [ C.istr ' ' s(1) ' ' s(2:end) Ci.istr ];
         C.op = C.op + Ci.op;
         C.hc = C.hc + Ci.hc;
      end
   end

end

