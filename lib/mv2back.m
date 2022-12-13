function lh=mv2back(varargin)
% function lh=mv2back([lh,n])
%
%    move graphics handles (children) of current axis to background
%
% Examples
%
%    mv2back(n)     move n top (line) handles
%    mv2back(lh)    move given (line) handles
%    mv2back(lh,n)  move given (line) handles by n handles to the back
%
% Wb,Jul28,05

% tags: get line

  if ~nargin
     eval(['help ' mfilename]);
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  n=varargin{1};
  if isaxis(n)
     ch=get(gcf,'children');
     for i=1:length(n), ch(find(ch==n(i)))=[]; end
     set(gcf,'children',[ch;n(:)]);
     return
  end

  if isscalar(n) && isfig(n)
     ah=get(gca,'children');
     if n>=length(ah)
        wblog('ERR - # of handles larger than handles available (%d)!', length(ah))
        return
     end
     set(gca, 'children', [ah(n+1:end); ah(1:end)]);
  else
     lh=varargin{1}(:); ia=[]; ie=[];
     ph=get(lh,'Parent'); if iscell(ph), ph=unique(cat(1,ph{:})); end
     if length(ph)>1
        wblog('ERR','more than one parent for given line handles (%g) !?',...
        length(ph)); ph=ph(1);
     end
     ah=get(ph,'children'); m=numel(ah); ah_=ah;

     q=[isnumeric(ah),isnumeric(lh)];
     if xor(q(1),q(2))
        if q(1), lh=double(lh); else ah=double(ah); end
     end

     [i1,i2,I]=matchvec(ah,lh);

     if length(i2)~=length(lh)
        wblog('WRN','not all handles appear in current axes (%g)',...
        length(I.ix2)); lh(I.ix2)=[];
     end
     ah(i1)=[];

     if nargin>1, n=varargin{2}; k=min(numel(ah),n)';
          ah=[ ah(1:k); lh; ah(k+1:end) ];
     else ah=[ ah(:); lh]; end
     set(ph,'children',ah);

   % keep [xy]mark handles to the back
   % DON'T! some ymarks may be put intentionally to the foreground! // Wb,Aug01,17
   % h=[ findall(ah,'tag','xmark'); findall(ah,'tag','ymark') ];
   % if ~isempty(h) 
   %    [i1,i2,I]=matchvec(ah,h);
   %    ah(i1)=[]; ah=[ah(:); h(:)];
   %    set(ph,'children',ah);
   % end

     if ~nargout, clear lh; end
  end

end

% -------------------------------------------------------------------- %

