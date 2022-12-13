function lookupfont(font)
% function lookupfont(font)
%
%   print ascii table of characters for given font
%   see also listfonts.
%
% Wb,May09,13

  if nargin~=1 || ~ischar(font)
     helpthis, if nargin || nargout
     wberr('invalid usage'), end, return
  end

  lf=listfonts;
  q=regexpi(lf,font); n=numel(cat(1,q{:}));
  if ~n
     wberr('invalid font ''%s''',font); 
  elseif n>1
     wberr('ambiguous font ''%s''',font); 
  end

  for i=1:numel(q)
     if ~isempty(q{i}), font=lf{i}; break; end
  end

ah=smaxis(1,1,'tag',mfilename,'fpos',[1200 510 710 620]);
header('%M :: %s',font); addt2fig Wb
setax(ah(1,1))

   kk=33:126;

   i=[92 94 95 123 125]-(kk(1)-1);
   kk(i)=-kk(i);

   fmt=['%3d(%c): {\\color{blue}\\fontname{' font '}%c}'];
   fm2=['%3d(\\%c): {\\color{blue}\\fontname{' font '}\\%c}'];

   for k=kk
      [j,i]=ind2sub([16 16],abs(k));
      if k>0
           text(i,j,sprintf(fmt,k, k, k),'FontSize',14);
      else text(i,j,sprintf(fm2,k,-k,-k),'FontSize',14);
      end

      if k==kk(1), hold on; end
   end

   axis([2.8 9 0 17]); box on

