function plotEdata(H,varargin)
% function plotEdata(H [,opts])
%
%    plot eigenspectrum of given QSpace,
%    expecting scalar operator upon input!
%
% see also @QSpace/plotQSpectra.m
% Wb,Feb15,13

% adapted from rnrg_armando0.m

   getopt('init',varargin);
      qfac=getopt('qfac',[]);
      qstr=getopt('qstr','');
      tstr=getopt('tstr','');
      if getopt('-gca'), ah=gca; else ah=getopt('ah',[]); end

   if isempty(ah)
        ah=getopt('get_last',[]);
   else getopt('check_error');

   H=QSpace(H);

   if ~nargin || ~isscalarop(H)
      helpthis, if nargin || nargout
      error('Wb:ERR','invalid usage'), end, return
   end

   q=isdiag(H); istr='';
   if ~q
      [ee,I]=eigQS(H); R_=H; q=2;
      if ~isreal(ee)
         error('Wb:ERR','\n   ERR got complex eigenvalues !??'); end
      H=QSpace(I.EK);
   else
      if q==1
         for i=1:numel(H.data), H.data{i}=diag(H.data{i}).'; end
         q=2;
      elseif q~=2, error('Wb:ERR',['\n   ' ...
         'ERR invalid input QSpace (not an operator!? [%g])'],q);
      end
   end

   if ~isempty(qfac)
      if ~iscell(qfac)
         if isvector(qfac), qfac=diag(qfac); end
         qfac={ qfac, zeros(max(size(qfac)),1) };
      else
         if numel(qfac)~=2 || ~isvec(qfac{2})
            error('Wb:ERR','\n   ERR invalid usage'); end
         if size(qfac,1)==1, qfac=qfac'; end
      end

      if norm(qfac{1}-diag(diag(qfac{1})))==0
           istr=sprintf('[%s]',vec2str(diag(qfac{1}),'-f'));
      else istr=sprintf('[%s]',mat2str(qfac{1})); end

      if norm(qfac{2})
         istr=sprintf('%s + [%s]',vec2str(qfac{2},'-f'));
      end
   end

   if iscell(qstr)
      qstr=sprintf(', %s',qstr{:}); qstr=qstr(3:end);
   end

   if isempty(tstr)
      if isempty(qstr), tstr=istr;
      else tstr=qstr; end
   end

if isempty(ah)
ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
end

setax(ah(1,1))

  dd=H.data; qq=H.Q{1}; nd=numel(dd);

  for i=1:nd
     if i<8, o={'o-'}; else o={'*-'}; end

     if ~isvector(dd{i})
        error('Wb:ERR','\n   ERR expecting vector data'); end
     if ~isfinite(dd{i})
        error('Wb:ERR','\n   ERR invalid data (got nan''s etc)'); end

     if isempty(qfac)
          q=qq(i,:);
     else q=(qfac{1}*qq(i,:)'+qfac{2})'; end

     ss{i}=sprintf('[%s]',vec2str(q,'-f'));
     hh(i)=plot(dd{i},o{:},'Color',getcolor(i),'Disp',ss{i});
     hold on
  end

  for i=1:numel(hh)
     o={ss{i},'VerticalAl','middle','Color',get(hh(i),'Color')};
     h=text(numel(dd{i})+0.5,dd{i}(end),o{:});
     tp(i,:)=[ get(h,'pos'), double(h) ];
  end

  sms(4);

  ytight(1.1); dy=diff(ylim)/20;
  x=xtight; xlim([0.5 x(2)+2]); dx=diff(xlim)/10;

  [x,is]=sort(tp(:,2)); tp=tp(is,:); kmax=1; l=0;

  for it=1:size(tp,1), k=1;
     while 1, l=l+1;
        if it>1
           x=abs(tp(1:it-1,1)-tp(it,1))/dx;
           y=abs(tp(1:it-1,2)-tp(it,2))/dy;
           r=sqrt(x.^2+y.^2); if min(r)>=1, break; end
        else break; end

        h=tp(it,4); p=get(h,'Pos'); p(1)=p(1)+dx; set(h,'Pos',p);
        tp(it,1:3)=get(h,'Pos'); k=k+1;
     end

     if k>1
        if k>kmax, kmax=k; end
     end
  end

  x=xlim; x(3)=kmax*dx+3; if x(3) > x(2), xlim(x([1 3])); end

  if ~isempty(tstr), title(tstr); end

  label('multiplet index i','energy E_i');

end

