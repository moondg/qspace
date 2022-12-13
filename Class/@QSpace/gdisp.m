function gdisp(A)
% function gdisp(A)
%
%    graphical information regarding given QSpace A.
%
% Wb,Feb09,11

  if nargin~=1 || numel(A)~=1
     fprintf(1,'\n'); eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  r=numel(A.Q); n=numel(A.data); ss=ones(n,r);
  for i=1:numel(A.data)
     s=size(A.data{i}); s(end+1:r)=1;
     ss(i,:)=s;
  end

  if isempty(A.info), nsym=0; na={1,1};
  else
     nsym=size(A.info.cgs,2);

     q=strread(A.info.qtype,'%s','whitespace',','); l=0; 
     if numel(q)~=nsym, A.info.qtype, nsym
     error('Wb:ERR','\n  ERR failed to interprete symmetries'); end

     cgs=A.info.cgs; sc=cell(n,nsym);

     for i=1:numel(cgs)
        s=size(cgs{i}); s(end+1:r)=1;
        sc{i}=s;
     end

     m=0;
     for i=1:nsym
        sc{1,i}=cat(1,sc{:,i});
        if norm(sc{1,i}-1)>1E-12, m=m+1; else sc{1,i}=[]; end
     end
     sc=sc(1,1:end);

     na={m,2};
  end

  nm=inputname(1); if isempty(nm), nm=''; end

ah=smaxis(na{:},'tag',mfilename); header('%M :: %s',nm); addt2fig Wb
axes(ah(1,1))

  S=mean(ss,1);
  if r>2
     [d,is]=sort(S); i2=is(end-1:end); ix=is(1:end-2);
  else i2=1:2; ix=[]; end

  plot(ss(:,i2(1)),ss(:,i2(2)),'o'); sms(2);
  label(sprintf('D_{%g}',i2(1)),sprintf('D_{%g}',i2(2)), ...
     sprintf('remaining dimensions D[%s]=[%s]',...
     vec2str(ix,'sep',','),vec2str(S(ix),'sep','x')))

if nsym==0, return, end

axes(ah(2,1))

   h=loglog(sqrt(prod(cat(2,sc{:}),2)), sqrt(prod(ss,2)),'o',...
     'Color',[1 1 1]*0.8); set(h(1),'disp','prod(all)'); hold on

   Dd=sqrt(prod(ss,2));

   Dc=cell(1,nsym); d=ones(1,nsym);
   for i=1:nsym, if isempty(sc{i}), continue; end
      Dc{i}=sqrt(prod(sc{i},2));
      d(i)=mean(Dc{i});
   end

   [d,is]=sort(d,'descend');

   for i=is, if isempty(sc{i}), continue; end
      h=loglog(Dc{i}, Dd,'o',...
        'Color',getcolor(i)); set(h(1),'disp',q{i})
   end
   xtight
   set(gca,'XTick',logspace(0,4,5),'YTick',logspace(0,4,5));

   xlabel('avg. CGC dim.');  %' \equiv prod(size(CGC{}))^{1/2}');
   ylabel('avg. data dim.'); %' \equiv prod(size(data[]))^{1/2}');
   sms(2);

   legdisp Location NorthEast

for i=1:nsym, dd=cat(1,sc{:,i}); if isempty(dd), continue; end; l=l+1;
axes(ah(l,2)); 
   plot(dd(:,i2(1)),dd(:,i2(2)),'o');
   s=mean(dd,1);
   label(sprintf('D_{%g}',i2(1)),sprintf('D_{%g}',i2(2)), ...
      sprintf('%s: remaining dimensions D[%s]=[%s]',q{i},...
      vec2str(ix,'sep',','),vec2str(s(ix),'sep','x')))
   sms(2);
end

end

