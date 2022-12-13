function display(G)
% function display(G)
% Wb,Nov17,11

  nG=numel(G);
  if ~nG, fprintf(1,'\n   empty SymOp object\n\n'); return; end
  if nG==1,
     fprintf(1,'\n   SymOp object: %s\n\n', get_info(G,1));
     return
  end

  s=size(G); r=numel(s); ii=cell(1,numel(s));
  isv=r==2 && sum(s~=1)==1;

  fprintf(1,'\n   SymOp object%s\n\n',...
     [sprintf(': %g',s(1)), sprintf('-by-%g',s(2:end))]);

  fmt=repmat({',%g'},1,r); fmt(find(s>9))={',%2g'};
  if isv, fmt=[ fmt{find(s>1)}(2:end) '. ' ];
  else fmt=cat(2,fmt{:}); fmt=['(' fmt(2:end) ')']; end

  for i=1:nG,
     if isv, is=sprintf(fmt,i); 
     else
       [ii{:}]=ind2sub(s,i);
       is=sprintf(fmt,ii{:});
     end

     fprintf(1,'     %s %s\n',is, get_info(G(i),0));
  end
  fprintf(1,'\n');

end

function s=get_info(G,vflag)

   if ~isempty(G.type), s=sprintf('%s  ',G.type); else s=''; end

   if ~isempty(G.op) || ~isempty(G.hc)
      s=sprintf('%s%gx%g : ',s,size(G.op));
      if ~isempty(G.istr), s=[s G.istr]; else s=[s '''''']; end
   else
      if vflag, s=' (empty)'; else s='[]'; end
   end

end

