function display(G,vflag)
% function display(G [,vflag])
% Wb,Nov17,11

  if nargin<2 || isempty(vflag), vflag=0;
  elseif ischar(vflag)
     if ~isempty(regexp(vflag,'^[a-zA-Z]')), vflag=0;
     elseif isequal(vflag,'-V'), vflag=3;
     elseif isequal(vflag,'-v'), vflag=2;
     elseif vflag, vflag=1; else vflag=0; end
  end

  nG=numel(G);
  if ~nG, fprintf(1,'\n   empty SymOp object\n\n'); return; end
  if nG==1
     [s,nrm2]=get_info(G,1);
     fprintf(1,'\n   SymOp object%-20s norm2=%g\n\n',s,nrm2);
     disp_op(G,vflag);
     return
  end

  s=size(G); r=numel(s); ii=cell(1,numel(s));
  isv=r==2 && sum(s~=1)==1;

  fprintf(1,'\n   SymOp object%-20s norm2\n\n',...
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

     [s,nrm2]=get_info(G(i),0);
     fprintf(1,'     %s %-22s %8g\n',is,s,nrm2);
     disp_op(G(i),vflag);
  end
  fprintf(1,'\n');

end

% -------------------------------------------------------------------- %

function [s,nrm2]=get_info(G,vflag)

   if ~isempty(G.op) || ~isempty(G.hc)
      nrm2=trace(G.op*G.op');
      s={ G.type, sprintf('%gx%g',size(G.op)), G.istr };
         if isempty(s{1}), s{1}=''; end
         if isempty(s{3}), s{3}=''''''; end
      s=sprintf('%s%s : %s',s{:});
   else
      nrm2=0;
      if vflag, s=' (empty)'; else s='[]'; end
   end

end

% -------------------------------------------------------------------- %
% Wb,Jul31,24

function disp_op(G,vflag)
   if vflag
      if isempty(G.op) || vflag<2, disp(G.op);
      else fprintf(1,'\n');
         if vflag<3
              disp(sparse(G.op));
         else disp(full(G.op)); end
      end
   end
end

% -------------------------------------------------------------------- %

