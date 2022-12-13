function display(P)
% function display(P)
% Wb,Jul31,08

   if P.n==0
      fprintf(1,'\n   Empty parameter set.\n\n');
      return
   end

   fprintf(1,'\n   parameter length  data (range)\n\n');
   for i=1:P.r
      q=P.data{i}; n=numel(q);
      if n<8,
           s=['[ ' mat2str2(P.data{i},'fmt','%8.5g') ' ]'];
      else s=sprintf('[ %8.5g .. %.5g  ]',q(1),q(end));
      end
      fprintf(1,'  %9s %5g    %s\n', P.vars{i},P.s(i), s);
   end
   fprintf(1,'\n   Total of %g combinations.\n\n',P.n);

end

