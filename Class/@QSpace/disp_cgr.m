function disp_cgr(A)
% function disp_cgr(A)
% see also plot_cgr()
% Wb,Jan13,12

  if nargin~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(A.Q), fprintf(1,'\n   (empty QSpace)\n'); return, end
  if ~gotCGS(A), fprintf(1,'\n   No CG Spaces (all abelian)\n'); return, end

  cgs=A.info.cgs; [n,m]=size(cgs); s=cell(1,m);
  fmt=repmat(' | %-12s',1,m); fmt=[ '%6d. ' fmt(3:end) '\n'];

  if isfield(cgs,'S')
     fprintf(1,'\n  cgs structure array:\n');
     for i=1:n
        for j=1:m, S=cgs(i,j).S;
           if ~isempty(S), s{j}=sprintf(' %3d',S);
           elseif numel(cgs(i,j).data)~=1
              s{j}=sprintf('[%dx%d]',size(cgs(i,j).data));
           else
              s{j}=sprintf('   %.4g',full(cgs(i,j).data));
              if isempty(find(s{j}=='.',1)), s{j}=[s{j} '.']; end
           end
        end
        fprintf(1,fmt,i,s{:});
     end
  else
     fprintf(1,'\n  cgs cell array:\n');
     for i=1:n
        for j=1:m, s{j}=sprintf(' %3d',size(cgs{i,j})); end
        fprintf(1,fmt,i,s{:});
     end
  end
  fprintf(1,'\n');

end

