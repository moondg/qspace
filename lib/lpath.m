function lpath()
% function lpath()
% Wb,Jan19,19

  pp=strread(path,'%s','delimiter',':');
  idx=1:numel(pp);

  R=getenv('MATLAB_ROOT');
  if ~isempty(R)
     for p = {'toolbox','examples','help','polyspace'}; p=[R '/' p{1}];
        i=find_string(pp,p);
        if ~isempty(i)
           pp{i(1)}=sprintf('%s/... (%g occurances)',p,numel(i));
           pp(i(2:end))=[]; idx(i(2:end))=[];
        end
     end
  end

  H=getenv('HOME'); np=numel(pp);

  fprintf(1,'\n');
  for i=1:np
     p=regexprep(pp{i},H,'~');
     if p(1)~='~', p=['   ' p]; end
     if i>1 && i<np && all(diff(idx(i-1:i+1))==1)
          fprintf(1,'     :  %s\n',p);
     else fprintf(1,'  %4g  %s\n',idx(i),p);
     end
  end
  fprintf(1,'\n');

end

function i=find_string(pp,p)
   I=regexp(pp,p);
   for i=1:numel(I), I{i}=~isempty(I{i}); end
   i=find([I{:}]);
end

