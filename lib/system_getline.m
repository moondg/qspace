function [s,i]=system_getline(cmd)
% function s=system_getline(cmd)
% Wb,Aug21,12

% [s,i]=perl('~/bin/ml_system.pl',cmd);
% if regexp(s,'shell.*error')
%    i=i+10;
% end

  [i,s]=system(cmd);

% NB! system (as well as evalc('! cmd') capture intermediate
% input from the MatLab prompt @#(*&! (eg. by typing command
% while previous command is still running => skip!
% affects: repHome, pwd, cd, ...
% Wb,Aug21,12

  if 1
     l=find(s==10,1);
     if ~isempty(l), l=[0,l,length(s)+1];
        for j=2:numel(l), sj=s(l(j-1)+1:l(j)-1);
           if regexp(sj,'shell-init.*error')
              if isdesktop>1, e=char(27);
                   wesc_={ [e '[35;1m'], [e '[0m'] };
              else wesc_={ '','' }; end
              fprintf(1,[
              '\n   WRN ' mfilename '() system command returned' ...
              '\n   ' wesc_{1} sj wesc_{2} '\n'])
              continue
           end
           break
        end
        s=sj;
     end
  else
     l=max(find(s(1:end-1)==10));
     if ~isempty(l), s=s(l+1:end); end
  end

end

