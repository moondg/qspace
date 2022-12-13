function A=sqrt(A)
% function B=sqrt(A)
%
%    Compute sqrt(A) for diagonal QSpace(s) A
%    (error is issued if A is not diagonal)
%
% Wb,Mar11,19

  neg=0; dmin=0;
  for k=1:numel(A)
     q=isdiag(A(k)); if ~q, error('Wb:ERR',...
       '\n   ERR invalid usage (diagonal QSpace required)'); end
     if q>1
        for i=1:numel(A(k).data)
           dd=A(k).data{i}; A(k).data{i}=sqrt(dd);
           [neg,dmin]=check_data(dd,neg,dmin);
        end
     else
        for i=1:numel(A(k).data)
           dd=diag(A(k).data{i}); A(k).data{i}=diag(sqrt(dd));
           [neg,dmin]=check_data(dd,neg,dmin);
        end
     end
  end

  if neg==2,  wblog('WRN','complex data encountered');
  elseif neg, wblog('WRN','negative data encountered (%.3g)',dmin);
  end

end

% -------------------------------------------------------------------- %

function [neg,dmin]=check_data(dd,neg,dmin)
  if ~neg
     if ~isreal(dd), neg=2;
     else i=find(dd<0);
        if ~isempty(i), neg=1;
           dmin=min(dmin,min(dd));
        end
     end
   end
end

% -------------------------------------------------------------------- %

