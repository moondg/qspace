function dstr=get_dim_str(dd)
% function dstr=get_dim_str(dd)
%
%    Print multiplet (state space) dimension where
%    second row (if present) specifies full state space dimension
%    thus printing D* (D); otherwise, assuming D*=D, this just prints D.
%
% Wb,Feb14,21

   s=size(dd);
   if numel(s)>2 || s(1)>2, wberr('invalid usage'); end
   if any(diff(dd,[],1))<0, wberr('invalid usage'); end

   if s(1)==1
        for i=s(2):-1:1, dstr{i}=sprintf('%g',dd(i)); end
   else for i=s(2):-1:1, dstr{i}=sprintf('%g (%g)',dd(1,i),dd(2,i)); end
   end

   if numel(dstr)==1, dstr=dstr{1}; end

end

