function A=uplus(A)
% overloading "+A" operator
% Wb,Apr04,18

  if numel(A)~=1
     error('Wb:ERR','\n   ERR invalid usage');
  end

  as=regexprep(A.istr,'^[^=]+=\s+','');
  i=regexp(as,'[^\s+-].*');
  if i>1
     q=regexprep(as(1:i-1),'[\s+]','');
     if mod(length(q),2)
          as=['- ' as(i:end)];
     else as=['+ ' as(i:end)]; end
  end

  A.istr=as;

end

