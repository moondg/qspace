function ops=check_stype_transpose(ops,stype,istr)
% function ops=check_stype_transpose(ops,stype,istr)
%
%    NB! Columns in HAM.ops or HAM.oez specify site type.
%    For a uniform system with one site type (e.g. empty stype=[])
%    HAM.ops may have been (accidentally) written as row vector
%    by defining HAM.ops(i). If so, the routine here applies a
%    transpose together with issuing a warning.
%
% Wb,Mar11,26

  q=stype(:);
  ntype=max([max(q), numel(unique(q)), 1]);

  s=size(ops); n=prod(s);
  if numel(s)~=2
     s=sprintf('x%s',s);
     wbdie('invalid usage (size(%s) = %s)',istr,s(2:end));
  end
  if s(2)>ntype
     if s(1)>ntype
        wbdie('invalid usage (size(%s) = %s)',istr,s(2:end));
     end
     if s(1)==1
          istr=['switching ' istr ' to column vector'];
     else istr=['switching ' istr ' to transpose'];
     end
     wblog('WRN',[istr '\n(hint: columns specify site type)']);
     ops=ops.';
  end

end

