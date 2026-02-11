function r=get_rank(sym)
% function r=get_rank(sym)
%
%    Get rank of given symmetry.
%
% Wb,Aug17,16

  if     regexp(sym,'SU\d+$'), N=str2num(sym(3:end));   r=N-1;
  elseif regexp(sym,'Sp\d+$'), N=str2num(sym(3:end))/2; r=N;
  elseif regexp(sym,'SO\d+$'), N=str2num(sym(3:end));
         if mod(N,2), r=(N-1)/2; else r=N/2; end
  elseif isequal(sym,'A4$'), r=1;
  else sym, wbdie('invalid symmetry');
  end

end

