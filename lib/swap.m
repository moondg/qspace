function [B,A]=swap(A,B)
% function [B,A]=swap(A,B)
% Wb,Feb23,21

  if nargout>1, return; end
  if nargout, wberr('invalid usage'); end

  nm=cell(1,2);
  for i=1:2, nm{i}=inputname(i);
     if isempty(nm{i}), wberr(...
       'invalid usage (failed to derive names for inputs)');
     end
  end

  assignin('caller',nm{1},B);
  assignin('caller',nm{2},A); clear B

end

