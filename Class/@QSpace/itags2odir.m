function odir=itags2odir(A)
% function odir=itags2odir(A)
% Wb,May15,17

  odir=zeros(size(A));
  for k=1:numel(A), odir(k)=itags_to_odir(A.info.itags); end

end

