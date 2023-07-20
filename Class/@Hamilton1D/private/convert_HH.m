function HS=convert_HH(HH)
% function HS=convert_HH(HH)
%
%    Convert numerical HH Hamiltonian encoding
%    to more flexible structure setup HS.
%
% See also convert_HS.m
% Wb,Nov21,20

% adapted from map_HH.m // Wb,Nov21,20
% outsourced from setup_mpo_full.m // Wb,Mar08,22

  if ~isnumeric(HH) || size(HH,2)~=5, wbdie('invalid usage'); end

  nH=size(HH,1);
  wblog(' * ','converting HH -> HS (%g entries)',nH); 

  for i=nH:-1:1
     if diff(HH(i,[1 3]))
          HS(i)=struct('iop',reshape(HH(i,1:4),[2 2])','hcpl',HH(i,5));
     else HS(i)=struct('iop',HH(i,1:2),'hcpl',HH(i,5));
     end
  end

end

