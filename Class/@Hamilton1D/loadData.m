function S=loadData(HAM,k,varargin)
% function S=loadData(HAM,k [,fields])
% Wb,Apr11,14

  L=numel(HAM.mpo);

  if k==0
     S=load_dmrg_info(HAM,varargin{:});
  else
     if k<1 || k>L
        wberr('index out of bounds (k=%g/%g)',k,L); 
     end
     S=load_dmrg_data(HAM,k,varargin{:});
  end

end

