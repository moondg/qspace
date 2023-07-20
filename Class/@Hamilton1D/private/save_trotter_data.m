function HAM=save_trotter_data(HAM,S,k)
% function save_trotter_data(HAM,S,k)
% see also @Hamilton1D/private/setup_mpo_rotter2.m
% Wb,Jun02,19

% adapted from save_dmrg_data.m // Wb,Jun02,19

  L=numel(HAM.mpo);

  if nargin>2 && ~isnumber(k), wbdie('invalid usage'); end
  if k<1 || k>L, wbdie(...
    'invalid usage (k out of bounds: %g/%g',k,numel(HAM.mpo)); end

  mat=HAM.user.trotter.mat;
  if ~isempty(mat), mat=sprintf('%s_%03g.mat',mat,k);
     save2(mat,'-q','-f','-struct','S');
  else
     if ~isfield(HAM.user.trotter,'data') || isempty(HAM.user.trotter.data)
        HAM.user.trotter.data=repmat(emptystruct(S),1,L);
     end
     HAM.user.trotter.data(k)=S;
  end

end

