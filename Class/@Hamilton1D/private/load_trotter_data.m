function q=load_trotter_data(HAM,k,varargin)
% function load_trotter_data(HAM,k [,vars])
% see also @Hamilton1D/private/setup_mpo_trotter2.m
% Wb,Jun02,19

  L=numel(HAM.mpo);

  if nargin<2 || ~isnumber(k), wbdie('invalid usage'); end
  if k<0 || k>L, wbdie(...
    'invalid usage (k out of bounds: %g/%g',k,numel(HAM.mpo)); end

  sflag=nargin>2;
  if k==0
     q=HAM.user.trotter.info;
  else
     mat=HAM.user.trotter.mat;
     if ~isempty(mat), mat=sprintf('%s_%03g.mat',mat,k);
        q=load(mat,varargin{:});
        sflag=0;
     else
        q=HAM.user.trotter.data(k);
     end
  end

  if sflag
     if numel(varargin)>1
          q=getfields(q,varargin{:});
     else q=getfield (q,varargin{:});
     end
  end

end

