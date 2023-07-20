function save_dmrg_info(HAM,varargin)
% function save_dmrg_info(HAM [,vars])
% see also @Hamilton1D/private/save_dmrg_data.m
% Wb,Apr14,14

  if nargin<1
     wbdie('invalid usage'); end

  if ~isempty(HAM.store), s=HAM.store;

     q=struct('HAM',HAM);
     eval(sprintf('global %s_info; %s_info=q;',s,s));

  elseif ~isempty(HAM.mat), mat=HAM.mat;

     save(sprintf('%s_info.mat',mat),'HAM');

  else wbdie('invalid storage specification'); end

  if nargin>1
     for i=1:numel(varargin)
        varargin{i}={inputname(i+1), varargin{i}};
     end
     addto_dmrg_info(HAM,varargin{:});
  end

end

