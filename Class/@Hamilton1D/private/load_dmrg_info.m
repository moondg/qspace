function S=load_dmrg_info(HAM,varargin)
% function S=load_dmrg_info(HAM [,fields])
% see also @Hamilton1D/private/load_dmrg_data.m
% Wb,Apr14,14

  if nargin<1, wbdie('invalid usage'); end

  tflag=0;
  for i=1:numel(varargin)
     if isequal(varargin{i},'-t'), tflag=1; varargin(i)=[]; break; end
  end

  if using_mem(HAM), s=HAM.store;
     S=struct; % initialize! // INIT_VARS_EVAL `i ml' 

     eval(sprintf('global %s_info; S=%s_info;',s,s));

  else mat=HAM.mat;

     f=sprintf('%s_info.mat',mat);
     if ~exist(f,'file')
        if tflag, S=[]; return; end
        wbdie('got non-existing info file %s !?',f);
     end

     S=load(f);
  end

  n=numel(varargin);
  if n
     if n==1
        if ~isfield(S,varargin{:}) && tflag, S=[];
        else S=getfield(S,varargin{:}); end
     elseif n>1
         S=getfields(S,varargin{:});
     end
  end

end

