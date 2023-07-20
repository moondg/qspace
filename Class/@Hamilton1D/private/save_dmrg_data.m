function save_dmrg_data(HAM,S,k)
% function save_dmrg_data(HAM,S,k)
% see also @Hamilton1D/private/setupStorage.m
% Wb,Apr10,14

  if nargin~=2 && ~isnumber(k), wbdie('invalid usage'); end
  if k<1 || k>numel(HAM.mpo)
     wbdie('invalid usage (k out of bounds: %g/%g',k,numel(HAM.mpo));
  end

  t={''};
  if isfield(S,'OP')
     if ~isempty(S.OP) && ~isempty(S.OP(1).info), t=S.OP(1).info.itags; end
  else
     if ~isempty(S.HK) && ~isempty(S.HK(1).info), t=S.HK(1).info.itags; end
  end
  if ~isempty(regexp([t{:}],'^E')), t
     wblog('WRN','storing E[lr] data !?'); wbstop
  end

  if ~isempty(HAM.store), st=HAM.store;

     fs=struct; % initialize! // INIT_VARS_EVAL `i ml' 
     eval(['global ' st '; fs=fieldnames(' st ');']);

     fx=setdiff(fieldnames(S), fs);
     for i=1:numel(fx), eval(sprintf('%s(%g).%s=[];',st,k,fx{i})); end

     fx=setdiff(fs,fieldnames(S));
     for i=1:numel(fx), eval(sprintf('S.%s=%s(%g).%s;',fx{i},st,k,fx{i})); end

     eval(sprintf('%s(%g)=S;',st,k));

  elseif ~isempty(HAM.mat)
     mat=sprintf('%s_%03g.mat',HAM.mat,k);
     save2(mat,'-q','-f','-struct','S');

  else wbdie('invalid storage specification'); end

end

