function save_dmrg_data(HAM,S,k)
% function save_dmrg_data(HAM,S,k)
% see also @Hamilton1D/private/setupStorage.m
% Wb,Apr10,14

  if nargin~=2 && ~isnumber(k), wberr('invalid usage'); end
  if k<1 || k>numel(HAM.mpo)
     wberr('invalid usage (k out of bounds: %g/%g',k,numel(HAM.mpo));
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

  if ~isempty(HAM.store), s=HAM.store;

     ss=struct; % initialize! // INIT_VARS_EVAL `i ml' 
     eval(['global ' s '; ss=fieldnames(' s ');']);

     ss=setdiff(fieldnames(S), ss);
     if ~isempty(ss), for q=ss
         eval(sprintf('%s(%g).%s=[];',s,k,q{1})); end
     end
     eval(sprintf('%s(%g)=S;',s,k));

  elseif ~isempty(HAM.mat)
     mat=sprintf('%s_%03g.mat',HAM.mat,k);
     save2(mat,'-q','-f','-struct','S');

  else wberr('invalid storage specification'); end

end

