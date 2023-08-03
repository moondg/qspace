function setupStorage(HAM)
% function setupStorage(HAM)
% see also @Hamilton1D/private/get_dmrg_data.m
% Wb,Apr10,14

  save_dmrg_info(HAM);

  L=numel(HAM.mpo); q=QSpace;

  I.AK=q;
  I.HK=q;

  if using_full_MPO(HAM)
       o={'hconj',use_Hconj(HAM)};
  else o={}; I.OP=q;
  end

  I.Psi=q;

  I.info=struct(...
     'itags','','time',[],'cpu',[],'isw',[],o{:},'ek',{{}}, ...
     'dsw',[],'Eg',[],'Rho',QSpace(0,1));
  I=repmat(I,1,L);

  if ~isempty(HAM.store), s=HAM.store;
     eval(sprintf('global %s; %s=I;',s,s));

  elseif ~isempty(HAM.mat), mat=HAM.mat; nfx=0;

     for k=1:L, qk=I(k);
        f=sprintf('%s_%03g.mat',mat,k);
        if exist(f,'file'), nfx=nfx+1; end
        save(f,'-struct','qk');
     end

     while 1, k=k+1;
        f=sprintf('%s_%03g.mat',mat,k);
        if exist(f,'file')
           system(['/bin/rm ' f]); nfx=nfx+1;
        else break; end
     end

     if nfx
          wblog('I/O','setting up %g files (%g removed)',L,nfx); 
     else wblog('I/O','setting up %g files',L); end
     wblog('mat','%s_*.mat',repHome(mat));

  else wbdie('invalid storage specification'); end

end 

