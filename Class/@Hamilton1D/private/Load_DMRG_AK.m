function [Ak,AA]=Load_DMRG_AK(HAM,k,AA)
% function [Ak,AA]=Load_DMRG_AK(HAM,k,AA)
%
%    Buffer AK data around kc if sweep is still active.
%    Intended for routines like calc_Chrial.m that may operate
%    on an active job, i.e., which is still running.
%
%    Suggested initialization [like call to getCurrentSite()]:
%    [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace;
%
% Wb,Mar08,21

% outsourced from calc_Chiral.m

% NB! AA is not persistent, so it gets cleared once the caller routine exits
  persistent kc dk k1 k2

  if nargin==1 && isequal(HAM,'--status')
     Ak=add2struct('-',kc,dk,k1,k2);
     return
  end

  N=numel(HAM.mpo);

  if nargin<3
     [kc,Ic]=getCurrentSite(HAM);

     if ~Ic.active
        k1=kc; k2=kc; dk=1;
     else
        if nargin>1, dk=k; else dk=12; end
        if dk>=N, dk=N;   k1=1;    k2=N;
        elseif kc<3,      k1=1;    k2=k1+dk-1;
        elseif kc>N-2,    k2=N;    k1=k2-dk+1;
        elseif Ic.sdir>0, k1=kc-2; k2=k1+dk-1;
        else              k2=kc+2; k1=k2-dk+1;
        end
     end
     Ak=kc; AA=Ic; return
  end

  nA=numel(AA);

  if isempty(kc), wbdie('invalid usage (kc not yet initialized)');
  elseif dk && k>=k1 && k<=k2
     if k2>nA
        if k2-k1>1, wblog(' * ',...
          'preloading sites %g:%g (having kc=%g/%g)',k1,k2,kc,N); end
        for l=k1:k2
           AA(l)=load_dmrg_data(HAM,l,'AK');
        end
     end
     Ak=AA(k);
  else
     if k<=nA && ~isempty(AA(k).Q)
        wblog('WRN','ignoring non-empty AA(%d)',k); end
     Ak=load_dmrg_data(HAM,k,'AK');
  end

end

