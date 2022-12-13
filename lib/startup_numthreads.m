function startup_numthreads(I)
% function startup_numthreads(I)
% Wb,Jul18,19

  if ~nargin, I=mlinfo; end

  n4={I.nthreads.OMP, I.nthreads.MKL, I.nthreads.QSP, I.nthreads.slots};
     n4(find(cellfun(@isempty,n4)))={-1};
  n4=[n4{:}];

  nK=max(n4(1:2));
  nQ=n4(3);
  nS=n4(4);

  if ~isempty(nS) && nS>1
     if nK<1 && nQ>1 && nQ<nS
        x=round(nQ/nS); if x>=1, nK=x; end
     elseif nK>1 && nQ<1 && nK<nS || nK*nQ>nS
        x=round(nK/nS); if x>=1, nQ=x; end
     end
  end

  n2=[max([nS,nK,nQ]),nK];
  if any(n2>0) || isdeployed()
     setNumThreads(n2,'-q'); % do not report `changes' during startup
  end

  if nQ~=n4(3)
     if nQ>=1, s=num2str(nQ); else s=''; end
     setenv('QSP_NUM_THREADS',s);
  end

end

