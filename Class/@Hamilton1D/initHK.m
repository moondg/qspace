function kc=initHK(HAM,varargin)
% function kc=initHK(HAM [,opts])
%
%    initialize Hamiltonian representation in given basis,
%    where the MPState A (stored with HAM) is assumed properly
%    normalized, with k=kc being the current site.
%
% Options
%
%   '-t'  double check orthogonality of MPState HAM->A
%
% Wb,Apr08,14

% NB! This uses the @Hamilton1D object HAM to update
% the Hamiltonian (either stored in global memory
% space or on the file system). The class object HAM
% itself is not altered, hence also not returned.
% => keep it this way! // Wb,Apr10,14

  if nargin<2, helpthis
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  getopt('init',varargin);
     tflag=getopt('-t');
  getopt('check_error');

  L=numel(HAM.mpo);
  kc=getCurrentSite(HAM);

  if kc<1 || kc>L, wberr('current site out of bounds (%g/%g)',kc,L);
  end

  X1=load_dmrg_data(HAM,1);
  for k=1:kc-1
     fprintf(1,'\r  initHK: k=%g/%g (>>) ...  \r',k,L);
     X2=load_dmrg_data(HAM,k+1);
     if tflag
        if ~isIdentityQS(contractQS(X1.AK,[1 3],X1.AK,[1 3]))
           wblog('WRN','A(%g/%g) not LR orthonormalized',k,kc); end
        [X1.AK,X2.AK]=ortho2site(X1.AK,X2.AK,'>>');
        save_dmrg_data(HAM,X1,k);
        save_dmrg_data(HAM,X2,k+1);
     end
     updateHK(HAM,k,'>>');
  end

  for k=L:-1:kc+1
     fprintf(1,'\r  initHK: k=%g/%g (<<) ...  \r',k,L);
     X2=load_dmrg_data(HAM,k-1);
     if tflag
        if ~isIdentityQS(contractQS(X2.AK,[2 3],X2.AK,[2 3]))
           wblog('WRN','A(%g/%g) not RL orthonormalized',k,kc); end
        [X2.AK,X1.AK]=ortho2site(X2.AK,X1.AK,'<<');
        save_dmrg_data(HAM,X1,k);
        save_dmrg_data(HAM,X2,k-1);
     end
     updateHK(HAM,k,'<<');
  end

end

