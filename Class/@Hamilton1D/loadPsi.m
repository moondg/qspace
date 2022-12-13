function A=loadPsi(HAM,kk)
% function A=loadPsi(HAM)
% Wb,Apr11,14

  L=numel(HAM.mpo); if nargin<2, kk=1:L; end

  nk=numel(kk); A=QSpace(1,nk);
  vflag=(nk>1);

  for i=1:nk, k=kk(i);
     if vflag
        fprintf(1,'\r   loading AK(%g/%g) ...   \r',k,L); end
     A(i)=load_dmrg_data(HAM,k,'AK');
  end
  if vflag, fprintf(1,'\r%60s\r',''); end

end

