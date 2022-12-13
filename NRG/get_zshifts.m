function zz=get_zshifts(Nz)
% function zz=get_zshifts(Nz)
% Wb,May02,13

  zz=linspace(0,1,Nz+1); zz(end)=[];
  if Nz<=2, return; end

  n=log(Nz)/log(2);
  if abs(n-round(n))<1E-12, n=round(n);
     zz=reshape(permute(reshape(zz,repmat(2,1,n)),n:-1:1),[],1);
  elseif mod(Nz,2)==0
     q=Nz; n=0; while 1
       q=q/2; if abs(q-round(q))<1E-12, n=n+1; else q=2*q; break; end
     end
     s=[q, repmat(2,1,n)]; p=(n+1):-1:1;
     zz=reshape(permute(reshape(zz,s),p),[],1);
  else
     wblog('WRN','even number of z-shifts preferred (got Nz=%g)',Nz);
  end

end

