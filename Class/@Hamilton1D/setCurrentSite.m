function HAM=setCurrentSite(HAM,k2)
% function HAM=setCurrentSite(HAM,kc])
%
%    Shift current site to specificied position.
%    by default: kc=1.
%
% Wb,Jan20,16

  if nargin<2, k2=1; end
  k1=getCurrentSite(HAM); L=length(HAM.mpo); e=0;

  if isempty(k1) || k1<1 || k1>L, e=1;
  else
     X1=load_dmrg_data(HAM,k1);
     odir=itags2odir(X1.AK);
     if odir, e=2; end
  end
  if e, k1, L
     wbdie('failed to identify current site'); 
  end

  if k2==k1
     wblog(' * ','already got current site kc=%g',k2);
     return
  end

  Nkeep = getopt(HAM.info,'sweep?','Nkeep',[]);
  isw   = getopt(HAM.info,'sweep?','isw',  []); tag=' * ';

  if isempty(isw), isw=0; tag='WRN';
     if isfield(HAM.info,'sweep')
        HAM.info.sweep.isw=isw;
     else
        rtol=1E-8;
        HAM.info.sweep=add2struct('-',isw,rtol);
     end
  end
  if isempty(Nkeep), tag='WRN';
     d=getDimQS(load_dmrg_data(HAM,fix(L/2),'AK'));
     Nkeep=max(d(1,:));
     HAM.info.sweep.Nkeep=Nkeep;
  end

  wblog(' * ','%sset current site kc=%g -> %g',iff(k2==1,'re',''),k1,k2);
  wblog( tag ,'using isw=%g, Nkeep=%g',isw,Nkeep);

  if k2>k1, dir=+1;
     for k=k1+1:k2
        X2=load_dmrg_data(HAM,k);
        [X1,X2,I]=update2Site(HAM,X1,X2,dir);
        wblog_iter(k-1,dir,L,I.svd2tr,I,isw);
        save_dmrg_data(HAM,X1,k-1); X1=X2;
     end
     X2.AK=fix_normPsi(X2.AK);
     save_dmrg_data(HAM,X2,k);
  else dir=-1;
     for k=k1-1:-1:k2
        X2=load_dmrg_data(HAM,k);
        [X2,X1,I]=update2Site(HAM,X2,X1,dir);
        wblog_iter(k+1,dir,L,I.svd2tr,I,isw);
        save_dmrg_data(HAM,X1,k+1); X1=X2;
     end
     X2.AK=fix_normPsi(X2.AK);
     save_dmrg_data(HAM,X2,k);
  end

end

% -------------------------------------------------------------------- %
function Ak=fix_normPsi(Ak)

   r=numel(Ak.Q);
   if r==3, return; end
   if r~=4 || isempty(regexp(Ak.info.itags{4},'Psi'))
      info(QSpace(Ak));
      error('Wb:ERR','unexpected current site'); 
   end

   X=contract(Ak,'!4*',Ak);
   [x,Ix]=eigQS(X);

   e=max(abs(x(:,1)-1)); n=trace(X);
   if e>1E-2  || abs(n-round(n))>1E-3 || n<0.9, q='WRN';
   elseif e>1E-4  || abs(n-round(n))>1E-6, q=' * ';
   else q=''; end
   if ~isempty(q), wblog(q,...
      'fixing normalization of states in Psi @ %.3g (nPsi=%g)',e,n);
   end

   Z=Ix.EK;
   for i=1:numel(Z.data)
      Z.data{i} = 1./sqrt(Z.data{i});
   end

   Q=contractQS(Ix.AK,2,diag(QSpace(Z)),1);
   Ak=contract(Ak,4,Q,1);

end

% -------------------------------------------------------------------- %

