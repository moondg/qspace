function Xk=update_Qop_gauge(Xk,Xp,HAM,kdir_)
% function Xk=update_Qop_gauge(Xk,Xp,HAM [,kdir])
%
%    Update operator counting total charge within block.
%    While this is not necessary for open BC, where charge can be
%    simply computed from the particle count, Qop must be computed
%    for periodic BCs. This adds ((g^2)/2) * Qop^2 to Hamiltonian
%    (see updateHK).
%
%    The switch to this operator-based setting is via
%    having HAM.param.info.gauge.gtype = '..op'
%    (see got_gauge.m which returns i=1 for this case).
%
% Wb,Aug27,23

   p=HAM.info.param; Ig=p.gauge;
   if ~isfield(Ig,'stype'), wbdie('missing field Ig.stype'); end

   Ak=Xk.AK; l=itags2odir(Ak);
   if ~l && nargin>3, kdir=kdir_;
      if kdir>0, l=2; else l=1; end
   else
      if l==2, kdir=+1; elseif l==1, kdir=-1; 
      else wbdie('invalid Ak (odir=%d)',l); end

      if nargin>3 && kdir~=kdir_
      wbdie('inconsistent kdir %d/%d',kdir,kdir_); end
   end

   k=0; L=numel(HAM.mpo);
   regexp(getitags(Ak,3),'^s(\d+)(?@k=str2num($1);)');
   if numel(k)~=1 || k<=0, wbdie('failed to extract k-index of site'); end

   Eloc=HAM.oez(1).op;
   i=1; if isempty(regexp(HAM.ops(i).info,'ferm.*occ'))
      wbdie('failed to identify local charge operator'); end
   qop=HAM.ops(i).op;

   if Ig.stype(k)<0
      qop=qop-HAM.info.param.NC*Eloc;
   end

   qop2=qop*qop';
   if kdir>0, ic='!2*'; else ic='!1*'; end

   if kdir>0 && k==1 || kdir<0 && k==L
      Q=[ contract(Ak,ic,{Ak,qop, '-op:^s'}) ...
          contract(Ak,ic,{Ak,qop2,'-op:^s'});
      ];
      Xk.Qop=repmat(Q,max(abs(Ig.stype)),1);
   else
      iop=abs(Ig.stype(k)); nQ=size(Xp.Qop,1);
      if kdir>0 && k==L || kdir<0 && k==1
         if iop>=1 && iop<=nQ
            iop=ones(1,nQ);
         else wbdie('invalid iop=%d (having k=%d/%d)',iop,k,L); 
         end
      elseif iop>=1 && iop<=nQ && ~mod(iop,1)
         iop=(1:nQ)==iop;
      else wbdie('iQop out of bounds (%d/%d)',iop,nQ);
      end

      for i=1:nQ
         Q=contract(Ak,ic,{Xp.Qop(i,1),Ak});
         if iop(i)
            Q=Q+contract(Ak,ic,{Ak,qop,'-op:^s'}); end
         Xk.Qop(i,1)=Q;

         Q=contract(Ak,ic,{Xp.Qop(i,2),Ak});
         if iop(i)
            Q=Q+contract(Ak,ic,{Ak,qop2,'-op:^s'}) ...
            + 2*contract(Ak,ic,{Xp.Qop(i,1),{Ak,qop,'-op:^s'}});
         end
         Xk.Qop(i,2)=Q;

         if itags2odir(Q,1)<0
            Q=permute(Q,'21');
         end
         ee=eigQS(Q); e=min(ee(:,1));
         if e<-1E-8, wbdie('got Qop^2 with eigenvalue %.3g',e); end
      end
   end

 % if ~nargout, n=inputname(1);
 %    if isempty(n)
 %       wbdie('invalid usage (cannot set 1st argument in caller'); end
 %    assignin('caller',n,Xk);
 %    clear Xk
 % end

end

