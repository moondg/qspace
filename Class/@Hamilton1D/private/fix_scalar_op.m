function HK=fix_scalar_op(HK,Xk,k,dir)
% function HK=fix_scalar_op(HK,Xk,k,dir)
% 
%    add 0-sectors to H if required.
%
% outsourced from initNRG()
% Wb,Apr20,14

   if ~isnumeric(dir), dir=check_dir(dir); end

   da=getDimQS(Xk.AK); da=da(1,:);
   dh=getDimQS(HK);    dh=dh(1,:);

   if dir>0
      if any(dh~=da(2))
         if Xk.info.odir~=2
            wbdie('invalid A-tensor'); 
         end
         wblog('WRN','adding 0-blocks to HK (k=%g: %g->%g)',k,dh(1),da(2));
         Q=SetCC(QSpace(contractQS(Xk.AK,[1 3],Xk.AK,[1 3])),[1]);
         HK=HK+1E-99*Q;
      end
   else
      if any(dh~=da(1))
         if Xk.info.odir~=1
            wbdie('invalid A-tensor'); 
         end
         wblog('WRN','adding 0-blocks to HK (k=%g: %g->%g)',k,dh(1),da(1));
         Q=SetCC(QSpace(contractQS(Xk.AK,[2 3],Xk.AK,[2 3])),[1]);
         HK=HK+1E-99*Q;
      end
   end

end

% -------------------------------------------------------------------- %

