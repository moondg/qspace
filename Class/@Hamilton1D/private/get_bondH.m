function [J,hconj]=get_bondH(HAM,k1)
% function [J,hconj]=get_bondH(HAM,k1)
%
%    get bond OPs for bond [ k1 -> (k1+1) ], as required when
%    calculating H|Psi> where H = HL + HR + sum (coupling terms) with
%    * HL = X1.XK
%    * HR = X2.XK
%    * coupling terms: check mpo{k1} => X1.OP(J) * X2.OP(J)
%    regarding Davidson, see also $MEX/david2SiteQS.cc
%
% Wb,Jan19,15

% outsourced from sweepPsi.m // Wb,Jan20,16

 % look for coupling terms out of k1 to the right (k1+1)
   M=HAM.mpo(k1).M(:,3:end,:);
   [I,J,K]=ind2sub(size(M),find(M));

   ko=HAM.mpo(k1).iop(:,2); hconj=repmat(-1,size(ko));

   for i=1:numel(ko), x=ko(i);
      if x
         k=floor(x); io=(x-k)*2^10;
         w=[HAM.mpo(k).stype];
         hconj(i)=HAM.ops(io,w).hconj;
      end
   end

   hconj=hconj(J+2);

   if any(hconj<0), HAM.mpo(k1).iop, J, hconj
      wberr('failed to determine hconj from iop !?');
   end

   if numel(J)>1
      [J,i]=unique(J); hconj=hconj(i);
   end

end

