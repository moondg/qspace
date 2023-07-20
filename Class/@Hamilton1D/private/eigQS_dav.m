function [ee,Ig]=eigQS_dav(H,NPsi,G0);
% function [ee,Ig]=eigQS_dav(H,NPsi,G0);
%
%    do not just take lowest NPsi multiplets, since this may restribute
%    kept multiplets within the existing(!) ones, and therefore also may
%    eliminate actually targeted symmetry sectors
%
% Wb,Sep16,18

% outsourced from update_psi_2site.m // Wb,Dec20,21

   if nargin<3 || isempty(G0)
      [ee,Ig]=eigQS(H,'Nkeep',NPsi); return
   end

   [ee,Ig]=eigQS(H); n=0;
   [i1,i2,Im]=matchIndex(H.Q{1},G0.Q{1},'-s');

   if ~isempty(Im.ix1) || ~isempty(Im.ix2), H, G0
      wbdie('got symmetry sector mismatch !?');
   end
   if ~isequal(H.Q{1},Ig.EK.Q{1}) || ~isequal(H.Q{1},Ig.AK.Q{1})
      H, EK=QSpace(Ig.EK), AK=QSpace(Ig.AK)
      wbdie('got reordered symmetry sectors !?');
   end

   if ~isequal(i1,i2), G0=getsub(G0,i2); end
   for k=1:numel(i1), l=length(G0.data{k}); n=n+l;
      Ig.EK.data{k}=Ig.EK.data{k}(1:l);
      Ig.AK.data{k}=Ig.AK.data{k}(:,1:l);
   end
   Ig.NK=n;

   if n~=NPsi, wblog('WRN',...
     'total number of multiplets changed %g->%g !?',NPsi,n);
   end
end

