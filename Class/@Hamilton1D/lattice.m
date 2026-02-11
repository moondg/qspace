function [type,L,perBC,H]=lattice(HAM)
% function [type,L,perBC,H]=lattice(HAM)
% 
%    get lattice type; if outside the checked set of patterns below,
%    the return values are type='unknown' and perBC=nan.
% 
% Wb,Jul31,17

  if isempty(HAM) || ~isfield(HAM.info,'HH')
     wbdie('invalid usage (missing field HAM.info.HH)'); end

  L=numel(HAM.mpo);
  perBC=getfield2(HAM,'info','param','perBC',{0});
  naklt=getfield2(HAM,'info','param','naklt',{0});

  HH=HAM.info.HH;
  H=sparse(HH(:,1),HH(:,3),HH(:,end),L,L);
  if nargout>3
     i=find(HH(:,1)~=HH(:,3));
     H=H+sparse(HH(i,3),HH(i,1),HH(:,end),L,L);
  end

  for i=L-1:-1:0
     nd(i+1)=numel(find(diag(H,i)));
  end

  type='unknown';

  s1=full(diag(H,1));
  if norm(s1(2:end-1))==0 && all(nd(4:end)==0)
     type='ichain';
     if naklt, type=[type '/AKLT']; end
  elseif all(nd(3:end-1)==0)
     type='chain';
     if naklt, type=[type '/AKLT']; end
  elseif norm(s1(2:2:end))==0
     if ~perBC && all(nd(4:end)==0)
        type='ladder';
     elseif perBC && all(nd([4, 6:end])==0)
        type='iladder';
     end
  end

end

