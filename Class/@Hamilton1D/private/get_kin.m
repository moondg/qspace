function kin=get_kin(itags)
% function kin=get_kin(itags)
% Wb,Mar30,15

  k=get_kidx(itags);

% NB! through local state space, kout may appear multiple times
  kin=unique(k(find(k>0)));

end

