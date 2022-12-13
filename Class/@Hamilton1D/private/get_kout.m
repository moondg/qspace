function kout=get_kout(itags)
% function kout=get_kout(itags)
% Wb,Mar30,15

  k=get_kidx(itags);

% NB! through local state space, kout may appear multiple times
  kout=-unique(k(find(k<0)));

end

