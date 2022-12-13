function G=ctranspose(G)
% overloading the ' operator
% Wb,Nov17,11

% G=builtin('ctranspose',G); % don't!

  for k=1:numel(G), Gk=G(k); m=length(Gk.istr);
    i=regexp(Gk.istr(1:end-1),'[ .*+-]');
    if isempty(i) || i(1)==m
       if Gk.istr(end)==''''
            Gk.istr=Gk.istr(1:end-1);
       else Gk.istr=[ Gk.istr '''' ]; end
    elseif m>4 && isequal(Gk.istr(1),'[') && isequal(Gk.istr(end-1:end),']''')
         Gk.istr=Gk.istr(2:end-2);
    else Gk.istr=[ '[' Gk.istr ']''']; end

    Gk.op=Gk.op';
    Gk.hc=Gk.hc';

    if Gk.type=='+', Gk.type='-';
    elseif Gk.type=='-', Gk.type='+'; end

    G(k)=Gk;
  end

end

