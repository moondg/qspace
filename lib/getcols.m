function ncols=getcols()
% function ncols=getcols()
% Wb,Nov11,11

% NB! after repeated calls, this may take 10 sec for each call !(*Y_!
% Wb,May02,13
  persistent nc

  if isempty(nc)
     if isdesktop>1
        [e,s]=system('command -v getcols.pl');
        if e, nc=80;
        else
           nc=system('getcols.pl -q 2>/dev/null');
           if isempty(nc) || nc<64, nc(2)=80;
              wblog('WRN','got ncols=%g -> %g',nc);
              nc=nc(2);
           end
        end
     else
        nc=80;
     end
  end

  ncols=nc;

end

