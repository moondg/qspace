function dispMPO(HAM,kc)
% function dispMPO(HAM,kc)
% Wb,Aug28,15

  q=HAM.mpo(kc); d=size(q.M);
  ioff=size(HAM.oez,1);

  for i3=1:d(3)
     s=cell(d(1),1); M=q.M(:,:,i3); l=0;
     for i1=1:d(1)
        s{i1}=sprintf([repmat(' %6.4g',1,d(2))],M(i1,:));
     end

     I=find(M(2:end,1));
     if ~isempty(I), q.M, wberr('invalid mpo.M(:,1,%g) !?',i3); 
     end
     J=find(M(2,[1 3:end]));
     if ~isempty(J), q.M, wberr('invalid mpo.M(2,:,%g) !?',i3); 
     end
     if i3==1 && norm(diag(M(1:2,1:2))-1), q.M
        wberr('invalid mpo.M(1:2,1:2,1) : diag()!=Id !?',i3); 
     end

     if M(1,2), l=l+1;
        if i3<1+ioff, wberr('invalid mpo.M(1,2,1:%g)',ioff); end
        s{l}=[s{l} sprintf('     %% (1,2) got local term: op(%g)',i3-ioff)];
     end

     J=find(M(1,3:end));
     for i=1:numel(J), j=J(i); l=l+1;
        [k1,i1]=get_kiop(q.iop(j+2,2,1));
        [k2,i2]=get_kiop(q.iop(j+2,2,2)); s{l}=[s{l} ...
          sprintf('     %% start of coupling (site %g-%g; op %g)',...
          k1,k2, i3-ioff)
        ];
        if k1~=kc || i1~=i3-ioff, q.iop
           wberr('inconsistent iop !?');
        end
     end

     if ioff==2 && i3==2, l=l+1;
        s{l}=[s{l} '     % operators propagated using fermionic parity op'];
     elseif l<2, l=2;
     end

     I=find(M(3:end,2));
     for j=1:numel(I), i=I(j); l=l+1;
        [k0,i0]=get_kiop(q.iop(i+2,1)); s{l}=[s{l} ...
          sprintf('     %% end of coupling (site %g-%g; op %g/%g)',...
          k0, kc, i0, i3-ioff)
        ];
     end

     [I,J]=find(M(3:end,3:end));
     for ip=1:numel(I), i=I(ip); j=J(ip); l=l+1;
        [k1,i1]=get_kiop(q.iop(i+2,1,1));
        [k2,i2]=get_kiop(q.iop(i+2,2,2));
        [k3,i3]=get_kiop(q.iop(j+2,2,1)); s{l}=[s{l} ...
          sprintf('     %% propagate coupling at (%g,%g): (site %g-%g; op %g)',...
          i+2,j+2, k1,k2,i1)
        ];
        if i3>2, wberr('invalid propagating term M(3:end,3:end,%g) !?',i3);
        elseif i1~=i3, disp(q.iop);
           wberr('inconsistent iop (%g,1) <> (%g,2) !?',i+2,j+2);
        elseif k1~=k3, l=l+1; s{l}=[s{l} sprintf(['     ' ...
           '%% NB! operator from site %g added with %g'], k1,k3)
        ];
        end
     end

     ss=strvcat(s{:});
     fprintf(1,'\nHAM.mpo(%g).M(:,:,%g) = \n\n',kc,i3);
     disp(ss);
  end

  fprintf(1,'\n');

end
