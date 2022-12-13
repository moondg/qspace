function kk=itag2iter(A)
% function k=itag2iter(A)
%
%    extract site position to k from A.info.itags
%    with expected itag format '[non-number]##[non-number]'
%
% Wb,Aug01,17

% tags: itags2iter, itag2pos, itags2pos

  kk=zeros(size(A)); n=numel(A);

  for i=1:n
     it=A(i).info.itags; r=numel(it); e=0;
     x=regexprep(it,'^[^\d]*(\d+)[^\d]*$','$1');
     for j=1:r
        x{j}=str2num(x{j});
        if isempty(x{j}), e=j; break; end
     end

     if ~e || e>2, x=[x{:}];
        if r==2
           if diff(x)==0, kk(i)=x(1); else e=2; end
        elseif r==3
           if diff(x(1:2))==0, e=-e; kk(i)=x(1);
           elseif ~e && (x(3)==x(1) || x(3)==x(2)), kk(i)=x(3);
           else e=3; end
        else e=4;
        end
     end

     if e>0, kk(i)=nan; disp(it)
        wblog('WRN','unexpected itags');
     end
  end

end

