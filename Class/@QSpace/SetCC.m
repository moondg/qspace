function [A,t]=SetCC(A,idx)
% function [A,t]=SetCC(A [,idx])
%
%    Set "dagger" for specified index-labels (all if idx is not specified).
%    If any A.data{k} is complex, A.data will also be complex conjugated.
%
% Wb,Aug04,12

  if nargin<2, cflag=1; idx=1:numel(A(1).Q); else cflag=0; end
  n=numel(idx);

  if ~n
     if ~isempty(A)
        wbdie('invalid empty QSpace'); end
     t={};
     return
  end

  c='*'; 

  for k=1:numel(A)
     [i,t]=gotITags(A(k));
     if n>i || max(idx)>i
        wbdie('input index set out of bounds or not unique');
     end
     for j=1:n, i=idx(j);
        if ~isempty(t{i}) && t{i}(end)==c
             if numel(t{i})>1, t{i}=t{i}(1:end-1); else t{i}=''; end
        else t{i}=[t{i},c];
        end
     end
     A(k).info.itags=t;

     if cflag
        for i=1:numel(A(k).data)
           if ~isreal(A(k).data{i}), A(k).data{i}=conj(A(k).data{i}); end
        end

        cgr=A(k).info.cgr;
        for i=1:numel(cgr), q=cgr(i).qdir;
           if ~isempty(q)
               if ~ischar(q), q
                  wbdie('invalid cgr.qdir'); end
               ip=find(q=='+');
               in=find(q=='-'); q(ip)='-'; q(in)='+';
               cgr(i).qdir=q;
           end
        end
        A(k).info.cgr=cgr;
     end
  end

end

