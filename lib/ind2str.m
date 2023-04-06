function s=ind2str(idx)
% function s=ind2str(idx)
%
%    Index to string;
%    e.g. [1, 4 5 6 7 8, 10] results in '[1, 4:8, 10]'
%
% Wb,Sep06,15

  if ~isvector(idx),
     wbdie('invalid usage (vector expected)'); 
  end

  n=numel(idx);
  if n<3
     if n, 
          s=sprintf(', %g',idx); s=['[ ' s(3:end) ' ]'];
     else s='[]'; end
     return
  end

  if any(diff(idx)<0)
       sflag=1; idx=sort(idx);
  else sflag=0;
  end

  idx=reshape(idx,1,[]);

  d2=diff([2*idx(1)-idx(2), idx, 2*idx(end)-idx(end-1)],2);

  ig=[find(d2), n+1];

  i=find(ig>1 & ig<n);
  j=find(~d2(ig(i)-1) & d2(ig(i)+1)); ig(j)=[];

  ng=numel(ig); s=cell(1,ng); i1=1;

  for k=1:ng, i2=ig(k);
     i=i1:i2-1;
     if i2-i1<3
        s{k}=sprintf(', %g',idx(i)); s{k}=s{k}(3:end);
     else
        di=diff(idx(i));
        if norm(diff(di)), di, wblog('ERR','!?');
        else di=di(1); end
        if di==1
             s{k}=sprintf('%g:%g',idx(i(1)),idx(i(end)));
        else s{k}=sprintf('%g:%g:%g',idx(i(1)),di,idx(i(end)));
        end
     end
     i1=i2;
  end

  if numel(s)==1, s=s{1};
  else
     s(2,:)={', '}; s{end}=''; s=[s{:}];
     s=['[ ' s ' ]' ];
  end

end

