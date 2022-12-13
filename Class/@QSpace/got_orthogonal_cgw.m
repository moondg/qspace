function [iu,e]=got_orthogonal_cgw(A,nm)
% function got_orthogonal_cgw(A)
%
%    (double) check whether A got orthogonal cgw spaces
%
% Wb,Feb19,16

  iu=1; e=0;
  [q,I,d]=uniquerows(cat(2,A.Q{:})); if all(d<2), return; end

  cgw=reshape({A.info.cgr.cgw},size(A.info.cgr));
  iu=2;

  dd=ones(size(cgw));
  for i=1:numel(cgw)
      cgw{i}=mpfr2dec(cgw{i});
      if isempty(cgw{i}), cgw{i}=1; end
      dd(i)=numel(cgw{i});
  end

  [n,m]=size(cgw);

  for i=find(sum(diff(dd,[],2).^2,2))'
      l=max(dd(i,:));
      for j=1:m, cgw{i,j}(end+1:l)=0; end
  end

  iom=find(d>1);
  for k=iom, W={};
     for i=I{k}, W{i}=mkron(cgw{i,:}); end
     W=cat(1,W{:});

     E=W*W'; e=norm(E-eye(size(E)),'fro');

     if e>1E-12, iu=0;
        if nargout<2
           if size(W,1)<10, W, E, end

           if nargin<2, nm=inputname(1); end
           if ~isempty(nm)
                wblog('WRN','%s got non-orthonormal cgw (@ e=%.3g) !?',nm,e);
           else wblog('WRN','cgw not orthonormal (@ e=%.3g) !?',e);
           end
        end
        return
     end
  end

end

