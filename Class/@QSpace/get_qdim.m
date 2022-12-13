function [D,dd]=get_qdim(A,k)
% function [dd,D]=get_qdim(A,dim)
%
%   get cgr dimension for given dimension
%
% Wb,Nov27,15

  cgr=A.info.cgr;

  dd=ones(size(cgr));
  for i=1:numel(cgr)
     if ~isempty(cgr(i)) && numel(cgr(i).size)>=k
        dd(i)=cgr(i).size(k);
     end
  end

  D=prod(dd,2);

end

