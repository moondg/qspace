function q=isOrtho(A)
% function a=isOrtho(A)
%
%    This routine looks for single outgoing index. If found,
%    it returns the location of this leg within the index order,
%    or zero otherwise.
%
% Wb,Jan15,15

  q=zeros(size(A));

  for k=1:numel(A)
     t=getitags(A(k));
     i=find(~isemptycells(regexp(t,'\*$')));

     if isfield(A(k).info,'cgr')
        qdir=cat(1,A(k).info.cgr.qdir);
        if size(qdir,1)>1
           if norm(diff(qdir,[],1),'fro'), qdir
              wbdie('got qdir mismatch');
           end
           qdir=qdir(1,:);
        end
        if ~isequal(find(qdir=='-'),i), t, qdir
           wbdie('got qtag/qdir mismatch');
        end
     end
     if numel(i)==1 , q(k)=i; end
  end

end

