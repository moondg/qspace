function xx=getscalar(A)
% function x=getscalar(A)
%
%    Assuming that A is the result out of full contraction,
%    this results in a scalar that can be extracted here
%    with safety checks.
%
% Wb,Aug28,08

  xx=zeros(size(A));
  for k=1:numel(A), Ak=A(k);
     if numel(Ak.data)==1 && numel(Ak.data{1})==1
        if ~isempty(Ak.Q) && size(Ak.Q{1},1)>1
        wblog('WRN','severe QSpace inconsistency !?'); end

        x=Ak.data{1};

        if ~isempty(Ak.info) && isfield(Ak.info,'cgr') && ~isempty(Ak.info.cgr)
           cgr=Ak.info.cgr; 
           for i=1:numel(cgr)
              cgr(i).cgw=norm(cgr(i).cgw);
           end
           x=x*prod(cat(1,cgr.cgw));
        end

     elseif (isempty(Ak.Q) || isempty(Ak.Q{1})) && ...
            (isempty(Ak.data) || isempty(Ak.data{1}))
        x=0;
     else
        display(A)
        wbdie('invalid scalar QSpace');
     end

     q=[abs(x), normQS(A(k))];
     if all(q>1E-20)
        e=abs(diff(q)/max(q)); if e>1E-8
           wbdie('severe inconsistency in scalar value');
        end
     end
     xx(k)=x;
  end

end

