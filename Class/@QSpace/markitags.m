function A=markitags(A,idx)
% function A=markitags(A,idx)
%
%    mark specified indices in idx for QSpace A.
%
% Wb,Aug12,23

  for k=1:numel(A)
     if isempty(A(k).Q), continue; end
     tk=A(k).info.itags;
     if any(idx>numel(tk)), wbdie('invalid usage (index out of bounds)'); end
     for i=idx
        q=regexp(tk{i},'(''*)(\**)$','tokens');
        if isempty(q), tk{i}=[ tk{i} '''' ];
        else
           q=q{1}; q=[ length(q{1}), length(q{2}) ];
           t=tk{i}(1:end-sum(q));
           if mod(q(1)+1,2), t(end+1)=''''; end
           if mod(q(2),  2), t(end+1)='*';  end
           tk{i}=t;
        end

     end
     A(k).info.itags=tk;
  end
end

