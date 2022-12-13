function ia =isAbelian(A,varargin)
% function i=isAbelian(A [,opts])
% Wb,Feb20,13

  getopt('init',varargin);
     lflag=getopt('-l');
  getopt('check_error');

  ia=zeros(size(A));

  for k=1:numel(A), Ak=A(k);
     if ~isfield(Ak.info,'cgr') || isempty(Ak.info.cgr)
        ia(k)=1; continue;
     end

     if lflag
        q=strread(Ak.info.qtype,'%s','delimiter',',');
        if isempty(q), ia(k,1)=2; continue; end
        q=regexp(q,'\<A\>'); n=numel(q); 
        for i=1:n
           switch numel(q{i})
              case 1, q{i}=1;
              case 0, q{i}=0;
              otherwise error('Wb:ERR','\n   ERR unexpected symmetry');
           end
        end
        ia(k,1:n)=cat(2,q{:});
     else
        q=regexprep(Ak.info.qtype ,'\<A\>,*','');
        if isempty(q), ia(k)=2; continue; end
     end
  end

  if isfield(A(1).info,'cgr') && ~isempty(A(1).info.cgr)
  if lflag && size(ia,2)~=size(A(1).info.cgr,2)
     wblog('WRN','mismatch between size(ia,2) and rank (%g/%g)',...
     size(ia,2),numel(A(1).info.cgr));
  end
  end

end

