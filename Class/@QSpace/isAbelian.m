function ia =isAbelian(A,varargin)
% function i=isAbelian(A [,opts])
% Wb,Feb20,13

  getopt('init',varargin);
     lflag=getopt('-l');
  getopt('check_error');

  ia=zeros(size(A));
  apat='\<APZ\d*';

  for k=1:numel(A), Ak=A(k);
     if ~isfield(Ak.info,'cgr') || isempty(Ak.info.cgr)
        ia(k)=1; continue;
     end

     if lflag
        q=strread(Ak.info.qtype,'%s','delimiter',',');
        if isempty(q), ia(k,1)=2; continue; end
        q=regexp(q,apat); n=numel(q); 
        for i=1:n
           switch numel(q{i})
              case 1, q{i}=1;
              case 0, q{i}=0;
              otherwise wbdie('unexpected symmetry');
           end
        end
        ia(k,1:n)=[q{:}];
     else
        q=regexprep(Ak.info.qtype ,[apat ',*'],'');
        if isempty(q), ia(k)=2; continue; end
     end
  end

  if lflag && isfield(A(1).info,'cgr') && ~isempty(A(1).info.cgr)
     q=[size(ia,2), size(A(1).info.cgr,2) ];
     if diff(q)
        wblog('WRN','mismatch between size(ia,2) and rank (%g/%g)',q);
     end
  end

end

