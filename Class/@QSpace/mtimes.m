function C=mtimes(A,B,ia,ib)
% function C=mtimes(A,B,ia,ib)
% overloading * operator (matrix multiply)

  if isnumeric(A) || isnumeric(B)

     if nargin~=2, wbdie('invalid number of arguments'); end
     if     isnumeric(A), fac=A; C=B;
     elseif isnumeric(B), fac=B; C=A; end

     for k=1:numel(C), n=length(C(k).data);
         data=C(k).data; for i=1:n, data{i}=fac*data{i}; end
         C(k).data=data;
     end

  else
     if numel(A)~=1 || numel(B)~=1, wbdie('invalid usage'); end

     if nargin<2, wbdie('invalid number of arguments'); end
     if isempty(A.Q) && isempty(A.data) C=QSpace(); return, end
     if isempty(B.Q) && isempty(B.data) C=QSpace(); return, end

     ra=length(A.Q); if nargin<3, ia=[]; end
     rb=length(B.Q); if nargin<4, ib=[]; end

     if nargin<3 && ra==3 && rb==3
        qa=uniquerows(A.Q{3}); qb=uniquerows(B.Q{3});
        if size(qa,1)==1 && size(qb,1)==1
           if ~isequal(qa,qb)
              wbdie('invalid usage (operator mismatch)'); end

           ta=isempty(regexp(A.info.itags{3},'\*$'));
           tb=isempty(regexp(B.info.itags{3},'\*$'));
           if ta~=tb, ia='23';
           else ia='13*'; wblog('WRN',...
               'using dagger with first operator in mprod (*)');
           end
           ib='13';
        end
     end

     if isempty(ia)
        if isstruct(A), A=QSpace(A); end
        if ra<2 || ra>3 || ~isop(A)
           wbdie('invalid usage (A must be operator)'); end
        ia=2;
     end
     if isempty(ib)
        if isstruct(B), B=QSpace(B); end
        if rb<2 || rb>3 || ~isop(B)
           wbdie('invalid usage (B must be operator)'); end
        ib=1;

        if ra>2 && rb>2
           wbdie('invalid usage (got two rank-3 IROPs)');
        end
     end

     if ~ischar(ia) && ia>ra || ~ischar(ib) && ib>rb
        wbdie('index out of bounds (A: %d/%d; B: %d/%d)',ia,ra,ib,rb);
     end

     if ra==3 && rb==2
          p=[ 1 3 2 ];
     else p=[]; end

     C=contractQS(A,ia,B,ib,p);
     if numel(C.Q)==3, C.info.otype='operator'; end

     C=class(C,'QSpace');

  end

end

