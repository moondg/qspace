function C=times(A,B,opt)
% overloading .* operator
% Options
%
%   --skip-cgc   C inherits CGC space of tensor A *as it is*
%                (yet only for the symmetry sectors shared with B).
%
% Wb,Jul04,12

  cgflag=1;

  if nargin==3 && isequal(opt,'--skip-cgc') cgflag=0;
  elseif nargin~=2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  i=[ isnumeric(A), isnumeric(B) ];
  if all(i==0)
     if     isscalar(A), A=A.data{1}; i(1)=1;
     elseif isscalar(B), B=B.data{1}; i(2)=1; end
  end

  if any(i)
     if     i(1), fac=A; C=B;
     elseif i(2), fac=B; C=A; end

     for k=1:numel(C)
         n=length(C(k).data);
         for i=1:n, C(k).data{i}=fac*C(k).data{i}; end
     end
  else
     if cgflag
        i=find([gotCGS(A), gotCGS(B)]);
        if numel(i)==1, error('Wb:ERR',['\n   ' ... 
           'ERR invalid usage (got operator with and without CGC spcaces !??)']);
        elseif  isempty(i), cgflag=0; end
     end

     QA=cat(2,A.Q{:});
     QB=cat(2,B.Q{:});

     if numel(A.Q) ~= numel(B.Q) || size(QA,2) ~= size(QB,2)
        error('Wb:ERR','\n   ERR incompatible QSpaces when using .*');
     end

     if size(uniquerows(QA),1)~=size(QA,1) ...
     || size(uniquerows(QB),1)~=size(QB,1)
        error('Wb:ERR','\n   ERR QSpaces do not have unique Q-labels');
     end

     [Ia,Ib,I]=matchIndex(QA,QB,'-s');

     C=A;
     C.data=A.data(Ia);
     for i=1:numel(C.Q), C.Q{i}=A.Q{i}(Ia,:); end

     for i=1:numel(Ia)
         C.data{i}=C.data{i} .* B.data{Ib(i)};
     end

     if cgflag
        error('Wb:ERR','\n   ERR invalid usage (got non-abelian CGC data');
     end
  end

end

