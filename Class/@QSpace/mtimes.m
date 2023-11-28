function C=mtimes(A,B,ia,ib)
% function C=mtimes(A,B [,ia,ib])
% overloading * operator (matrix multiply)
%
%    Options ia and ib are only relevant, if this routine
%    is explicitly called as mtimes(..). By default,
%    this assumes operators, ia=2 and ib=1 (matrix multiply).
%
% Alternative usage
% function C=mtimes(A,B,'--op')
%
%    In case that both A and B are irops (rank-3 operators),
%    by default, the routine above contracts(!) the operator index,
%    assuming it can be contracted, i.e., S'*S also performs the
%    dot product on the irop index.
%
%    In order to prevent this, the option '--op' keeps the
%    operator indices open, and if present on both A and B,
%    contracts them into a single irop index. If the operator
%    indices have different directions, a 1j symbols inserted
%    for this (therefore note that for non-trivial 1j symbols
%    such as for half-integer spins, this defines the sign
%    of the operator overall. This way, for example,
%    mtimes(F1',F2,'--op') keeps the irop index.
%
% Wb,Jul04,12 ; Nov19,23

  if isnumeric(A) || isnumeric(B)

     if nargin~=2, wbdie('invalid number of arguments'); end
     if     isnumeric(A), fac=A; C=B;
     elseif isnumeric(B), fac=B; C=A; end

     for k=1:numel(C), n=length(C(k).data);
         data=C(k).data; for i=1:n, data{i}=fac*data{i}; end
         C(k).data=data;
     end

  else
     if numel(A)~=1 || numel(B)~=1
        wbdie('invalid usage (got QSpace arrays)'); end
     if nargin<2, wbdie('too few input arguments'); end

     if isempty(A.Q) && isempty(A.data) C=QSpace(); return, end
     if isempty(B.Q) && isempty(B.data) C=QSpace(); return, end

     nargs=nargin-2; use_op=0;
     if nargs==1
        if isequal(ia,'--op'), use_op=1; nargs=0;
        elseif ischar(ia), wbdie('invalid usage (option %s)',ia); end
     end

     ra=length(A.Q); if nargs<1, ia=[]; end
     rb=length(B.Q); if nargs<2, ib=[]; end

     if ~nargs && ~use_op && ra==3 && rb==3
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
     end

     if ~ischar(ia) && ia>ra || ~ischar(ib) && ib>rb, wbdie(...
        'index out of bounds (A: %d/%d; B: %d/%d)',ia,ra,ib,rb);
     end

     if ra==3 && rb==2, p=[1 3 2]; else p=[]; end
     C=contractQS(A,ia,B,ib,p);

     rC=numel(C.Q);
     if rC==3, C.info.otype='operator'; end

     C=class(C,'QSpace');
     if rC==4 && use_op, C=fuse_op_indices(C); end
  end

end

% -------------------------------------------------------------------- %
function A=fuse_op_indices(A)

  t2=A.info.itags([2 4]); c2=zeros(1,2);
  for i=1:2
     if ~isempty(t2{i}) && t2{i}=='*', c2(i)=1;
        if length(t2{i})==1, t2{i}='';
        else t2{i}=t2{i}(1:end-1); end
     end
  end

  if isequal(t2{:}), t2=t2{1}; else t2=''; 
     wblog('WRN','option --op encountered different operator itags');
  end

  q=diff(c2);
  if q % mixed operator directions, e.g., F1'*F2)
     if q>0
         U1=getIdentity(A,2,'-0');
         A2=getIdentity(U1,2,A,4);
         A=contract(A,'24',contract(U1,'2*',A2,1),'12');
     else
         U2=getIdentity(A,4,'-0');
         A2=getIdentity(A,2,U2,2);
         A=contract(A,'24',contract(U2,'2*',A2,2),'21');
     end
  else
     A2=getIdentity(A,2,A,4);
     if c2(1)>0 
          A=contract(A,'24',A2,'12');
     else A=contract(A,'24',A2,'12*');
     end
  end
end

% -------------------------------------------------------------------- %
