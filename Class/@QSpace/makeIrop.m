function A=makeIrop(A,t_op)
% function A=makeIrop(A [,t_op])
%
%    make A an irreducible operator with a proper 3rd irop index,
%    which thus is only intended for rank-2 scalar operators.
%    Operators already in rank-3 format are checked for the
%    right format.
%
% Options
%
%    t_op specifies itag for the irop leg, only if added
%    (empty by default).
%
% See also: fixScalarOp, appendSingletons.
% Wb,Feb12,16

  if nargin<2, t_op=''; end
  for k=1:numel(A), A(k)=makeIrop_1(A(k),t_op); end

end

% -------------------------------------------------------------------- %
function A=makeIrop_1(A,t_op)

   r=numel(A.Q); if ~r, return; end
   if r==2
      if isempty(A.info.cgr)
         A.Q{3}=A.Q{1}-A.Q{2};
         A.info.itags{3}='*';
      else
         E3=getIdentityQS(A,getvac(A));
         A=QSpace(contractQS(A,2,E3,'3*'));
      end
      if ~isempty(t_op), A=setitags(A,t_op,3); end
   elseif r~=3
      wberr('invalid usage (got rank-%g op)',r);
   end

   if datasize(A,'-s')<1E4
      try
         x=contractQS(A,'12',A,'21');
      catch me
         A, wberr('invalid operator'); 
      end
   else
      [q1,D1,d1]=getQDimQS(A,1);
      [q2,D2,d2]=getQDimQS(A,2); [i1,i2,Ix]=matchIndex(q1,q2);
      if ~isempty(i1)
         if isequal(D1(i1,:),D2(i2,:)) || isequal(d1(i1,:),d2(i2,:))
            A, wberr('invalid operator A(%g)'); 
         end
      end
   end
end

% -------------------------------------------------------------------- %

