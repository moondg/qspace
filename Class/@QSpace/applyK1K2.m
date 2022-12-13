function [S1,S2]=applyK1K2(S,H)
% function [S1,S2]=applyK1K2(S,H)
%
%    auxilliary routine that applies energy kernels
%    K1 (principal value) and K2 (remainder for diagonal values).
%    e.g. required for static magnetic susceptibility
%    (K2 accounts for correction due to degeneracies).
%
% Wb,Jul05,12

  if isempty(S) || isempty(S.data)
     S1=S; S2=S; return
  end

  if none(numel(S.Q)==[2 3]) || numel(H.Q)~=2, error('Wb:ERR',...
    '\n   ERR invalid usage'); end
  if norm(H.Q{1}-H.Q{2})>1E-12, error('Wb:ERR',...
    '\n   ERR invalid usage (H expected a scalar operator)'); end

  [Ia,Ib,I]=matchIndex(S.Q{1},H.Q{1},'-s');
  [Ja,Jb,J]=matchIndex(S.Q{2},H.Q{2},'-s');

  if ~isempty(I.ix1) || ~isempty(J.ix1), error('Wb:ERR',...
    '\n   ERR failed to identify all Q-sectors in S'); end
  if ~isequal(Ia,Ja) || ~isequal(Ia,1:length(Ia)), error('Wb:ERR',...
    '\n   ERR invalid match'); end
  for i=1:length(H.data)
     s=size(H.data{i}); if numel(s)>2 || s(1)~=1, error('Wb:ERR',...
    '\n   ERR invalid H (expecting diagonal representation)'); end
  end

  S1=S;
  S2=S;

  for i=1:numel(S.data)
     e1=H.data{Ib(i)}; n=length(e1);
     e2=H.data{Jb(i)}; m=length(e2);

     q=repmat(e2,n,1) - repmat(e1',1,m); l=find(abs(q)<1E-12);
     q=1./q; q(l)=0;

     S1.data{i}=S1.data{i}.*q;

     q(:)=0; q(l)=1;
     S2.data{i}=S2.data{i}.*q;
  end

end

