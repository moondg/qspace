function [q,EE,qq,dz,DZ]=getQSpectra(H,ws)
% function [q,EE,qq,dz,DZ]=getQSpectra(H,ws)
%
%   get symmetry-resolved energy spectra.
%
% Input
%
%   H  diagonal scalar operator (with diagonals stored only)
%   ws specifies which symmetry to pick
%
% Output
%
%   q   unique set of labels for symmetry ws
%   EE  cell array of energies, one cell for each symmetry in q
%   qq  full set of symmetry labels (expanded to same size as EE)
%   dz  full set of multiplet degeneracies
%   DZ  full set of multiplet degeneracies (expanded to same size as EE)
%
% Wb,Apr12,12

   Q=H.Q; if numel(Q)~=2 || ~isequal(Q{:})
      wbdie('invalid (scalar) operator'); end
   if isdiag(H,'-d')~=2, wbdie('got non-diagonal input operator !?'); end

   if isreal(ws), ws=getsym(H,'-I',ws); end

   DZ=getzdim(H);
   DZ=expand2cell(DZ(:,:,1),H.data);
   QQ=expand2cell(Q{1},H.data);

   H=struct(H);

   Q=round(1E6*H.Q{1})*1E-6; n=size(Q,1);
   [q,I,D]=uniquerows(Q(:,ws.j));

   m=numel(I); EE=cell(m,1); dz=cell(m,1); qq=cell(m,1);

   for j=1:m
      [EE{j},is]=sort(cat(2,H.data{I{j}}));
      dz{j}=cat(1,DZ{I{j}}); dz{j}=dz{j}(is,:);
      qq{j}=cat(1,QQ{I{j}}); qq{j}=qq{j}(is,:);
   end

   DZ=dz;
   for i=1:numel(dz)
    % if norm(diff(dz{i}))>1E-12, keyboard, end % Wb,Jul06,18
    % e=std(dz{i},[],1);
    % if any(e>1E-12), wblog('WRN',...
    %    'got varying dz within symmetry sector !? (e=%.3g)',max(e));
    %  % NB! e.g for the same abelian charge sector,
    %  % may have several different spin sectors! // Wb,Apr20,16
    %  % keyboard
    % end
      dz{i}=mean(dz{i},1);
   end

   dz=cat(1,dz{:});

   i=find(abs(dz-round(dz))<1E-12); dz(i)=round(dz(i));

end

% -------------------------------------------------------------------- %

function x=expand2cell(x,data)

   x=mat2cell(x,ones(size(x,1),1),size(x,2));
   for i=1:numel(data)
      x{i}=repmat(x{i},numel(data{i}),1);
   end

end

% -------------------------------------------------------------------- %

