function [H,HC]=fuse_data(HC)
% function [H,HC]=fuse_data(HC)
% Wb,Aug14,16

% outsourced from update_psi_2site.m // Wb,Dec20,21

% NB! certain off-diagonal blocks may drop out of HPsi
% for later Davidson iterations in HC
% => fix by adding zero blocks that permits the use of cell2mat()
% => fix missing zero-blocks given that diagonal entries must exist
% returning fixed HC as 2nd output argument

  for k=1:numel(HC.data), C=HC.data{k}; [m,n]=size(C);
     for i=1:m
     for j=i+1:n, if isempty(C{i,j})
        C{i,j}=zeros([ size(C{i,i},1), size(C{j,j},2) ]);
        if ~isempty(C{j,i}), wberr('invalid HC !?'); end
        C{j,i}=C{i,j}';
     end, end
     end
     HC.data{k}=C;
  end

  if nargout<2, return; end

  H=struct(HC);
  for i=1:numel(H.data), Hi=cell2mat(H.data{i});
     if ~isreal(Hi), ee=diag(Hi);
         e = [ norm(imag(ee)), norm(real(ee)) ];
         if e(2)>1, e=e(1)/e(2); else e=e(1); end
         if e>1E-12, wberr('got complex energy !?'); end
         if e, Hi=0.5*(Hi+Hi'); end
     end
     H.data{i}=Hi;
  end
  H=QSpace(H);

end

