function plotH(HAM,varargin)
% function plotH(HAM [,opts])
% Wb,Feb16,17

  H5=HAM.info.HH;

  [q,I,D]=uniquerows(H5(:,[2 4])); n=numel(I);

  if n<=1
     ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
  else
     ah=smaxis(ceil((n+1)/2),2,'tag',mfilename); header('%M'); addt2fig Wb
     for i=1:n, setax(ah(i+1))
        Hi=sparse(H5(I{i},1),H5(I{i},3),H5(I{i},end));
        mp(Hi,'cmap','-gca');
        axis equal tight
     end
  end

setax(ah(1,1))

  H=sparse(H5(:,1),H5(:,3),H5(:,end));
  mp(H,'cmap','-gca');
  axis equal tight

end

