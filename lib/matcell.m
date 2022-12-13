function v=matcell(v)
% function v=matcell(v)
% convert a vector or array elementwise into a cell array 
% similar to cellstr.
% Wb,Feb26,10

% tags: mat2cell, vec2cell, veccell

  v=mat2cell(v,ones(1,size(v,1)),ones(1,size(v,2)));

end

