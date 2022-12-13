function [pp,j]=getfield(P,fld)
% Function pp=getfield(P,fld)
%
%    Get field <fld>, either specified by name
%    in which case this becomes equivalent to P.<fld>
%    or by index.
%
% Wb,May17,22

  if ischar(fld)
     for i=1:P.r
        if isequal(P.vars{i},fld)
           pp=P.data{i}; return
        end
     end
     error('Wb:ERR','invalid field ''%s'' for given PSet',fld); 
  elseif isnumeric(fld)
     if any(fld<1 | fld>P.r)
        error('Wb:ERR','field index out of bounds for PSet'); end
     pp=P.data(fld); if numel(pp)==1, pp=pp{1}; end
  end

end

