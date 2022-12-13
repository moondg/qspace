function s=emptystruct(s)
% Function s=emptystruct(s)
%
%    return structure with same set of fields
%    yet all fields initialized to [].
%
% Wb,Jan29,09

  s=fieldnames(s)'; s(2,:)={[]}; s=struct(s{:});

end

