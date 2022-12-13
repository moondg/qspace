function P=setfield(P,fld,val)
% Function P=setfield(P,fld,val)
%
%    Set field <fld> to specified value. If fld is specified
%    - by name, a new field is added if no field with such name exists yet
%    - as index, it must be within range.
%
% Wb,May17,22

  if nargin<3 || ~isnumeric(val), helpthis
     if nargin || nargout, wberr('invalid usage'), end
     clear P; return
  end

  if ischar(fld)
     for i=1:P.r
        if isequal(P.vars{i},fld)
           P.data{i}=val; P=update_size(P);
           return
        end
     end

     P=subsasgn(P,struct('type','.','subs',fld),val);

  elseif isnumeric(fld)
     if numel(fld)~=1,    error('Wb:ERR','invalid usage (fld)'); end
     if fld<1 || fld>P.r, error('Wb:ERR','field index out of bounds'); end
     P.data{fld}=val; P=update_size(P);
  end

end

