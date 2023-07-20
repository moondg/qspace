function k=get_kidx(itags)
% function k=get_kidx(itags)
% 
%    return k index specified with itags
% 
% See also QSpace/itags2int.m and QSpace/get_odir.m
% Wb,Mar30,15

  if ~nargin
     helpthis, if nargout, wbdie('invalid usage'), end
     return
  end

  if isstruct(itags) && isfield(itags,'itags')
     itags=itags.itags;
  elseif ischar(itags), itags={itags};
  elseif isa(itags,'QSpace') && numel(itags)==1
     itags=itags.info.itags;
  elseif ~iscell(itags) || isempty(itags) || ~ischar(itags{1})
     wbdie('got invalid input itags'); 
  end

  for i=1:numel(itags)
     q=str2num(regexprep(itags{i},'^[^\d]*(\d+)[^\d]*$','$1'));
     if isempty(q), itags{i}=nan;
     elseif isempty(regexp(itags{i},'\*$')), itags{i}=+q;
     else itags{i}=-q; end
  end

  k=[itags{:}];

end

