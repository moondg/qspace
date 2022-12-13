function t=getqdir(varargin)
% function t=getqdir(A [,'-s'])
%
%    qdir in numeric format (+- => [+1,-1])
%    NB! input arguments are directly handed over to getitags.m
%
% Options
%
%    '-s'  return qdir as string (default: numeric; see above)
%
% Wb,Mar30,15

  if nargin>1 && ischar(varargin{end}) && isequal(varargin{end},'-s')
       sflag=1; l=nargin-1;
  else sflag=0; l=nargin;
  end

  t=getitags(varargin{1:l});

  if ~iscell(t), t={t}; end
  for i=1:numel(t)
     if isempty(regexp(t{i},'\*$')), t{i}=+1; else t{i}=-1; end
  end

  t=[t{:}];
  if sflag
     t(t>0)='+'; t(t<0)='-'; t=char(t);
  end

end

