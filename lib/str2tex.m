function s=str2tex(s,varargin)
% function s=str2tex(s [,pat])
% Wb,Feb18,19

% outsourced from param2str.m

  s=regexprep(s,'\\*\<(Lambda|Gamma|Delta|alpha|delta|omega|sigma)\>','\\$1');
  s=strrep(s,'_','\_');
  s=regexprep(s,'\\*\<epsd\>','\\epsilon_d');

  if nargin<2, return
  elseif nargin==2 && iscell(varargin{1}), pp=varargin{1};
     if mod(numel(pp),2), wberr('invalid usage'); end
     if     size(pp,1)==1, pp=reshape(pp,2,[])';
     elseif size(pp,2)~=2, wberr('invalid usage');
     end
  elseif mod(nargin-1,2), wberr('invalid usage');
  else
     pp=reshape(varargin,2,[])';
  end

  if ~all(cellfun(@ischar,pp),'all');
     wberr('invalid usage (pairs of string pattern expected)');
  end
  for i=1:size(pp,1)
     s=regexprep(s,pp{i,1},pp{i,2});
  end
end

