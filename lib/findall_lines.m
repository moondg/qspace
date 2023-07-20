function [lh,tt]=findall_lines(varargin)
% function lh=findall_lines([ah|fh])
%
%    Find all line objects in specfied axis handles (default: gca)
%    including ErrorBars etc. This also returns a single row vector,
%    and therefore is directly suitable for ranges in for loops.
%
% Wb,Oct16,20

  if nargin>1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  elseif nargin, ah=varargin{1};
  else ah=gca; end

  tt=get(get(gca,'Children'),'Type');
  if iscell(tt), tt=unique(tt); else tt={tt}; end
  nt=numel(tt);

  hh={};
  for i=numel(tt):-1:1
     if ~isempty(regexpi(tt{i},'line|errorbar'))
        hh{end+1}=findall(gca,'Type',tt{i})';
     end
  end

  lh=[hh{:}];

end

