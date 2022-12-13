function wbdie(varargin)
% function wbdie(varargin)
%
%    standardize error format
%
% see also wberr.m
% Wb,Jan19,20

  if nargin, msg=sprintf(varargin{:}); else msg=''; end

  msg=[ regexprep([10 msg],[' *' 10 ' *'],[10 '   ERR ']), 10];

  S=dbstack;
  if numel(S)>1, dispstack(S(2:end)); end

  S=struct('message',msg,'identifier','Wb:ERR', 'stack',S(min(2,end)));
  error(S);

end

