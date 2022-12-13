function args=args_to_cmd(varargin)
% function [cmd=]args_to_cmd(args,...)
%
%    All input arguments are expected to be strings that either
%    contain direct assigments, or may also contain scripts
%    which themselves set variables.
%
%    All assignments are separated by semicolon,
%    and a trailing semicolon is enforced.
%
%    If no output is requested, the command is already evaluated
%    in the caller routine.
%
% Format options
%
% <var>:=<value>
%    set default value, if variable does not exist yet
%    this simply translates into setdef('<var>',<value>)
%
% Wb,Jul05,19

  args={};
  for i=1:numel(varargin)
     if iscell(varargin{i}), args=[args, varargin{i}];
     else args{end+1}=varargin{i}; end
  end
  args=strjoin(args,'; ');

  args=regexprep(args,'\s*[;]\s*','; ');
  args=regexprep(args,'[\s,;]*$',';','emptymatch');

  if nargin==1
     s=inputname(1);
     if ~isempty(s), s=['   % ' s]; end
  else s=''; end

  fprintf(1,'>> %s%s\n',args,s);

  args=regexprep(args,'(\w+):=([^,;])+','setdef(''$1'',$2)');

  if ~nargout
         evalin('caller',args); clear args
  end

end

