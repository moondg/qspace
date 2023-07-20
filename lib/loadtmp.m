function S=loadtmp(varargin)
% Function: S=loadtmp([opts,][usual arguments for matlab load])
%
%    load variables from temporary file
%    as data exchange from other MatLab session.
%
%    NB! if no return argument is specified,
%    the variables are loaded into the workspace.
%
% Options (first argument only)
%
%   '-<#>' or #   use temporary file with # = 1 .. 9.
%
% See also save2tmp.m
% Wb,Nov12,09

  f=[ getenv('HOME') '/Matlab/tmp.mat' ];
  q=[];

  if nargin, q=varargin{1};
     if ischar(q) && ~isempty(regexp(q,'^-[0-9]+$'))
        n=str2num(q(2:end));
        if n>9 || n~=round(n), wbdie('invalid tid (%s)',q);
        else q=n; end
     elseif isnumeric(q)
        if numel(q)~=1 || q<1 || q>9 || q~=round(q), q
           wbdie('invalid usage');
        end
     else wbdie('invalid usage'); end

     if ~isempty(q), varargin(1)=[]; 
        f=strrep(f,'tmp.mat',sprintf('tmp%g.mat',q));
     end
  end

  if nargout<1
     cmd=['load ' f ' ' sprintf(' %s', varargin{:}) ];
     evalin('caller', cmd); 
  else
     S=load(f,varargin{:});
  end

end

