function [D,sym,err]=get_dir(sym,varargin)
% function [D,sym]=get_dir(sym [[,'-c',] subdir])
% Options
%
%    '-c'  keep cell structure
%
% Wb,Jan30,15

% sym=regexprep(sym,'(S[pU])\(*(\d+)\)','$1$2\n');
  sym=check_sym(sym);

  D=getenv('RC_STORE');
  D=strread(D,'%s','delimiter',':');

  nD=numel(D);

  if nD==1 && ~exist([D{1} '/' sym '/CStore']);
     getRC(sym,'--ping');
  end

  ix=[]; DX={}; err=0;
  for i=1:nD
     D{i}=[ D{i} '/' sym ];
     if ~exist(D{i},'dir')
        if i==1
           DX{end+1}=D{i};
        end
        ix(end+1)=i;
     end
  end

  if ~isempty(ix), D(ix)=[];
     if isempty(D)
        if nargout>2
           err=-1;
        else
           wbdie('invalid RC_STORE (symmetry %s not yet setup)',sym);
        end
     end
  else
     err=numel(DX);
     for i=1:numel(DX)
        wblog('WRN','non-existing %s',repHome(DX{i}));
     end
  end

  if nargin>1 && isequal(varargin{1},'-c')
       cflag=1; varargin(1)=[];
  else cflag=0; end

  nargs=numel(varargin);
  if ~nargs
     if ~cflag && numel(D)==1, D=D{1}; end
     return
  end

  for i=1:nargs, if ~ischar(varargin{i})
     wbdie('invalid usage'); end
  end

  Dx=sprintf('/%s',varargin{:}); ix=[];

  for i=1:numel(D)
     D{i}=[D{i},Dx];
     if ~exist(D{i},'dir')
        ix(end+1)=i;
     end
  end

  if ~isempty(ix)
     if numel(ix)==numel(D)
        wbdie('invalid RC_STORE (./%s does not exist)',Dx(2:end)); end
     D(ix)=[];
  end

  if ~cflag && numel(D)==1, D=D{1}; end

end

