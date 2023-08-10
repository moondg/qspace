function wbdie(varargin)
% function wbdie(varargin)
%
%    standardize error format
%
% see also deprecated wberr.m
% Wb,Jan19,20

  msg=''; vflag=1; k=2;
  if nargin, l=0;
     for i=1:nargin, q=varargin{i};
        if ~isempty(q) && ischar(q) && q(1)=='-'
           if isequal(q,'-v'), vflag=2;
           elseif isequal(q,'-q'), vflag=0;
           else
              fprintf(1,'\n   wbdie(): ignoring invalid option ''%s''\n',q);
           end
           l=i;
        elseif isnumber(q), k=2+q; l=i;
        else break; end
     end
     if l<nargin
        msg=sprintf(varargin{l+1:nargin});
     end
  end

  if ~isempty(msg)
     msg=[ regexprep([10 msg],[' *' 10],[10 '   ERR ']), 10];
     use_col=wblog('--hl-check');
     if use_col
        msg=regexprep(msg,'(ERR|WRN|invalid usage)([^\n\r]*)',...
        [char(27) '[31m$1' char(27) '[38;5;9m$2' char(27) '[0m']);
     end
  end

  S=dbstack('-completenames');
  if vflag && numel(S)>1, 
     L=getcols(); if L<60 || L>99, L=80; end
     L=repmat('─',1,L-1);
     if use_col, L=[ char(27) '[38;5;8m' L char(27) '[0m']; end
     fprintf(1,L);
     dispstack(S(k:end));
  end

  S=struct('message',msg,'identifier','Wb:ERR', 'stack',S(min(2,end)));
  error(S);

end

