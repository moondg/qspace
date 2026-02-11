function wbdie(varargin)
% function wbdie(varargin)
%
%    standardize error format
%
% see also deprecated wberr.m
% Wb,Jan19,20

  msg=''; vflag=1; k=2; l=0;
  use_col=1;

  if nargin
     for i=l+1:nargin, q=varargin{i};
        if ~isempty(q) && ischar(q) && q(1)=='-', l=i;
           if isequal(q,'-v'), vflag=2;
           elseif isequal(q,'-q'), vflag=0;
           elseif isequal(q,'-k'), use_col=0;
           else
              fprintf(1,'\n   wbdie(): ignoring invalid option ''%s''\n',q);
           end
        elseif isnumber(q), l=i;
           if q<-1
                use_col=0;
           else k=2+q; end
        else break; end
     end
     if l<nargin
        msg=sprintf(varargin{l+1:nargin});
     end
  end

  S=dbstack('-completenames');

  if use_col
     if numel(S)>32, use_col=0;
     elseif vflag || ~isempty(msg)
        use_col=wblog('--hl-check');
     end
  end

  if ~isempty(msg)
     msg=[ regexprep([10 msg],[' *' 10],[10 '   ERR ']), 10];
     if use_col
        msg=regexprep(msg,'(ERR|WRN|invalid usage)([^\n\r]*)',...
        [char(27) '[31m$1' char(27) '[38;5;9m$2' char(27) '[0m']);
     end
  end

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

