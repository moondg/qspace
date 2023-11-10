function s=itags_to_str(tt,use_col,varargin)
% function s=itags_to_str(tt [,use_col,...])
% Wb,Feb24,16

% outsourced from info.m
  if nargin>1 && ~isempty(use_col)
     if isnumeric(use_col) && use_col
        [use_col,e1,em]=syntax_hl('QS:info',varargin{:});
     elseif isstr(use_col)
        [use_col,e1,em]=syntax_hl( use_col, varargin{:});
     else wbdie('invalid usage'); end
     if use_col
        if ischar(e1)
             use_col=1; e1={e1};
        else use_col=3; end
     end
  else use_col=0; end

  if isempty(tt), s='';
  elseif iscell(tt), notags=1;
     for i=1:numel(tt) % permit mark ' and conj flag *
        if ~isempty(tt{i}) && isempty(regexp(tt{i},'^''?\*?$'))
           notags=0; break
        end
     end
     if notags
        for i=1:numel(tt) % permit mark ' and conj flag *
           if isempty(tt{i}), tt{i}='+';
           elseif tt{i}(end)~='*', tt{i}=['+' tt{i}];
           else tt{i}=['-' tt{i}(1:end-1)]; end
        end
     end

     if use_col
        for i=1:numel(tt), t=tt{i}; q=0;
           if ~isempty(t) && t(end)=='*',  t(end)=[]; q=bitor(q,2); end
           if ~isempty(t) && t(end)=='''', t(end)=[]; q=bitor(q,1); end
           if use_col==1 && bitand(q,2)
              t=[t '*']; q=bitand(q,1);
           end
           if q % make sure at least empty string '' is shown
              if isempty(t), t=''''''; end
              tt{i}=[e1{q} t em];
           elseif isempty(t), tt{i}=''''''; end
        end
     end
     if notags
        s=['{ ' tt{:} ' }'];
     else
        tt=reshape(tt,1,[]); tt(2,1:end-1)={', '};
        s=['{ ' tt{:} ' }'];
     end
  elseif ischar(tt)
     s=['{ ' regexprep(tt,'[,;| ]+',',') ' }'];
  else 
     disp(tt), wbdie('invalid itags');
  end

end

