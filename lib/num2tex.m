function s=num2tex(x,varargin)
% Function: num2tex(x [,'%..e'])
%
%    print number in %e format (%.4G)
%    replacing the E+NN part by 10^{NN}
%
%    NB! use capital %.3E to enforce exponent.
%
% Wb,Mar07,13 : extended to vectors and matrices
% Wb,Oct10,05

  getopt('init',varargin);
     sep  =getopt('sep',{', ', '\n'});

     if getopt('-th'),    xflag=1;
     elseif getopt('-x'), xflag=2;
     else                 xflag=0; end

     len  =getopt('len',[]);
  fmt=getopt('get_last','%.4G');

  n=numel(x); 

  if xflag==1
     if n~=1 || x<0, wbdie('invalid usage'); end
     if     x==1, s='st';
     elseif x==2, s='nd';
     elseif x==3, s='rd'; else s='th'; end

     s=sprintf('%g^{%s}',x,s);
     return
  end

  s=cell(size(x));
  for i=1:n
     s{i}=num2tex_1(x(i),fmt,xflag,len);
  end

  if ~iscell(sep), sep={sep}; end
  if ~isempty(strfind(sep{1},'\')), sep{1}=sprintf(sep{1}); end

  if n==1, s=s{1};
  elseif n>1 && numel(find(size(x)>1))==1
     s=join_to_string(s,sep{1});
  else
     if numel(sep)==1, sep={sep{:},'\n'}; end
     if ~isempty(strfind(sep{2},'\')), sep{2}=sprintf(sep{2}); end
     for i=1:size(x)
        s{i}=join_to_string(s(i,:),sep{1});
     end
     s=join_to_string(s(:,1),sep{2});
  end

end

% -------------------------------------------------------------------- %

function s=join_to_string(s,sep)
  s=reshape(s,1,[]);
  s(2,:)={sep}; s{end}='';
  s=cat(2,s{:});
end

% -------------------------------------------------------------------- %

function s=num2tex_1(x,fmt,xflag,len)

  if x==0, s=sprintf('%g',x); return; end

  if ~ischar(x)
     i=find(fmt=='e');
     if ~isempty(i) && abs(log10(abs(x)))<3
        fmt(i)='g';
     end
     s=sprintf(fmt,x);
  else s=x; end

  s=regexprep(upper(s),'E[+-]00','');

  hasE=any(find(s=='E'));
  hasdot=any(find(s=='.'));

  if hasdot && ~xflag
     s=regexprep(s,'[.0]*([E+-])0*','$1');
     if ~hasE, s=regexprep(s,'[.0]*$',''); end
  end

  if hasE
     s=regexprep(s,'E\+0*([0-9]*)','\\cdot10^{$1}');
     s=regexprep(s,'E\-0*([0-9]*)','\\cdot10^{-$1}');
     if ~xflag
        s=regexprep(s,'^1\.0*\\cdot','');
        s=regexprep(s,'^1\\cdot','');
     end
  end

  if ~isempty(len)
     l=length(regexprep(regexprep(s,'\\[a-z]*','.'),'[{}]',''));
     if len && l<len
        s=[repmat(' ',1,len-l), s];
     end
  end

end

% -------------------------------------------------------------------- %

