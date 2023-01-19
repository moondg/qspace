function qs=toStr(A,varargin)

  sA=size(A); n=prod(sA); ii=1:n;

  if ~isempty(inputname(1))
  vname=inputname(1); else vname='ans'; end

  switch nargin
     case 1
     case 2, if ischar(varargin{1})
                  vname=varargin{1};
             else ii=varargin{1}; end
     case 3, ii=varargin{1}; vname=varargin{2}; 
  otherwise, eval(['help ' mfilename]); return;
  end

  if length(sA)==2 && sA(2)==1, sA=sA(1); end

  if n==1
     qs{1}=sprintf('%s = %s', vname, toStr_1(A));
  elseif n>1
     for i=1:length(ii)
        idx=vec2str(ind2sub_aux(sA,ii(i)),'sep',',');
        qs{i}=sprintf('%s(%s) = %s',vname,idx, toStr_1(A(ii(i))));
     end
  else
     qs{1}=sprintf('%s = (empty QSpace)', vname);
  end

  nl=sprintf('\n');
  for i=1:length(qs), qs{2,i}=nl; end
  qs=cat(2,qs{:});

  for i=length(qs):-1:1
  if qs(i)>' ', break; end, end

  qs=qs(1:i);

end

% -------------------------------------------------------------------- %
function kk = ind2sub_aux(s,k)

  k=k-1; kk=zeros(size(s));

  for i=1:length(s)
  kk(i)=mod(k,s(i)); k=(k-kk(i))/s(i); end

  kk=kk+1;

end

% -------------------------------------------------------------------- %

function qs=toStr_1(A)

  nd=length(A.data); nq=length(A.Q);
  fmt=sprintf('%%-%ds', 8*(nd>1)+5*length(A.Q));

  for i=1:nd
      sa=size(A.data{i}); s=prod(sa);

      for k=1:nq, qq(k,:)=A.Q{k}(i,:); end
      qq=mat2str2(qq,'fmt','%2g','rowsep','; ');

      if nd<2 && s<24 && length(sa)==2 && size(sa,2)<6
         if sa(1)==1
              qs=sprintf('{%s} [%s ]\n', qq, mat2str2(A.data{i},'fmt','%4g'));
         else qs=sprintf('[ // {%s}\n%s\n]\n\n',...
              qq, mat2str2(A.data{i},'fmt','%4g'));
         end
         continue
      end

      if nd<2, s=''; else s=sprintf('%2d.  ',i); end
      qs=sprintf(fmt, sprintf('%s%s   ', s, ...
      vec2str(sa,'fmt','%2d','sep',' x ')));
      qs=[qs, sprintf('[%s ]  ',qq)];

      sa=prod(sa);
      if sa==1
         qs=[qs sprintf('%8g\n', A.data{i}(1))];
         continue;
      end

      sa=sa*8;

      if     sa<2^10, s=sprintf('%g b',sa);
      elseif sa<2^20, s=sprintf('%.1f k',sa/2^10);
      else            s=sprintf('%.1f M',sa/2^20);
      end

      qs=[qs, sprintf('%7s\n',s)];
  end

  qs=qs(1:end-1);

end

% -------------------------------------------------------------------- %
