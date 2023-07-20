function yn=avgdata(y,n,lflag)
% Function: y_new = avgdata(y,n [,'-l'])
% 
%    averages data over n points.
%    If '-l' is specified, length of data is preserved.
% 
% Wb,Nov17,99 ; Wb,Apr05,13

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  if nargin<2, n=2; end
  if nargin<3, lflag=0;
    elseif isequal(lflag,'-l'), lflag=1;
    else wbdie('invalid usage'); 
  end

  tflag=size(y,1)<size(y,2);
  if tflag, y=y'; end

  if n<=1, yn=y;
     if n~=1, wblog('WRN','got n=%g !??',n); end
     return
  end

  if ~lflag
     yn = y(1:end-n+1,:);
     for i=2:n, yn=yn+y(i:end-n+i,:); end
     yn=yn/n;
  else
     w=0:(1/n):1; w=[ w(2:end) fliplr(w(2:end-1)) ];
     w=w/sum(w);

     if size(y,1)<2*n, wbdie(...
       'invalid usage (too few data points for n=%g (%g)',n,size(y,1));
     end

     yn = w(n)*y(n:end-n+1,:);
     for i=1:n-1, yn = yn + w(n+i)*[ 
        y(n+i:(end-n+1)+i,:) + ...
        y(n-i:(end-n+1)-i,:)];
     end

     we=zeros(n-1,3*n-3); l=2*n-2; k=[n, 2*n-1, 3*n-2];
     for i=1:n-1
        we(i,i:i+l)=w;
     end
     we=[ we(:,1:k(1)-1), we(:,k(1):k(2)-1)+flipud(we(:,k(2):k(3)-1)') ];

     yn=[ flipud(fliplr(we)*y(1:2*n-2,:))
         yn
         we*y(end-2*n+3:end,:)
     ];
  end

  if tflag, yn=yn'; end

end

