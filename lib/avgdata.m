function ym=avgdata(y,varargin)
% function ym=avgdata(y [,m, opts])
% 
%    Average data over m points (default m=2).
%
% Options
%
%    '-l'  preserve length (by default, length is reduced by (m-1))
% 
% Wb,Nov17,99 ; Wb,Apr05,13 ; Wb,Mar15,25

  if ~nargin
     if ~helpthis(nargout,varargin{:}), wbdie('invalid usage'); end
     return
  end

  lflag=0; m=2;

  if nargin>1, m=varargin{1};
     if nargin>2
        getopt('init',varargin(2:end));
           lflag=getopt('-l');
        getopt('check_error');
     end
     if m==1, ym=y; return; end
     if mod(m,1) || m<2, wbdie('invalid usage (m=%g)',m); end
  end

  transp=(diff(size(y))>0);
  if transp, y=y.'; end

  [N,n]=size(y);
  if m>N, wbdie('too few data points for m=%g (N=%g)',m,N); end

  if ~lflag
     ym = y(1:N-m+1,:);
     for i=2:m, ym=ym+y(i:N-m+i,:); end
     ym=ym/m;
  elseif N>2
     if 2*m>N+1, m=floor(N/2); end

     w=linspace(0,1,m+1);
     w=[ w(2:end) flip(w(2:end-1)) ];
     w=w/sum(w);

     ym = w(m)*y(m:N-m+1,:);
     for i=1:m-1
        ym = ym + w(m+i)*( y(m+i:(N-m+1)+i,:) + y(m-i:(N-m+1)-i,:) );
     end

     p=min(3,m-1); nb=2*m;
     Y1=[ zeros(m-1,n);  y(1:nb,:)    ];
     Y2=[ y(N-nb+1:N,:); zeros(m-1,n) ];
     for j=1:n
        pj=polyfit(-1:nb-2,Y1(m:end,j),p); Y1(1:m-1,j)=polyval(pj,-m:-2);
        pj=polyfit(-nb+2:1,Y2(1:nb,j),p);  Y2(nb+1:end,j)=polyval(pj,2:m);
     end

     i1=m;      i1m=i1+m-2; y1 = w(m)*Y1(i1:i1m,:);
     i2=nb-m+2; i2m=nb;     y2 = w(m)*Y2(i2:i2m,:);
     for i=1:m-1
        y1 = y1 + w(m+i)*( Y1(i1+i:i1m+i,:) + Y1(i1-i:i1m-i,:) );
        y2 = y2 + w(m+i)*( Y2(i2+i:i2m+i,:) + Y2(i2-i:i2m-i,:) );
     end

     ym=[y1; ym; y2];
  else
     ym=y;
  end

  if transp, ym=ym.'; end

end

