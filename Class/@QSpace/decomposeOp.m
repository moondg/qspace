function [cc,err,I]=decomposeOp(H,X,varargin)
% function [cc,err,I]=decomposeOp(H,{O1,O2,...}[,R][,opts])
%
%    decompose operator H into operator (sets) specified in X={O1, O2, ...};
%    if R is specified, operator overerlap is calculated
%    w.r.t. density matrix R.
%
% Wb,Apr16,13

  if nargin<2  || isempty(X)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     qflag=getopt('-q');
  R=getopt('get_last',[]);

  if isQSpace(X), X=X(:);
  else
     for i=1:numel(X)
        if ~isQSpace(X{i}), error('Wb:ERR',...
          '\n   ERR invalid usage (X must contain QSpaces only)'); end
        X{i}=X{i}(:);
     end
     X=cat(1,X{:});
  end
  nops=numel(X);

  if ~isempty(R)
     if numel(R)~=1, error('Wb:ERR','\n   ERR invalid projector'); end
     e=norm(R-R'); if e>1E-14
        wblog('WRN got nonsymmetric projector !?? (%.3g)',e);
     end
     q=trace(R);
     if abs(q-1)>1E-12
        if abs(q-1)>0.1
           wblog('got unnormalized density matrix (%.3g) !??',q); end
        R=R*(1/q);
     end
  end

  S=zeros(nops,nops); b=zeros(nops,1);

  for i=1:nops
      if isempty(X(i)), error('Wb:ERR',['\n   ' ... 
        'ERR invalid operator basis (got empty QSpace)']); end
      if isempty(R)
         b(i)=getscalar(QSpace(contractQS(X(i),[1 2],H,[1 2])));
      else 
         b(i)=trace(R*(X(i)'*H));
      end
  end

  for i=1:nops
  for j=i:nops
      if isempty(R)
           S(i,j)=getscalar(QSpace(contractQS(X(i),[1 2],X(j),[1 2])));
      else S(i,j)=trace(R*(X(i)'*X(j))); end
      if j>i, S(j,i)=S(i,j);
      elseif i==j && abs(S(i,j))<1E-12
         error('Wb:ERR','\n   ERR invalid operator basis (empty QSpace?)'); 
      end
  end
  end

  if ~isreal(S), error('Wb:ERR',...
     '\n   ERR invalid usage (got complex QSpaces!?)');
  end

  if qflag
     warning off MATLAB:nearlySingularMatrix
  end

  cc=S\b; 

  if qflag
     warning on MATLAB:nearlySingularMatrix
  end

  if nargout>1
     Hc=cc(1)*X(1); for i=2:nops, Hc=Hc+cc(i)*X(i); end
     Q=H-Hc;

     if isempty(R)
          err=norm(Q)/norm(H);
     else err=trace(R*(Q'*Q)) / trace(R*(H'*H));
     end

     if nargout>2, I=add2struct('-',S,b,Hc,R); end
  end

end

