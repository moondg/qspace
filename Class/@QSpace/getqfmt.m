function [qfmt,Q3]=getqfmt(A)
% function [qfmt,Q3]=getqfmt(A)
%
%    Get format string to print q-labels in A.Q{l}.
%    This uses compact notation for SU(N), Sp(2n), etc
%    as along as q-labels stay with in hexadecimal range (q<16),
%    and extend number format otherwise.
%
% Wb,Apr02,15

% merged with private/get_q_fmt.m // Wb,Aug13,24

  if ~isempty(A.info), qtype=A.info.qtype;
  elseif isempty(A.Q), qfmt=''; return
  else qtype=''; end

  rA=numel(A.Q);
  if rA, nQ=size(A.Q{1},2); else nQ=0; end

  sym=strread(qtype,'%s','delimiter',',;')';
  if isempty(sym)
     qfmt=strjoin(repmat({'%2g'},1,nQ),' ');  if rA>1
     qfmt=strjoin(repmat({qfmt},1,rA),' ;');  end
     return
  end

  QQ=cat(1,A.Q{:}); gotQ=~isempty(QQ); l=1; r=0; ext=0;
  Qmin=min(QQ,[],1);
  Qmax=max(QQ,[],1);

  nsym=numel(sym);
  qfmt=repmat({''},2,nsym);

  for j=1:nsym, l=l+r; r=1; compact=0;
     i=regexp(sym{j},'\d+$');
     if ~isempty(i)
        w=sym{j}(1:i-1); n=str2num(sym{j}(i:end));
        switch w
           case 'SU', compact=1; r=n-1;
           case 'Sp', compact=1; r=n/2;
           case 'SO', compact=1; r=floor(n/2); 
           case {'Z','P'},       r=1;
           otherwise wbdie('unexpected symmetry %s',sym{j});
        end
     end
     if compact
        if all(Qmax(l:l+r-1)<16)
             qfmt{1,j}=repmat('%X',1,r);
        else qfmt{1,j}=strjoin(repmat({'%2g'},1,r),' '); ext=1;
        end
     elseif Qmax(l)<10  && Qmin(l)>=0,  qfmt{1,j}='%2g';
     elseif Qmax(l)<100 && Qmin(l)>-10, qfmt{1,j}='%3g';
     else                               qfmt{1,j}='%4g';
     end
  end

  if l+r-1~=nQ, wbdie('qset mismatch (nQ=%g/%g)',l,nQ); end

  if ext, qfmt(2,1:nsym-1)={', '}; else qfmt(2,1:nsym-1)={' '}; end
  qfmt=[qfmt{:}];

  if nargout>1
     Q3=permute(cat(3,A.Q{:}),[3 2 1]);
     qfmt=strjoin(repmat({qfmt},1,rA),' ; ');
  end

end

