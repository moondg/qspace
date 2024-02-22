function [A,I]=squeeze(A,varargin)
% function [A,I]=squeeze(A [,OPTS])
%
%    remove singleton dimensions (similar to matlab's routine)
%
% Options
%
%   '-start'  only singletons at the front
%   '-end'    only singletons at the end (default: all)
%
% Wb,Apr24,08 ; Wb,Dec01,20

  getopt('init',varargin);
     if getopt('-start'), ks='start';
     elseif getopt('-end'), ks='end'; else ks=[]; end
  ks=getopt('get_last',ks);

  nA=numel(A);
  for k=1:nA, [A(k),I{k}]=squeeze_1(A(k),ks); end
  if nA==1, I=I{1}; end

end

% -------------------------------------------------------------------- %
% only indizes with q=0 (scalars) are considered singletons

function [A,ks]=squeeze_1(A,ks)

  Q=A.Q; r=numel(Q); qnrm=zeros(1,r);
  for i=1:r, qnrm(i)=norm(Q{i},'fro'); end

  if isempty(ks), ks=find(qnrm==0);
  elseif ~isnumeric(ks)
     if isequal(ks,'start')
        for i=1:r,    if qnrm(i), ks=1:i-1; break; end; end
     elseif isequal(ks,'end')
        for i=r:-1:1, if qnrm(i), ks=i+1:r; break; end; end
     else ks, wbdie('invalid usage'); end
     if any(ks<1 || ks>r), ks, wbdie('index out of range (r=%g)',r); end
  else
     ks=unique(ks);
  end

  if isempty(ks), return; end

  qdir=getqdir(A); qk=qdir(ks); nk=numel(ks);
  if any(qnrm(ks))
     wbdie('invalid usage (specifed non-scalar dimension)');
  end

  if nk==1
     v=getvac(A,'-1d'); if qk>0, ic='1*'; else ic=1; end
     A=contract(A,ks,v,ic);
     return
  end

  E=getvac(A); Z=getIdentity(E,'-0');
  if nk==2
     if qk(1)>0; ic='12*'; else ic='12'; end
     if diff(qk)
          A=contract(A,ks,E,ic);
     else A=contract(A,ks,Z,ic); end
     return
  end

  E3=getIdentity(E,1,E,1); B=E3;
  for i=4:nk, B=contract(B,i-1,E3,1); end

  if qk(nk)>0, B=contract(B,nk,Z,'1'); end
  for i=find(qk<0), if i<nk
     B=contract(B,i,Z,'1*',initperm(nk,'--last2',i)); end
  end

  A=contract(A,ks,B,[sprintf('%g',1:nk) '*']);

end

% -------------------------------------------------------------------- %

