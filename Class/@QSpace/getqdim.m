function [dd,s]=getqdim(A)
% function [dd,ss]=getqdim(A)
%
%   get symmetry specific length of q-labels
%
% Wb,Sep18,11

  if nargin<1 || numel(A)~=1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(A.Q), dd=[]; return; end

  m=size(A.Q{1},2);
  if isempty(A.info) || isempty(A.info.cgr)
     dd=ones(1,m);
     if nargout, s=repmat({'A'},1,m); end
     return
  end

  q=A.info.qtype;
  s=strread(q,'%s','delimiter',',');
  n=numel(s); dd=nan(1,n);

  for i=1:n
     switch s{i}
        case 'A', dd(i)=1;
        case 'SU2', dd(i)=1; % s{i}='SU(2)';
        case 'SU3', dd(i)=2; % s{i}='SU(3)';
        case 'SU4', dd(i)=3; % s{i}='SU(3)';
        otherwise error('Wb:ERR','\n   ERR unknown symmetry ''%s''',s{i}); 
     end
  end

  if sum(dd)~=m, dd
     error('Wb:ERR','\n   ERR symmetry inconsistency (%g/%g)',m,sum(dd));
  end

end

