function fmt=getqfmt(A,varargin)
% function fmt=getqfmt(A [,opts])
%
%    get format string for a single Q{i}(j,:)
%
% Wb,Dec14,15

% see also display -> get_q_fmt(qtype,r,m)

  if isempty(A), fmt=''; return; end
  getopt('init',varargin);
     sep  =getopt('sep',' ');
     bflag=getopt('-b');
  getopt('check_error');

  d=0; fmt=reshape(getsym(A,'-c'),1,[]);

  Q=abs(cat(1,A.Q{:}));

  for i=1:numel(fmt), s=fmt{i};
     if ~ischar(s), wbdie('got invalid symmetry (string required) !?'); end
     if regexp(s,'^SU\d+$')
        r=str2num(s(3:end))-1; d=d+r;
        fmt{i}=repmat('%X',1,r);
     elseif regexp(s,'^Sp\d+$')
        r=str2num(s(3:end))/2; d=d+r;
        fmt{i}=repmat('%X',1,r);
     else
        d=d+1;
        if any(Q(:,d)>10), fmt{i}='%3g'; else fmt{i}='%2g'; end
     end
  end

  if d~=size(Q,2), wbdie('qset mismatch (len=%g/%g)',d,size(Q,2)); end

  fmt(2,1:end-1)={sep};
  fmt=[fmt{:}]; if bflag, fmt=['(' fmt ')']; end

end

