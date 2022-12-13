function [p,pstr,jj]=subsref(P,SS)
% Function [p,pstr,ip]=subsref(P,S)
%
%    Overloading the subsref() routine for PSet, as in 
%       P(i)   returns the parameters for given index 1<=i<=P.n
%       P()    returns the parameters for internal iterator i=P.i
%
%    Note that this routine does not alter P in any way. In particular,
%    this does not keep track of time stamps. For the latter, e.g.,
%    together with loops, use
%
%       while next(P), [p,pstr]=P(); structexp(p);
%          ...
%       end
%
%   Alternatively, may access specific parameter index
%
%       [iz,nz]=P(ip,'B');
%
%   NB! using col-major, i.e. first index (parameter) is fastest
%   such that OUT with OUT(i) <- P(i) can be simply reshaped using P.s.
%
% Output
%    
%    p     specific parameters setting for given index ip
%          returned as structure (array, if length(ip)>1)
%    pstr  info string that can be used for logging (for scalar ip only).
%    ip    sub index for given parameters (using ind2sub)
%
% Wb,Jul31,08 ; Wb,Mar13,19

  S=SS(1); ns=length(SS);

  if isequal(S.type,'()'), nargs=length(S.subs);
     if ~nargs, S.subs{1}=P.i;
     else
        if nargs>2 || ~isnumeric(S.subs{1}) || ...
           nargs>1 && ~ischar(S.subs{2}), S.subs
           error('Wb:ERR','invalid usage');
        end
     end

     k=reshape(S.subs{1},1,[]); nk=length(k); k_=k;
     jj=zeros(P.r,nk);

     if any(k>P.n) || isempty(k)
        error('Wb:ERR','index out of range (%g/%g)',max(k),P.n); end

     jj=cell(1,numel(P.s)); [jj{:}]=ind2sub(P.s,k_);
     jj=cat(2,jj{:});

     if nargs>1, f=S.subs{2};
        p=0; pstr=0;
        for i=1:P.r, if isequal(P.vars{i},f)
            p=jj(:,i); pstr=P.s(i); break;
        end,end
        return 
     end

     s=P.vars;
     for i=size(jj,1):-1:1
         for j=1:P.r, s{2,j}=P.data{j}(jj(i,j)); end
         p(i)=struct(s{:});
     end

     if nargout<=1
        if ns>1, p=builtin('subsref',p,SS(2:end)); end
        return
     end

     if nk>1 || ns>1, wblog('WRN',...
        'pstr only specified for SINGLE input index (%g)',nk); end

     s(4,:)=s(2,:);
     s(2,:)=mat2cell(jj(1,:),1,ones(1,P.r));
     s(3,:)=mat2cell(P.s,    1,ones(1,P.r));

     for i=1:P.r
        if s{3,i}>1
             s{1,i}=sprintf('%s(%g/%g)=%.4g  ', s{:,i});
        else s{1,i}=sprintf('%s=%.4g  ', s{[1 4],i}); end
     end

     pstr=[s{1,:}]; pstr=pstr(1:end-2);

     if numel(k_)~=1, return; end

     if sum(P.s>1)>1, pstr=[ sprintf('%g/%g:  ',k_,P.n), pstr ]; end

  elseif isequal(S.type,'.') && ischar(S.subs)

     if builtin('isfield',struct(P),S.subs)
        p=builtin('subsref',P, S);
        if ns>1, p=builtin('subsref',p,SS(2:end)); end
        return
     end

     n=S.subs;
     for i=1:P.r
        if isequal(P.vars{i},n), p=P.data{i};
           if ns>1, p=builtin('subsref',p,SS(2:end)); end
           return
        end
     end

     error('Wb:ERR','var ''%s'' not declared within PSet',S.subs); 

  else
     p=builtin('subsref',P,SS);
  end

end

