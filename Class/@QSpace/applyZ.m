function A=applyZ(A,idx,varargin)
% function A=applyZ(A,idx [,iq,opts])
%
%    apply fermionic signs to given index / indices (idx)
%    with fermionic charge (or parity Z2) assumed at position
%    iq in the Q-labels (iq may be skipped, if only one
%    abelian symmetry label is present).
%
% Options
%
%   'qfac',.. : factor to apply on charge labels, to obtain
%    actual charge (must be integer; enforced); e.g. note that
%    single-flavor spinless fermions have an artificial factor
%    of 2 on their charge labels to permit simple symmetry around
%    half-filling; hence in this case: applyZ(...,'qfac',1/2);
%
% Wb,Sep17,18

  getopt('init',varargin);
     qfac=getopt('qfac',0);
  iq=getopt('get_last',[]);

  n=numel(idx); icheck=1;

  for k=1:numel(A)
     if ~isempty(A(k).data), Ak=A(k);
        if icheck, icheck=0;
           r=getsym(Ak,'-r'); d=getsym(Ak,'-d');
           if isempty(iq), i=find(r==0);
              if numel(i)==1, iq=sum(d(1:i));
              else error('Wb:ERR','\n   ERR invalid usage (missign iq)');
              end
           elseif numel(iq)~=1 || isempty(find(cumsum(d)==iq)) || r(iq)
              iq, error('Wb:ERR','\n   ERR invalid index iq !?'); 
           end
        end

        r=[numel(Ak.Q), size(Ak.Q{1},2)];
        if any(idx>r(1)) || iq>r(2)
           error('Wb:ERR','\n   ERR index out of bounds (%s/%d; %d/%d) !?',...
           vec2str(idx,'-f'),r(1),iq,r(2));
        end

        if n==1, q=Ak.Q{idx}(:,iq); else
           q=cell(1,n); for i=1:n, q{i}=Ak.Q{idx(i)}(:,iq); end
           q=sum([q{:}],2);
        end

        if qfac
           q=q*qfac; e=norm(q-round(q));
           if e>1E-8, error('Wb:ERR',...
             '\n   ERR got non-integer charge labels !?');
           end
        end

        for i=find(mod(q',2)), Ak.data{i}=-Ak.data{i}; end
        A(k)=Ak;
     end
  end

end

