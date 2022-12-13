function [qd,qa]=get_defining(sym)
% function [qd,qa]=get_defining(sym)
%
%    Get defining (and adjoint) representation for sym.
%
%    (if 2nd return argument is specified, this returns
%    the adjoint representation, which results out
%    of the tensor product of the defining representation
%    with its dual).
%
% Wb,Oct05,15

  if regexp(sym,'SU\d+$'), N=str2num(sym(3:end));
     q=repmat('0',1,N-1); q(1)='1';
     q_=fliplr(q);
  elseif regexp(sym,'Sp\d+$'), N=str2num(sym(3:end))/2;
     q=repmat('0',1,N); q(1)='1';
     q_=q;
  else sym, error('Wb:ERR','\n   ERR invalid symmetry');
  end

  qd=LoadRData(sym,q);

  if nargout>1
     I3=LoadMData(sym,[q,',' q_],'-q');
     jj=I3.J; jj(sum(jj.^2,2)==0,:)=[];

     if size(jj,1)>1
        if regexp(sym,'Sp\d+$')
           i=find(jj(:,1)==2); if numel(i)==1, jj=jj(i,:); end
        end
     end

     if size(jj,1)~=1, error('Wb:ERR',...
       '\n   ERR failed to determine adjoint representation'); end

     qa=LoadRData(sym,sprintf('%g',jj));
  end

end

