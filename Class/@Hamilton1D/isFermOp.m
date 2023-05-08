function isf=isFermOp(HAM,varargin)
% function q=isFermOp(HAM,op1,op2,...)
%
%    Determine whether given operator is fermionic operator
%    returning non-zero for fermionic operators, and zero otherwise
%
%    NB! the fermionic character is determined based on HAM.oez
%    hence this is only intendend for local operators.
%
% Wb,Jun05,19

  nargs=numel(varargin); isf=cell(1,nargs);
  for k=1:nargs, Ak=varargin{k};
     if ~isQSpace(Ak), wberr('invalid usage'); end
     isf{k}=zeros(size(Ak)); NK=numel(Ak);
  end

  if size(HAM.oez,1)>1
     Z=HAM.oez(2,:);
     if numel(Z)>1, wberr(...
       'not yet implemented for different site-types (stype)'); end
     Z=Z.op;

     E2=getIdentity(Z,Z); E2=setitags(E2,{'a','b','ab'});
     Z2=contract(E2,'!3*',{Z,'-op:b',{Z,'-op:a',E2}});

     dd=Z2.data;
     for i=1:numel(dd), 
        z=dd{i}; dd{i}=[ z(1), norm(z-z(1)*eye(size(z))) ];
     end
     dd=cat(1,dd{:}); e=[ norm(dd(:,2)), norm(abs(dd(:,1))-1) ];
     if any(e>1E-12), wberr('unexpected parity operator'); end

     [e,I,D]=uniquerows(round(dd(:,1)));
     if ~isequal(e,[-1;1]), wberr('unexpected parity operator'); end
     Q2=Z2.Q{1};

     fop=zeros(size(Q2,1),1); fop(I{1})=1;

     for k=1:nargs, Ak=varargin{k}; fk=nan(size(Ak)); nk=numel(Ak);
        for i=1:nk, Aki=Ak(i);
           r=numel(Aki.Q); if ~r, fk(i)=0; continue; end

           if r==2
            % if non-abelian rank-2 operators are scalars => qop=0
            % only for abelian symmetries, qop may be non-zero
            % in this case qop can be determined by difference
            % (for non-abelian symmetries this does not hold!)
            % Wb,Jun26,21
              qop=uniquerows(Aki.Q{1}-Aki.Q{2});
           elseif r==3
              qop=uniquerows(Aki.Q{3});
           else wberr('invalid rank-%g operator (arg %g, %g)',r,k+1,i);
           end

           [ia,ib,Im]=matchIndex(Q2,qop);
           if ~isempty(Im.ix2), wberr(...
               'invalid local operator (%g/%g) !?',k+1,i); end
           f=unique(fop(ia));
           if numel(f)~=1, wberr(['invalid local operator ' ...
               '(%g/%g: got mixed parity !?)'],k+1,i); end
           fk(i)=f;
        end; isf{k}=fk;
     end
  end

  if nargs==1, isf=isf{1};
  elseif all(NK==1), isf=[isf{:}]; end

end

