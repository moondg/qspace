function [ss,Is]=calc_SdotS(HAM,varargin)
% function [ss,Is]=calc_SdotS(HAM [,dx])
%
%    calculate local (rank-2 ops) or nearest-neighbor (rank-3 ops)
%    expectation values (at distance dx; default: dx=1, i.e. NN)
%
% Options
%
%    'z',...   fermionic signs (specify index in QIDX to use)
%
% Wb,May12,17

  getopt('init',varargin);
     zz=getopt('z',[]);
  dx_in=getopt('get_last',[]);

  L=numel(HAM.mpo);
  [kc,Ic]=getCurrentSite(HAM);

  RR=QSpace(1,L);

  if isfield(HAM.info,'lops')
     Is.istr='using HAM.info.lops';
     ops=HAM.info.lops; nops=numel(ops);
  else
     Is.istr='using HAM.ops';
     ops=QSpace(size(HAM.ops)); nops=numel(ops);
     for i=1:nops, ops(i)=HAM.ops(i).op; end
  end

  ops=untag(ops); 
  rop=zeros(1,nops); for i=1:nops, rop(i)=rank(ops(i)); end

  dx=ones(1,nops);
  dx(find(rop<3))=0;

  zz(end+1:nops)=0;

  if isfield(HAM.info,'param') ...
  && isfield(HAM.info.param,'perBC') && HAM.info.param.perBC
     dx(find(dx))=2;
  else 
     q=abs(HAM.info.HH(:,3)-HAM.info.HH(:,1));
     q=full(sparse(ones(size(q)),q+1,ones(size(q))));
     if numel(q)>2 && q(3)>2*q(2)
        dx(find(dx))=2;
     end
  end

  if ~isempty(dx_in)
     if numel(dx_in)~=nops, wberr(...
       'invalid usage [numel(dx) does not match nops=%g]',nops); end

     i=find(dx_in==0);
     if any(rop(i)~=2), wberr('got local operators of rank~=2 !?'); end

     dx_=dx; dx=dx_in;
  end

  Ak=load_dmrg_data(HAM,kc,'AK');
  NPsi=0; tk=getitags(Ak);

  if numel(tk)>3
     q=getDimQS(Ak); d=q(:,4); NPsi=d(end); r=numel(tk);
     q=sprintf(',%s',tk{:}); s=sprintf('%s; %g states',q(2:end),NPsi);
     if r==4 && ~isempty(regexp(tk{end},'Psi')), q=iff(d(1)~=1,'s','');
          wblog('NB!','targeting %g multiplet%s (%s)',d(1),q,s);
     else wblog('WRN','got rank-%g QSpace (%s)',r,s);
     end
     if NPsi<=1, wberr('got NPsi=%g !?',NPsi); end

     EPsi=getIdentity3(Ak,4,'PSI','--rho'); % tag 'PSI', '--rho' => normalize

     ss=cell(L,nops);
  else
     ss=nan(L,nops);
  end

  if kc>2
     fprintf(1,'%32s (kc=%g)\r','',kc);
     for k=kc:-1:2, if k<kc,
        Ak=load_dmrg_data(HAM,k,'AK'); end
        fprintf(1,'\r   %s << %d/%d ...  \r',mfilename,k,kc);

        if     k<kc, Q1={Ak,RR(k+1)};
        elseif NPsi, Q1={Ak,EPsi};
        else         Q1= Ak; end

        RR(k)=contract(Ak,'!1*',Q1);
     end
     fprintf(1,'\n');
  end

  SS=QSpace(nops,max(dx)); Ak=QSpace; tk={};
  rho=QSpace; Rflag=0;

  for k=1:L
     fprintf(1,'\r   %s >> %d/%d ...  \r',mfilename,k,L);
     q=load_dmrg_data(HAM,k); Al=Ak; Ak=q.AK; tl=tk; tk=getitags(Ak);

     if k==1, dk=getDimQS(Ak);
        if dk(1,end)>1, Rflag=1; rho=QSpace(L,2);
        elseif isfield(q.info,'rho'), rho=q.info.rho(end);
        end
     end
     if Rflag && isfield(q.info,'rho')
        rho(k,2)=q.info.rho(end);
     end

     if numel(tk)>3 && k~=kc
        r=numel(tk); q=sprintf(',%s',tk{:});
        wberr('got rank-%g QSpace (%s)',r,s(2:end));
     end

     QL=Ak; QR=Ak;

     if     k>kc, Qk=contract(Ak,RR(k-1)); QL=Qk;
     elseif k<kc, Qk=contract(Ak,RR(k+1)); QR=Qk;
     elseif NPsi, Qk={Ak,EPsi}; QL=Qk; QR=Qk;
     else   Qk=Ak; end

     if Rflag
        rho(k,1)=contractQS(Ak,'!3',Qk,'*');
     end

     RR(k)=contract(Ak,'!2*',Qk);

     for i=1:nops
        if dx(i)==0
           X=contractQS(Ak,'*',{Qk,ops(i),['-op:' tk{3}]});
        elseif k>1
           q=ops(i); if zz(i), q=applyParity(q,1,zz(i)); end
           X=contractQS(Ak,'*',{SS(i,1),'*',{QR,q,['-op:' tk{3}]}});
        else X=QSpace;
        end
        if NPsi
           if ~isempty(X.data)
              if numel(X.data)~=1, X, wbdie('unexpected X.data'); end
              ss{k,i}=X.data{1};
           end
        else
           ss(k,i)=getscalar(QSpace(X));
        end
     end

     for i=1:nops
        if dx(i)
           if k>1
              for j=1:dx(i)-1
                 Q=Ak; if zz(i), Q=applyParity(Q,3,zz(i)); end
                 if k==kc && NPsi, Q={Q,EPsi}; end
                 SS(i,j)=contract(Ak,'!2*',{SS(i,j+1),Q});
              end
           end
           j=dx(i);
           SS(i,j)=contract(Ak,'!2*',{QL,ops(i),['-op:' tk{3}]});
        end
     end

  end
  if kc>0
       fprintf(1,'\n');
  else fprintf(1,'\r%60s\r',''); end

  if iscell(ss) && max(dx)<L
     for i=1:nops
        for k=1:dx(i)
           if isempty(ss{k,i}), ss{k,i}=zeros(size(ss{dx(i)+1,i})); end
        end
     end

     ss=permute(cat(3,ss{:}),[3 1 2]); s=size(ss);
     ss=reshape(ss,[L nops, prod(s(2:end))]);

  elseif NPsi && ~iscell(ss)
     ss=ss/NPsi;
  end

  if nargout>1
     Is=add2struct(Is,ops,rop,dx,kc,rho,NPsi);
  end

end

