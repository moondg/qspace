function [ss,Iout]=calc_ScorS(HAM,varargin)
% function [ss,Iout]=calc_ScorS(HAM [,k0])
%
%    calculate static correlation function vs. distance
%    relative to (center) site k0; if the operator is a scalar
%    operator of rank-2 (such as Sz), it's expectation values
%    along the chain are also calculated.
% 
%    The operators used are taken 
%    in the following order of preference (whichever exists first)
%       1) HAM.info.xops (NB! tdDMRG uses lops, hence make sure
%          that xops are properly defined; `x'ops by their name,
%          refer to operators to be used with correlations)
%       2) HAM.info.lops: using split_ops(lops)
%       3) HAM.ops(:).op (default if 1 & 2 are not set)
%
% Options
% 
%   '-k'   stop with keyboard at end of routine
%   '-v'   verbose mode (report progress)
% 
% Wb,Jul06,19

% adapted from calc_SdotS
% load Wb190705-db8_tst_tdDMRG.mat

  L=numel(HAM.mpo); ss=[];

  getopt('init',varargin);
     kflag=getopt('-k');
     vflag=getopt('-v');
  k0=getopt('get_last',[]);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);

  RR=QSpace(1,L); lops=1;

  if isfield(HAM.info,'xops')
     istr='using HAM.info.xops';
     ops=HAM.info.xops(:); nops=numel(ops);
  elseif isfield(HAM.info,'lops')
     istr='using HAM.info.lops';
     [ops,n]=split_ops(HAM.info.lops,'--trans','--s2f');
     if any(diff(n))
        error('Wb:ERR','\n   ERR unexpected local op set');
     end
     [nops,lops]=size(ops);
  else
     istr='using HAM.ops';
     nops=numel(HAM.ops); ops=QSpace(nops,1);
     for i=1:nops, ops(i)=HAM.ops(i).op; end
  end
  ops=untag(ops); 

  if isempty(k0) && isfield(HAM.user,'trotter')
     k0=getfield2(HAM.user(1).trotter,'info','ops','--def',[]);
     if ~isempty(k0)
        if iscell(k0) && numel(k0)==2
             k0=k0{2}; if ~isnumber(k0), k0=[]; end
        else k0=[]; end
     end
  end
  if isempty(k0), k0=floor((L+1)/2); end

  rop=zeros(1,nops);
  for i=1:nops, rop(i)=rank(ops(i)); end

  isf=reshape(isFermOp(HAM,ops),1,[]);
  if any(isf), Zop=HAM.oez(2).op; end

  isc=zeros(1,nops);
  for i=1:nops, if isscalarop(ops(i),'-l'), isc(i)=1; end, end
  iloc=find(isc);

  kmax=min(L,max([k0,kc])+3);

  SL=QSpace(1,nops); XL=QSpace(1,L);
  SR=QSpace(1,nops); XR=QSpace(1,L);

  if vflag, fprintf(1,'\n'); end
  for k=kc:k0-1
     if vflag, wblog('>> ','overlap %g/%g (L=%g) \r\\',k,k0,L); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
     XL(k)=contract(Ak,'!2*',Q);
  end
  for k=kc:-1:k0
     if vflag, wblog('<< ','overlap %g/%g (L=%g) \r\\',k,k0,L); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k<kc, Q={Ak,XR(k+1)}; else Q=Ak; end
     XR(k)=contract(Ak,'!1*',Q);
  end
  if vflag, fprintf(1,'\n'); end

  k=k0;
  [Ak,AA]=Load_DMRG_AK(HAM,k,AA);

  if k>kc, xAk=contract(XL(k-1),Ak); else xAk=Ak; end
  if k<kc, Akx=contract(Ak,XR(k+1)); else Akx=Ak; end
  for i=1:nops
     SL(i)=contract(Ak,'!2*',{xAk,ops(i),'-op:^s'});
     q=ops(i); if isf(i), q=Zop*q; end
     SR(i)=contract(Ak,'!1*',{Akx,q,     '-op:^s'});
  end

  for k=[ k0:L, k0-1:-1:1]
     Ak=Load_DMRG_AK(HAM,k,AA);

     if k>k0
        if vflag, wblog('>> ','1/2 calculate correlator %g/%g \r\\',k,L); end
        xAk=Ak; if k>kc, xAk=contract(XL(k-1),xAk); end

        for j=1:lops
        for i=1:nops
           q=ops(i,j); if isf(i), q=Zop*q; end
           Q={SL(i),{Ak,q,'-op:^s','*'}}; if k<kc, Q={Q,XR(k+1)}; end
           x=contract(Ak,'*',Q);
           ss(k,i,j)=getscalar(x);
        end, end

        for i=1:nops
           if isf(i), Q={Ak,Zop,'-op:^s'}; else Q=Ak; end
           SL(i)=contract(Ak,'!2*',{SL(i),Q});
        end
     elseif k<k0
        if vflag
           if k==k0-1, fprintf(1,'\n'); end
           wblog('<< ','2/2 calculate correlator %g/%g \r\\',k,L);
        end

        Akx=Ak; if k<kc, Akx=contract(Akx,XR(k+1)); end

        for j=1:lops
        for i=1:nops
           Q={{Ak,ops(i,j),'-op:^s','*'},SR(i)}; if k>kc, Q={XL(k-1),Q}; end
           x=contract(Ak,'*',Q);
           ss(k,i,j)=getscalar(x);
        end, end

        for i=1:nops
           if isf(i), Q={Ak,Zop,'-op:^s'}; else Q=Ak; end
           SR(i)=contract(Ak,'!1*',{Q,SR(i)});
        end
     else
        ot3=['-op:',Ak.info.itags{3}];
        for i=1:nops,         q0=setitags(ops(i,1),ot3); qj=q0;
        for j=1:lops, if j>1, qj=setitags(ops(i,j),ot3); end
           Q={xAk,{qj,'!2*',q0},'*','-op:^s'}; if k<kc, Q={Q,XR(k+1)}; end
           x=contract(Ak,'*',Q);
           ss(k,i,j)=getscalar(x);
        end, end
     end
     if k>=k0
          if k>=kc, XL(k)=contract(Ak,'!2*',xAk); end
     else if k<=kc, XR(k)=contract(Ak,'!1*',Akx); end
     end

     if ~isempty(iloc)
        if k>=k0
             Q=xAk; if k<kc, Q={Q,XR(k+1)}; end
        else Q=Akx; if k>kc, Q={XL(k-1),Q}; end
        end

        for j=1:lops
        for i=1:numel(iloc)
           q=fixScalarOp(ops(iloc(i),j));
           x=contract(Ak,'*',{Q,q,'*','-op:^s'});
           ss(k,nops+i,j)=getscalar(x);
        end, end
     end
  end
  if vflag, fprintf(1,'\n\n'); end

  Iout=add2struct('-',istr,ops,isf,isc,iloc,kc,k0);
  if kflag, wbstop, end

end

