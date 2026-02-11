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
     beta =getopt('beta',10);
     Rc=getopt('Rc',[]);
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
     k0=getfield2(HAM.user(1).trotter,'info','ops',{[]});
     if ~isempty(k0)
        if iscell(k0) && numel(k0)==2
             k0=k0{2}; if ~isnumber(k0), k0=[]; end
        else k0=[]; end
     end
  end
  if isempty(k0), k0=ceil(L/2); end

  rop=zeros(1,nops);
  for i=1:nops, rop(i)=rank(ops(i)); end

  isf=reshape(isFermOp(HAM,ops),1,[]);
  if any(isf), Zop=HAM.oez(2).op; end

  isc=zeros(1,nops);
  for i=1:nops, if isscalarop(ops(i),'-l'), isc(i)=1; end, end
  iloc=find(isc);

  kmax=min(L,max([k0,kc])+3);

  k=kc;
  [Ak,AA]=Load_DMRG_AK(HAM,k,AA);

  IPsi=-99;
  if rank(Ak)==4, Akc=Ak;
     dg=getDimQS(Ak); nPsi=dg(1,end); NPsi=dg(end);
     if NPsi>1
        if isempty(Rc)
           if beta<Inf && nPsi>1, IPsi=-1; else IPsi=1; end
        elseif isint(Rc), IPsi=Rc;
        elseif isequal(Rc,'all'), IPsi=1:nPsi;
        else IPsi=0;
           if ~isa(Rc,'QSpace'), wbdie('invalid Rc'); end
           q=normQS(Rc); if abs(q-1)>1E-12
              wblog('WRN','invalid trace(Rc) = %.4g',q); 
              Rc=Rc/normQS(Rc);
           end
        end
     end

     Ik=load_dmrg_data(HAM,k,'info');
     if isfield(Ik,'Eg'), Eg=Ik.Eg(end); 
        if isnumeric(Eg)
           if nPsi~=1, wbdie('unexpected Eg for nPsi=%d',nPsi); end
        else
           nd=numel(Eg.data);
           if all(IPsi>0), l=0;
              eps=max(abs(diag(Eg,'-d'))); eps=1e-12*max(1,eps);
              for i=1:nd
                  n=length(Eg.data{i}); j=l+1:l+n;
                  Eg.data{i}=double(single(Eg.data{i})) + j*eps; l=l+n;
              end
           end
        end
     else
        wblog('WRN','missing field HAM(%d).info.Eg',kc);
        Eg=getIdentity(Ak,4);
     end
  end

for iPsi=IPsi

  if iPsi>0 && isnumeric(Eg)
     Rc=1/sqrt(NPsi); Ak=Rc*Akc;
     AA(k)=Ak;
  elseif iPsi>=-1
     if iPsi>0
        wblog(' * ','computing correlations based on state %d/%d',iPsi,nPsi);

        e0=sort(diag(Eg,'-d')); e0=e0(iPsi);
        for i=1:nd, j=find(Eg.data{i}==e0,1);
           if ~isempty(j)
              Rc=getsub(Eg,i);
                 n=length(Eg.data{i});
                 q=zeros(n,1); q(j)=1;
              Rc.data{1}=q; break
           end
        end
        Rc=Rc/normQS(Rc);
     elseif iPsi==-1
        eg=sort(diag(Eg,'-d')); dE=max([1e-3,diff(eg)]); e0=eg(1);
        Rc=getIdentity(Eg,2);   bfac=beta/dE;
        for i=1:nd
           Rc.data{i}=exp(-bfac*(Eg.data{i}-e0))';
        end
        Rc=Rc/trace(Rc);

     end
     Ak=contract(Akc,4,Rc,1);
     AA(k)=Ak;
  end

  SL=QSpace(1,nops);
  SR=QSpace(1,nops);

  if kc<k0, XL=QSpace(1,L); XL(k)=contract(Ak,'!2*',Ak); end
  if kc>k0, XR=QSpace(1,L); XR(k)=contract(Ak,'!1*',Ak); end

  if vflag, fprintf(1,'\n'); end
  for k=kc+1:k0-1
     if vflag, wblog('>> ','overlap %g/%g (L=%g) \r\\',k,k0,L); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     XL(k)=contract(Ak,'!2*',{XL(k-1),Ak});
  end
  for k=kc-1:-1:k0
     if vflag, wblog('<< ','overlap %g/%g (L=%g) \r\\',k,k0,L); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     XR(k)=contract(Ak,'!1*',{Ak,XR(k+1)});
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
  if numel(IPsi)>1, SS{iPsi}=ss; end
end

  if vflag, fprintf(1,'\n\n'); end

  if isvar('SS') && numel(SS)~=1, ss=SS; end
  if nargout>1
     Iout=add2struct('-',istr,ops,isf,isc,iloc,Rc,kc,k0); end
  if kflag, wbstop, end

end

