function [Ht,Iout]=tuneHamQtot(H0,Qtot,varargin)
% function [Ht,Iout]=tuneHamQtot(H0,Qtot,Qop [,opts])
%    
%    Adjust `chemical potential' (dmu) in every symmetry label
%    such that Ht = H0+dmu*NQ has low-energy expectation value
%    equal or close to Qtot.
%
%     - here Qtot can be an arbitrary float
%
%     - symmetries that have NaN-entries in Qtot are not tuned,
%       i.e., ignored, abelian and non-abelian alike
%
% Wb,Aug18,18

% previously for non-abelian symmetries also Qtot==0 enries
% were ignored, i.e., set to NaN

  getopt('init',varargin);
     nrep=getopt('nrep',2);
     Ufac=getopt('Ufac',[]);
  Qop=getopt('get_last',[]);

  if numel(H0.data)==1 && isequal(H0.Q{1},0)
     if ~Qop, wbdie(...
     'invalid usage (Qop needs to be specified hqving Q=0)'); end
  end

  if ~isempty(Qop), wtune=1; else wtune=2; end

  EK0=get_EK0(H0);
  dE=mean(diff(sort(cat(1,EK0.data{:}))));
  if dE<1E-3
     if dE<=0, wbdie('invalid dE=%g',dE); 
     else wblog('WRN','got small dE=%g',dE); end
  end

  Ek=cat(1,EK0.data{:});
  Q=EK0.Q{1}; q=sum((Q-repmat(Qtot,size(Q,1),1)).^2,2);
  i=(q==min(q)); 

  if any(Ek(i)<min(Ek)+1E-6)
   % already got lowest (closest) symmetry sector correct - done!
   % WRN! the tuning below may lead to collaps to single symmetry
   % sector after truncation due to dense spectrum at low energies
   % i.e., when having dE >> level spacing since dE necessarily is a
   % rather crude estimate that does not care about model specifics
   % => avoid tuning if there is no need! // Wb,Feb02,23
     Ht=H0;
     if nargout>1, EKt=EK0; mu=[]; xd=[]; yd=[];
        Iout=add2struct('-',Qtot,wtune,Qop,EK0,EKt,dE);
        if wtune, Ufac=0;
             Iout=add2struct(Iout,Ufac);
        else Iout=add2struct(Iout,mu,xd,yd); end
        Iout.wtune=-Iout.wtune;
     end
     return
  end

  if wtune==2
     if isempty(Ufac) || iscell(Ufac), x=dE/25;
        if isempty(Ufac)
             Ufac=x;
        else Ufac=x*Ufac{1}; end
     end
  elseif ~isempty(Ufac)
     wbdie('invalid usage (got Ufac with Qop)');
  end

  r=getsym(H0,'-r'); nsym=numel(r);
  s=r; s(find(r==0))=1;

  Q=H0.Q{1}; n2=size(Q,2);
  if sum(s)~=n2, wbdie('got nq=%g/%g',n2,sum(s)); end

  if ~isequal(size(Qtot),[1,n2]), wbdie('invalid Qtot'); end
  nd=numel(H0.data); iq=cumsum([0,s]); 

  for k=1:nsym, j=iq(k)+1:iq(k+1);
     qj=Qtot(j); Qj=uniquerows(Q(:,j)); dQ=Qj-repmat(qj,size(Qj,1),1);
     x=sum(dQ.^2,2); l=find(x==min(x)); 

     if r(j)>1 && ~isfinite(norm(qj)) && any(~isnan(qj))
        wbdie('invalid usage (got partial NaN for non-abelian)');
     elseif wtune==1 && r(j), wbdie(['invalid usage ' ... 
        '(Qop may only be specified for abelian symmetry)']);
     end
  end

  if all(isnan(Qtot)), wbdie([
     'invalid usage (0/%g symmetries to tune)\n' ...
     'hint: use flag -f or Qtot~=0'],nsym); 
  end

  if wtune==2
   % NB! this case acts differently as compared to when Qop is
   % specified (see below) where the latter case acts like shifting
   % the chemical potential in a linear fashion assuming U1 symmetry.
   % Here by building Qop from scratch, this can only make use of the
   % total symmetry labels for given block as found in H0.Q.
   % Here for the case no Qop: H0 -> H0 + HU with HU = U*Id
   %  * Qop is initialized to diagonal matrices
   %    with U = Ufac*norm(symmetry labels - Qtot)^2
   %  * by using norm2(), this acts like an interaction (Casimir) that
   %    favors Qtot, like Hubbard (U/2)*(\hat{n}-nd)^2 to favor filling nd.
   % => the correction added to H0 here with no Qop this needs to be
   %    substracted in the caller (NRG) routine, prior to reentering
   %    here for the next iteration, since otherwise iteraction at
   %    earlier iterations gets included multiplet times!
   %    (this is in contrast to no Qop, in which case only the local
   %    chemical potential at a given iteration/site is tuned).
   % #UNDO_TUNE_QOP // Wb,Jan31,23
     Qop=repmat(getIdentity(H0,2),1,nsym);
     Q=Qop(1).Q{1};

     q=size(Ufac);
     if q(1)~=1 || q(2)>1 && q(2)~=nsym
        wbdie('invalid usage (size mismatch of Ufac (%d/%d)',q(2),nsym); 
     end

     for i=1:nd, s=size(Qop(1).data{i});
        for k=1:nsym, j=iq(k)+1:iq(k+1);
           U=Ufac(min(end,k))*sum((Q(i,j)-Qtot(j)).^2);
           Qop(k).data{i}=diag(repmat(U,1,s(1)));
        end
     end

     Ht=H0;
     for k=1:numel(Qop), Ht=Ht+Qop(k); end

     if nargout>1
        EKt=get_EK0(Ht);
        for k=1:nsym, Qop(k)=diag(Qop(k),'-c'); end
        Iout=add2struct('-',Qtot,wtune,Qop,EK0,EKt,dE,Ufac);
     end
     return
  end

  nQ=numel(Qop);

  if ~isQSpace(Qop), wbdie('invalid usage (invalid QSpace Qop)'); end
  for k=1:nQ, Q=Qop(k).Q; l=numel(Q);
     if l~=2, wbdie('invalid Qop (got rank-%d)',l);
     elseif ~isequal(Q{:}), wbdie('invalid Qop (non-scalar)'); end
  end

  nx=5; xd=linspace(-1,1,nx); yd=zeros(1,nx);

  Ht=H0;
  for irep=1:nrep
     for k=1:nQ
        j=find(sum(Qop(k).Q{1}.^2,1));
        j=j(find(~isnan(Qtot(j))));
        if isempty(j), wbdie('invalid Qop (fails to address Qtot)'); end

        for i=1:nx, [Ex,Ix]=eigQS( Ht + xd(i) * Qop(k) );
           Rx=Ix.EK; Rq=Rx; Q=Rq.Q{1}; E0=min(Ex(:,1));
           Qop(k)
           for l=1:numel(Rx.data)
              x = exp( -(Rx.data{l}(1)-E0) / dE);
              Rq.data{l} = x*sum( (Qtot(j) - Q(l,j)).^2 );
              Rx.data{l} = x;
           end
           yd(i) = trace(Rq)/trace(Rx);
        end

        p=polyfit(xd,yd,2);
        x=-p(2)/(2*p(1));
        mu(irep,k)=x;

        Ht = Ht + x*Qop(k);
     end
  end

  if nargout>1
     EKt=get_EK0(Ht);
     Iout=add2struct('-',Qtot,wtune,Qop,EK0,EKt,dE,mu,xd,yd);
  end

end

% -------------------------------------------------------------------- %

function EK=get_EK0(HK)

  [ee,Ie]=eigQS(HK);
  EK=Ie.EK; nd=numel(EK.data);
  for i=1:nd, EK.data{i}=min(EK.data{i}); end

end

% -------------------------------------------------------------------- %

