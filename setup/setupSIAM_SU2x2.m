% setup operator set for standard SIAM problem (single channel)
% using particle-hole and spin SU(2) symmetry (if B=0)
% Wb,Sep19,11
% -------------------------------------------------------------------- %

% adapted from setupSIAM.m
  global param

  setdef('U',0.08,'Gamma',0.01,'B',0,'Lambda',2,'N',50);
  if isvar('Uepsd'), q=-U/Uepsd; else q=-U/2; end
  setdef('epsd',q);

  Bflag=(numel(B)==1 && B~=0);

  setdef('SUNflag',0);
  if Bflag && SUNflag, wbdie(...
    'invalid usage (cannot use SUNflag together with B=%g',B);
  end

  istr=''; sym='';
  param=add2struct('-',istr,sym,U,epsd,Gamma,B,Lambda,N);

  co={'-x'};
  if exist('z','var'), co(end+1:end+2)={'z',z}; param.z=z; end
  if exist('ALambda','var') && ~ALambda
     co(end+1)={'-w'};
     if isfield(param,'ALambda') param=rmfield(param,'ALambda'); end
  end

  clear ZFLAG*

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

  setdef('NC',1,'SYM',''); Aflag=0;

  os={'FermionS','','',{'NC',NC,'-v'}};
  if NC>1, os(2:3)={',SUNchannel', sprintf(',SU%g',NC)}; end

  if SUNflag
     os{1}=os{1}(1:end-1);
     if NC~=0.5
          os={'Fermion','SUNchannel',sprintf(',SU%g',2*NC),{'NC',2*NC,'-v'}};
     else os={'Fermion','Acharge','A',{'NC',1,'-v'}};
     end

     [FF,Z,IS]=getLocalSpace(os{1},os{2}, os{4}{:});

     istr=sprintf('spinless SU%g Anderson model (SIAM)',2*NC);

  elseif NC==1 && isequal(SYM,'Acharge') ...
     || isequal(SYM,'A') && isequal(IS.sym,'Acharge')

     [FF,Z,SS,IS]=getLocalSpace(os{1},'Acharge', os{4}{:});
     Aflag=2;

  elseif epsd~=-U/2 && Bflag || isequal(SYM,['A,A' os{3}])
     [FF,Z,SS,IS]=getLocalSpace(os{1},['Aspin,Acharge' os{2}], os{4}{:});
     Aflag=1;

  end

   if exist('TK0')==1, TK=TK0;
   elseif NC==1, TK=TKondo;
   end

  if isempty(istr)
     if NC==1
          istr='one-channel Anderson model (SIAM)';
     else istr=sprintf('%g-channel Anderson model (SIAM)',NC); end
  end
  param.istr=istr;

  if Aflag, NF=FF;
     for i=1:numel(FF)
        FF(i)=fixScalarOp('-l',FF(i));
        NF(i)=FF(i)'*FF(i);
     end
     Sop=SS(end);

     HU=skipzeros(NF(1)*NF(2));
     H0=skipzeros(epsd*(NF(1)+NF(2)) + U*HU - B*SS(3));
     E0=IS.E;

  elseif SUNflag

     NF=QSpace(contractQS(FF,'13*',FF,'13'));
     HU=0.5*skipzeros(NF*NF-NF);
     H0=skipzeros(epsd*NF + U*HU);

     E0=IS.E;

     V=getIdentity(E0,'-0');
     X=getIdentity(E0,V,2); X=getsub(X,matchIndex(X.Q{3},[1 1]));
     Sop=contract(X,2,V,'2*',[1 3 2]) / sqrt(2);

     [~,is]=sortrows(Sop.Q{1},'descend');
     Sop=getsub(Sop,is);
     Sop.info.otype='operator'; Sop_=Sop; n=numel(Sop.data);

     for i=n:-1:1
        Sop(i)=getsub(Sop(1),i);
        X=getLocalSpace('Spin',Sop(i).info.cgr(1).type,Sop(i).Q{1}(1,:));
        Sop(i).data{1}=X.data{1};
        Sop(i)=Sop(i) * (norm(X)/norm(Sop(i)));
     end

     clear A V X

  elseif epsd~=-U/2 && ~Bflag || isequal(SYM,['SU2,A' os{3}])
     [FF,Z,SS,IS]=getLocalSpace(os{1},['SU2spin,Acharge' os{2}], os{4}{:});
     E0=IS.E;

     NF=QSpace(contractQS(FF,'13*',FF,'13'));
     HU=0.5*skipzeros(NF*NF-NF);
     H0=skipzeros(epsd*NF + U*HU);

  else, ZFLAG=2;
     if NC>1
          os(2:3)={'SpNchannel', sprintf('Sp%g',2*NC)};
     else os(2:3)={'SU2charge', 'SU2'}; end

     if Bflag || isequal(SYM,['A,' os{3}]), if ~Bflag, Bflag=2; end
          os{2}=['Aspin,'   os{2}];
     else os{2}=['SU2spin,' os{2}]; end

     [FF,Z,SS,IS]=getLocalSpace(os{1},os{2},os{4}{:});
     E0=IS.E;

     HU=0.5*skipzeros(IS.Q2-E0);
     if Bflag
        H0=skipzeros(U*HU-B*SS(3));
     else
        H0=skipzeros(U*HU);
     end
  end

  if isempty(H0), H0=0*E0; end

  FN=comm(FF,HU);

  q=getvac(Z);
  A0=QSpace(permuteQS(getIdentityQS(q,1,Z,1),[1 3 2]));

  SYM=getsym(H0); param.sym=IS.sym;

  [ff,~,Iwc]=getNRGcoupling(Gamma,Lambda,N,co{:});

% -------------------------------------------------------------------- %
  FC=FF; op1=[];

  if Bflag
     op1=[ FC(1), FN(1) ];
     op2=[ FC([1,1]) ];    zflags=ones(1,numel(FC));
  elseif SUNflag
     op1=[ FC(1), FN(1) ];
     op2=[ FC(1), FC(1) ]; zflags=ones(1,2*numel(FC));
  else
     op1=[ FC(1), FN(1), SS(end) ];
     op2=[ FC(1), FC(1), SS(end) ]; zflags=[1 1 0];
  end

  if isset('ZFLAG') && ZFLAG>1
     ac_cmd='[ac,Ic]=getGC_bulla(ox, ax(:,1), ax(:,2) + 0.5*ax(:,1),''-iB'');';
  else
     ac_cmd='[ac,Ic]=getGC_bulla(ox, ax(:,1), ax(:,2),''-iB'');';
  end

  Z0=Z;

% -------------------------------------------------------------------- %
% add first n0add Wilson sites to H0 and A0 // Wb,Jul01,13
% -------------------------------------------------------------------- %

  if isset('n0add') && n0add>0

     ZF=FF;
     for i=1:numel(FF), ZF(i)=contractQS(Z,2,FF(i),1); end

     for it=1:n0add
        wblog(' * ','adding %g/%g Wilson site to A0 ...',it,n0add);

        if isset('ZFLAG') && ZFLAG>1
           if mod(it+n0add-1,2)==(ZFLAG-2), q1=FF; else q1=ZF; end
           q2=q1;
        else
           q1=FF; q2=ZF;
        end

        for i=1:numel(q1)
           q1(i)=contractQS(A0,[1 3],contractQS(A0,3,q1(i),2),[1 3]);
        end

        A0_=A0;
        A0=QSpace(permuteQS(getIdentityQS(A0,2,IS.E,1),[1 3 2]));

        H0_=H0;
        H0=QSpace(contractQS(A0,[1 3],contractQS(H0,2,A0,1),[1 3]));

        for i=1:numel(q2)
           q1(i)=contractQS(q1(i),1,A0,1,'-L');
           q2(i)=contractQS(A0,3, ff(it)*q2(i),1,'-L'); % f' :: LRs(o)
           if numel(q2(i).Q)==2
                q=QSpace(contractQS(q1(i),[1 3],q2(i),[1 3]));
           else q=QSpace(contractQS(q1(i),[1 2 4],q2(i),[1 4 3]));
           end

           if isset('ZFLAG') && ZFLAG>1
                H0=H0+q;
           else H0=H0+q+q'; end
        end
     end

     fprintf(1,'\n'); clear ZF q1 q2
     param.n0add=n0add;

     ff=ff(n0add+1:end);
  end

% -------------------------------------------------------------------- %

  inl(1); disp(param);

