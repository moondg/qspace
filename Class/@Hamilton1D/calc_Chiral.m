function [chi,Iout]=calc_Chiral(HAM,varargin)
% function [chi,Iout]=calc_Chiral(HAM [,k0])
%
%    calculate static chiral correlation function vs. distance
%    relative to (center) site k0.
% 
% Options
% 
%    -k       stop with keyboard at end of routine
%    -v       verbose mode (report progress)
%   --full    us full 3-Hubbard site permtuation operator 
%             rather than chiral operator based on spin operator
%             (the latter is the default)
% Wb,Mar04,21

% adapted from calc_ScorS
% mac>> lma DMRG_TLHubbard.mat; [cc,Ich]=calc_Chiral(HAM,'-k','-v');

  p=HAM.info.param;
  if ~isfield(p,'L') || ~isfield(p,'W')
     wbdie('invalid usage (got missing fields W or L in HAM.info.param'); 
  end

  L=p.L; W=p.W; N=L*W;
  if W<2 || N~=numel(HAM.mpo),
     wbdie('invalid system size YC%g x %g = %g !?',W,L,N);
  end

  getopt('init',varargin);
     kflag=getopt('-k');

     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end

     wchi=getopt('wchi',{});
     if isempty(wchi)
        if getopt('--full'), wchi='full';
        else wchi='spin'; end
     end

  k0=getopt('get_last',[]); och={};

  if iscell(wchi), och={'wchi',wchi}; else
     switch wchi
        case 'spin'
           n=numel(HAM.info.IS);
           if n>1, wbdie('invalid usage [got numel(IS)=%g]',n); end
           och={'wchi',{wchi, HAM.info.IS.sym}};
        case 'full'
        otherwise wchi, wbdie('invalid switch'); 
     end
  end

  [X3,I3]=get_Perm3(HAM.oez(1).op,'--chi-real',och{:});

  otags={'oc1','oc3'};
  ops=QSpace(1,2);
    ops(1)=I3.Ex;
    ops(2)=(2^3) * setitags(I3.chi,{'','',otags{:}});

  Chi={ { ops(1), ['-op:^s:' otags{1}],'*'}
        { ops(2),  '-Op:^s'               }
        { ops(1), ['-op:^s:' otags{2}]    }
  };

  if numel(k0)<=1
     [chi,Ix]=calc_Chiral_0(HAM,Chi,k0,vflag,kflag);
  elseif numel(k0)==2
     [chi,Ix]=calc_Chiral_2(HAM,Chi,k0,vflag,kflag);
  elseif isequal(k0,'--sym')
     [chi,Ix]=calc_Chiral_sym(HAM,Chi,vflag,kflag);
  elseif isequal(k0,'--plain')
     [chi,Ix]=calc_Chiral_plain(HAM,Chi,vflag,kflag);
  elseif iscell(k0) && isequal(k0{1},'--ysys') && numel(k0)==3
     k0=[k0{2:end}];
     [chi,Ix]=calc_Chiral_ysys(HAM,Chi,k0,vflag,kflag);
  elseif isequal(k0,'--yxops')
     [chi,Ix]=calc_Chiral_yxops(HAM,vflag,kflag);
  else wbdie('invalid usage (not yet implemented)');
  end

  if I3.real_chi && (~ischar(k0) || isempty(regexp(k0,'yxops')))
     chi=-chi;
  end

  Iout=add2struct('-',X3,I3,Ix,ops,L,W,N);

end

% -------------------------------------------------------------------- %
function a=area_triangle(HAM,k3)

  n=size(k3,1); a=nan(n,2);
  if size(k3,2)~=4, wbdie('invalid usage'); end

  for i=find(all(k3,2))', xy=HAM.info.XY(floor(k3(i,:)),:);
     a(i,:)=[
      + det([ xy(4,:)-xy(1,:); xy(2,:)-xy(1,:) ])
      - det([ xy(4,:)-xy(1,:); xy(3,:)-xy(1,:) ])
     ];
  end
end

% -------------------------------------------------------------------- %
function [ch1,Ix]=calc_Chiral_plain(HAM,Chi,vflag,kflag)

  p=HAM.info.param; L=p.L; W=p.W; N=L*W;
  wblog('==>','plain chiral correlations (YC%gx%g)',N,W,L);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0); ic='!1*';

  for k=kc:-1:1
     if vflag, wblog('<L|','loading A(%g/%g; %g) \r\\',k,kc,N); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);

     if k<kc, Q={XR(k+1),Ak}; else Q=Ak; end
     XR(k)=contract(Ak,ic,Q);
  end

  k1m=fix(W/2); ic='!2*'; k1=-1; l=fix(L/2);
  ch1=zeros(L,2); ch1str='';

  for k=1:N-1
     if vflag
        wblog('|L>','loading A(%g/%g) %-24s  \r\\',k,N,ch1str);
     end

     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k>=kc
        if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
        XL(k)=contract(Ak,ic,Q);
     end

     [j,i]=ind2sub([W L],k);
     if j==k1m
            k1=k; l=i;
     elseif k1<1, continue; end

     if k==k1
        Q=Ak; if k>kc, Q={XL(k-1),Q}; end
        SL(l,1)=contract(Ak,ic,{Q, Chi{1}{:} });
        SL(l,2)=contract(Ak,ic,{Q, Chi{3}{:} });            k3(l,1)=k;
     elseif k==k1+1
        SL(l,1)=contract(Ak,ic,{{Ak,SL(l,1)}, Chi{2}{:} }); k3(l,2)=k; 
        SL(l,2)=contract(Ak,ic,{ Ak,SL(l,2)}); 
     else
        for j=1:2, SL(l,j)=contract(Ak,ic,{Ak,SL(l,j)}); end
     end

     i=l-1; if i<1, continue; end
     if k==k1
        SL(i,1)=contract(Ak,ic,{ Ak,SL(i,1)});
        SL(i,2)=contract(Ak,ic,{{Ak,SL(i,2)}, Chi{2}{:} }); k3(i,3)=k;
     elseif k==k1+1
        SL(i,1)=contract(Ak,ic,{{Ak,SL(i,1)}, Chi{3}{:} });
        SL(i,2)=contract(Ak,ic,{{Ak,SL(i,2)}, Chi{1}{:} }); k3(i,4)=k;

        for j=1:2
           if k<kc
                ch1(i,j)=getscalar(contract(SL(i,j),XR(k+1)));
           else ch1(i,j)=trace(SL(i,j)); end
           ch1str=sprintf(' ch1(%g)=[ %8.3g %8.3g ]',k,imag(ch1(i,:)));
        end
     end
  end
  if vflag, fprintf(1,'\n'); end

  Ix=add2struct('-',kc,k1m,L,W,N);
  if kflag, wbstop, end

end

% -------------------------------------------------------------------- %
function [chi,Ix]=calc_Chiral_0(HAM,Chi,k0,vflag,kflag)

  p=HAM.info.param; L=p.L; W=p.W; N=L*W;

  if isempty(k0), k0=floor(N/2);
     if mod(L,2)==0
        k0=k0-W+floor(W/2);
     end
  end
  wblog('==>','chiral correlations rel. to k0=%g/%g (YC%gx%g)',k0,N,W,L);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);

  for k=1:k0
     if vflag, fprintf(1,'--> loading A(%g/%g; %g) ... \r',k,k0,N); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);

     if k>=kc
        if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
        XL(k)=contract(Ak,'!2*',Q);
     end
  end

  k2=k0+W+1; k2m=mod(k2,W);
  QR=QSpace(1,2); ic='!1*'; chi=zeros(L,1);

  for k=max(k2,kc):-1:1
     if vflag, wblog('<L|','loading A(%g/%g; %g) \r\\',k,k0,N); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);

     if k>k2 && k<=kc 
        if k<kc, Q={XR(k+1),Ak}; else Q=Ak; end
        XR(k)=contract(Ak,ic,Q);
     end
     if k>k2, continue; end

     j=mod(k,W);
     if j==k2m, k2=k;
        l=1+floor((k-1)/W);
        i1=1+mod(l,2); i2=3-i1;
     end

     if k==k2
        if k>k0
             Q=Ak; if k<kc, Q={Q,XR(k+1)}; end
        else Q={Ak,SR}; end
        QR(i1)=contract(Ak,ic,{Q, Chi{1}{:} });
     elseif k==k2-1
        QR(i1)=contract(Ak,ic,{{Ak,QR(i1)}, Chi{2}{:} });
     else
        QR(i1)=contract(Ak,ic,{Ak,QR(i1)});
     end

     if k==k2
        QR(i2)=contract(Ak,ic,{Ak,QR(i2)});
     elseif k==k0
        QR(i2)=contract(Ak,ic,{QR(i2),{Ak, Chi{3}{:} }});
     elseif k==k2-1 && k<k0-W
        Q=Ak; if k>kc, Q={XL(k-1),Q}; end
        QR(i2)=contract(Ak,'*',{QR(i2),{Q, Chi{3}{:} }});
        chi(l)=getscalar(QR(i2));
     end

     if k<k0, SR=contract(Ak,ic,{Ak,SR});
     elseif k==k0, SR=QR(i2); end
  end

  k1=k0-W; k1m=mod(k0,W); k2=k0+2*W; ic='!2*';
  QL=QSpace(1,2); SL=QSpace;

  for k=k1:N
     if vflag, wblog('|R>','loading A(%g/%g; %g) \r\\',k,k0,N); end
     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k<=k2, AA(k)=Ak; end

     j=mod(k,W);
     if j==k1m, k1=k;
        l=1+floor((k-1)/W);
        i1=1+mod(l,2); i2=3-i1;
     end

     if k==k1
        if k<k0
             Q=Ak; if k>kc, Q={XL(k-1),Q}; end
        else Q={SL,Ak}; end
        QL(i1)=contract(Ak,ic,{Q, Chi{3}{:} });
     else
        QL(i1)=contract(Ak,ic,{Ak,QL(i1)});
     end

     if k==k1
        QL(i2)=contract(Ak,ic,{{Ak,QL(i2)}, Chi{2}{:} });
     elseif k==k0+1
        QL(i2)=contract(Ak,ic,{QL(i2),{Ak, Chi{1}{:} }});
     elseif k==k1+1 && k>k0+W
        Q=Ak; if k<kc, Q={Q,XR(k+1)}; end
        QL(i2)=contract(Ak,'*',{QL(i2),{Q, Chi{1}{:} }});
        chi(l)=getscalar(QL(i2));
     end

     if k>k0+1, SL=contract(Ak,ic,{Ak,SL});
     elseif k==k0+1, SL=QL(i2); end
  end

  k1=k0-W; opsL=Chi; opsR=Chi;
  for i=1:3
     opsL{i}{2}=regexprep(opsL{i}{2},':oc',':ol');
     opsR{i}{2}=regexprep(opsR{i}{2},':oc',':or');
     if i==2
        opsL{i}{1}=itagrep(opsL{i}{1},'^oc','ol');
        opsR{i}{1}=itagrep(opsR{i}{1},'^oc','or');
     end
  end

  chi0=zeros(2,2,2);

  for x=0:1
     for l=1:2, for r=1:2
        ops=cell(2,3*W); m2=zeros(1,3*W);
        if l==1
             i=[1,2,1+W  ]; ops(1,i)=opsL; m2(i)=1;
        else i=[2+W,1+W,1]; ops(1,i)=opsL; m2(i)=1; end

        if r==1
             i=[1,2,1+W  ]+x*W; ops(2,i)=opsR; m2(i)=m2(i)+2;
        else i=[2+W,1+W,1]+x*W; ops(2,i)=opsR; m2(i)=m2(i)+2; end

        n=find(m2,1,'last');
        for i=n:-1:1, k=k1+i-1; Ak=AA(k);
           if i<n, Q={Ak,SR};
           elseif k<kc, Q={Ak,XR(k+1)}; else Q=Ak; end

           for j=1:2, if ~isempty(ops{j,i})
               Q={ Q, ops{j,i}{:} }; end
           end

           if i>1, SR=contract(Ak,'!1*',Q);
           elseif k>kc
                SR=QSpace(contractQS({XL(k-1),Ak},'*',Q));
           else SR=contract(Ak,'*',Q);
           end
        end
        chi0(x+1,l,r)=getscalar(SR);
     end, end
  end
  if vflag, fprintf(1,'\n'); end

  Ix=add2struct('-',k0,kc,k);
  Ix.chi0=chi0(:,:);

  if kflag, wbstop, end

end

% -------------------------------------------------------------------- %
function [chi,Ix]=calc_Chiral_2(HAM,Chi,k0,vflag,kflag)

  p=HAM.info.param;
  L=p.L; W=p.W; N=L*W;

  if abs(diff(k0))<2*W, wbdie(...
    'invalid usage (k0=[%g %g] values too close; W=%g)',k0,W); end
  if any(k0<=1 | k0>=N), wbdie(...
    'invalid usage (k0=[%g %g] / %g out of bounds)',k0,L); end
  wblog('==>','chiral correlations for k=[%g %g]/%g (YC%gx%g)',k0,N,W,L);

  k1=min(k0); k2=max(k0); k2m=k2+W+1;
  k3=repmat(-1,4,2); l=1; XL=QSpace(1,k1); ic='!2*';

  SL=QSpace(2,1); m=1;

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);

  for k=1:N
     if vflag, wblog('-->','loading A(%g/%g; %g) ... \r\\',k,k2m,N); end

     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k>=kc
        if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
        XL(k)=contract(Ak,ic,Q);
     end

     if k<k1, continue; end
     if k>=k2m && kc<N, ic='*'; end

     if k==k2, XC=SL; m=2; end

     for j=1:m
        if any(k==[k1,k2])
           if    k>=k2, Q={XC(j),Ak};
           elseif k>kc, Q={XL(k-1),Ak}; else Q=Ak; end

           SL(1,j)=contract(Ak,ic,{Q,Chi{1}{:} });
           SL(2,j)=contract(Ak,ic,{Q,Chi{3}{:} });
           if j==1, k3(l)=k+.1; l=l+1; end

        elseif any(k==[k1,k2]+1)
           SL(1,j)=contract(Ak,ic,{{Ak,SL(1,j)}, Chi{2}{:} });
           SL(2,j)=contract(Ak,ic,{ Ak,SL(2,j) });
           if j==1, k3(l)=k+.21; l=l+1; end

        elseif any(k==[k1,k2]+W)
           SL(1,j)=contract(Ak,ic,{ Ak,SL(1,j) });
           SL(2,j)=contract(Ak,ic,{{Ak,SL(2,j)}, Chi{2}{:} });
           if j==1, k3(l)=k+.22; l=l+1; end
        elseif any(k==[k1,k2]+W+1)
           SL(1,j)=contract(Ak,ic,{SL(1,j),{Ak, Chi{3}{:} }});
           SL(2,j)=contract(Ak,ic,{SL(2,j),{Ak, Chi{1}{:} }});
           if j==1, k3(l)=k+.3; l=l+1; end
        else
           for i=1:2
           SL(i,j)=contract(Ak,ic,{Ak,SL(i,j)}); end
        end
     end

     if k>=k2m && kc<N, break; end
  end
  if vflag, fprintf(1,'\n'); end

  chi=getscalar(SL);

  Ix=add2struct('-',k1,k2,kc,k,k3);

  Ix.trXL=trace(XL(k));
     e=norm(Ix.trXL-1);
     if e>1E-10, wblog('WRN','unexpected norm 1%+.2g',Ix.trXL-1); end

  Ix.atri=area_triangle(HAM,k3');
  a=Ix.atri(:);
  if any(~isfinite(a) | abs(a-1)>1E-3)
     wblog('WRN','got triangles with unexpected area');
     disp(Ix.atri);
  end

  if kflag, wbstop, end
end

% -------------------------------------------------------------------- %
% compute chiral corerlations symmetrically out from center
% -------------------------------------------------------------------- %
% WRN! rather expensive since the individual chiral correlations
% barely have any overlap, hence are effectively computed independent
% of each other (other than loading AK once for all of them)

function [chi,Ix]=calc_Chiral_sym(HAM,Chi,vflag,kflag)

  p=HAM.info.param; L=p.L; W=p.W; N=L*W;

  k0=floor(N/2);
     if mod(L,2)==0
        k0=k0-W+floor(W/2);
     end
  wblog('==>','sym. chiral correlations around k0=%g/%g (YC%gx%g)',k0,N,W,L);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);

  k1m=mod(k0,W); ic='!2*'; k1=-1; l=fix(L/2);

  SL=QSpace(l,2); XL=QSpace(1,N);
  SR=QSpace(l,2); XR=QSpace(1,N);
  k3=zeros(L,4);

  for k=1:k0-1
     if vflag, wblog('|L>','loading A(%g/%g; %g) \r\\',k,k0,N); end

     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k>=kc
        if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
        XL(k)=contract(Ak,ic,Q);
     end

     [j,i]=ind2sub([W L],k);
     if j==k1m
            k1=k; l=i;
     elseif k1<1, continue; end

     if k==k1
        Q=Ak; if k>kc, Q={XL(k-1),Q}; end
        SL(l,1)=contract(Ak,ic,{Q, Chi{1}{:} });
        SL(l,2)=contract(Ak,ic,{Q, Chi{3}{:} });            k3(l,1)=k;
     elseif k==k1+1
        SL(l,1)=contract(Ak,ic,{{Ak,SL(l,1)}, Chi{2}{:} }); k3(l,2)=k; 
        SL(l,2)=contract(Ak,ic,{ Ak,SL(l,2)}); 
     else
        for j=1:2, SL(l,j)=contract(Ak,ic,{Ak,SL(l,j)}); end
     end

     i=l-1; if i<1, continue; end
     if k==k1
        SL(i,1)=contract(Ak,ic,{ Ak,SL(i,1)});
        SL(i,2)=contract(Ak,ic,{{Ak,SL(i,2)}, Chi{2}{:} }); k3(i,3)=k;
     elseif k==k1+1
        SL(i,1)=contract(Ak,ic,{{Ak,SL(i,1)}, Chi{3}{:} });
        SL(i,2)=contract(Ak,ic,{{Ak,SL(i,2)}, Chi{1}{:} }); k3(i,4)=k;
     else
        for j=1:2, SL(i,j)=contract(Ak,ic,{Ak,SL(i,j)}); end
     end

     for i=1:l-2
        for j=1:2, SL(i,j)=contract(Ak,ic,{Ak,SL(i,j)}); end
     end
  end

  k2m=mod(k0+1,W); k2=N+1; ic='!1*';

  for k=N:-1:k0
     if vflag, wblog('<R|','loading A(%g/%g; %g) \r\\',k,k0,N); end

     [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
     if k<=kc
        if k<kc, Q={Ak,XR(k+1)}; else Q=Ak; end
        XR(k)=contract(Ak,ic,Q);
     end

     [j,i]=ind2sub([W L],k);
     if j==k2m
            k2=k; l=i;
     elseif k2>N, continue; end

     if k==k2
        Q=Ak; if k<kc, Q={Q,XR(k+1),Q}; end
        SR(l,1)=contract(Ak,ic,{Q, Chi{1}{:} });
        SR(l,2)=contract(Ak,ic,{Q, Chi{3}{:} });            k3(l,1)=k;
     elseif k==k2-1
        SR(l,1)=contract(Ak,ic,{{Ak,SR(l,1)}, Chi{2}{:} }); k3(l,2)=k;
        SR(l,2)=contract(Ak,ic,{ Ak,SR(l,2)});
     else
        for j=1:2, SR(l,j)=contract(Ak,ic,{Ak,SR(l,j)}); end
     end

     i=l+1; if i>L, continue; end
     if k==k2
        SR(i,1)=contract(Ak,ic,{ Ak, SR(i,1) });
        SR(i,2)=contract(Ak,ic,{{Ak,SR(i,2)}, Chi{2}{:} }); k3(i,3)=k;
     elseif k==k2-1,
        SR(i,1)=contract(Ak,ic,{{Ak,SR(i,1)}, Chi{3}{:} });
        SR(i,2)=contract(Ak,ic,{{Ak,SR(i,2)}, Chi{1}{:} }); k3(i,4)=k;
     else
        for j=1:2, SR(i,j)=contract(Ak,ic,{Ak,SR(i,j)}); end
     end

     for i=l+2:L
        for j=1:2, SR(i,j)=contract(Ak,ic,{Ak,SR(i,j)}); end
     end
  end
  if vflag
     s=whos('SL','SR','XL','XR','Ak');
     wblog(' * ','done @ %s',num2str2(sum([s.bytes]),'--bytes'));
  end

  n=fix(k0/W)-1; chi=zeros(L,2,2);
  for k=1:n, for l=1:2
     for i=1:2, for j=1:2
        k_=L-k+(2-l);
        x=(k_-k-1);
        chi(x,i,j)=getscalar(contract(SL(k,i),SR(k_,j)));
     end, end
  end, end
  chi=chi(:,:);

  if ~isreal(chi)
     for k=1:n, for i=1:2
        if kc<k0
             ch1(k,i)=trace(SL(k,i));
        else ch1(k,i)=getscalar(contract(SL(k,i),XR(k0))); end

        k_=L-k+1;
        if kc>k0
             ch1(k_,i)=trace(SR(k,i));
        else ch1(k_,i)=getscalar(contract(SR(k_,i),XL(k0-1))); end
     end, end
  else ch1=[];
  end

  Ix=add2struct('-',k1m,k2m,kc,k0,k3,ch1);
  Ix.atri=area_triangle(HAM,k3);

  if kflag, wbstop, end
end

% -------------------------------------------------------------------- %
function [chi,Ix]=calc_Chiral_ysys(HAM,Chi,k2m,vflag,kflag)

  p=HAM.info.param; L=p.L; W=p.W; N=L*W;

  if numel(k2m)~=2 || norm(k2m-round(k2m)) || ...
     any(k2m)>=W || diff(k2m)>=W-1 || diff(k2m)<2
     wbdie('invalid usage (invalid k2=[%g %g] / %g)',k2m,W);
  end
  wblog('==>','chiral correlations for k=[%g %g] (YC%gx%g)',k2m,W,L);

  XL=QSpace(1,N); ic='!2*'; chi=zeros(1,L-1);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);
  if kc>k2m(end), wblog('ERR',...
    'state not RL orthonormalized (kc=%g/%g)',kc,k2m(end)); end

  Chi_=Chi;
  Chi_{1}{2}='-op:^s:oc_1';
  Chi_{2}{1}=itagrep(Chi{2}{1},'oc','oc_');
  Chi_{3}{2}='-op:^s:oc_3';

  for ix=1:L-1
     for iy=1:W, k=iy+(ix-1)*W;
        if vflag, wblog('-->',...
           'loading A(%g: %g/%g, %g) ... \r\\',k,iy,W,ix); end
        [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
        if k>=kc
           if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
           XL(k)=contract(Ak,ic,Q);
        end
     end

     k1=k2m(1) + (ix-1)*W;
     k2=k2m(2) + (ix-1)*W;

     for k=k1:k1+W-1
        [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
        if k==k1;
           if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
           SL=contract(Ak,ic,{Q,Chi{1}{:}});
        elseif k==k1+1;
           SL=contract(Ak,ic,{{Ak,SL},Chi{2}{:}});
        elseif k<k2
           SL=contract(Ak,ic,{Ak,SL});
        elseif k==k2;
           SL=contract(Ak,ic,{{Ak,SL},Chi_{1}{:}});
        elseif k==k2+1;
           SL=contract(Ak,ic,{{Ak,SL},Chi_{2}{:}});
        else
           SL=contract(Ak,ic,{Ak,SL});
        end
     end

     k1=k1+W; k2=k2+W;

     for k=k1:k2
        [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
        if k==k1
           SL=contract(Ak,ic,{{Ak,SL},Chi{3}{:}});
        elseif k<k2
           SL=contract(Ak,ic,{Ak,SL});
        else
           if k<kc, Q={XR(k+1),Ak}; else Q=Ak; end
           SL=contract(Ak,'*',{{Q,SL},Chi_{3}{:}});
        end
     end
     chi(ix)=getscalar(SL);
  end
  if vflag, fprintf(1,'\n'); end

  Ix=add2struct('-',k2m,kc);

  if kflag, wbstop, end
end

% -------------------------------------------------------------------- %
function [ss,Ix]=calc_Chiral_yxops(HAM,vflag,kflag)

  p=HAM.info.param; L=p.L; W=p.W; N=L*W;

  if ~isfield(HAM.info,'xops')
     wblog('WRN','missing operators HAM.info.xops'); 
     chi=[]; Ix=[]; return
  else
     xops=HAM.info.xops; nops=numel(xops);
     wblog('==>','yxops correlations (YC%gx%g; %g xops)',W,L,nops);
  end

  XL=QSpace(1,N); ic='!2*'; ss=zeros(L,W,nops);

  [kc,Ic]=Load_DMRG_AK(HAM); AA=QSpace(1,0);
  if kc>W/2, wblog('ERR',...
    'state not RL orthonormalized (kc=%g/%g)',kc,W); end

  yidx=zeros(W,2);

  isf=isFermOp(HAM,xops);
  if any(isf), Zop=HAM.oez(2).op; end

  for ix=1:L
     for iy=1:W, k=iy+(ix-1)*W;
        if vflag, wblog('-->',...
           'loading A(%g: %g/%g, %g) ... \r\\',k,iy,W,ix); end
        [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
        if k>=kc
           if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
           XL(k)=contract(Ak,ic,Q);
        end
     end

     for iy=1:W
        k1=max(1,floor((W-(iy-1))/2));
        k2=k1+(iy-1);
        yidx(iy,:)=[k1 k2];

        k1=k1+(ix-1)*W;
        k2=k2+(ix-1)*W;

        for iop=1:nops, x=xops(iop);
           for k=k1:k2
              [Ak,AA]=Load_DMRG_AK(HAM,k,AA);
              if k1==k2    % X'.X on the same site
                 if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
                 if numel(x.Q)>2
                      x2=contractQS(x,'13*',x,'13');
                 else x2=contractQS(x,'1*', x,'1' );
                 end
                 SL=contract(Ak,'*',{Q,x2,'-op:^s'});
              elseif k==k1
                 if k>kc, Q={XL(k-1),Ak}; else Q=Ak; end
                 SL=contract(Ak,ic,{Q,x,'-op:^s','*'});
              elseif k<k2
                 if ~isf(iop), Q=Ak; else Q={Ak,Zop,'-op:^s'}; end
                 SL=contract(Ak,ic,{Q,SL});
              else
                 if isf(iop), x=contract(Zop,2,x,1); end
                 SL=contract(Ak,'*',{{Ak,SL},x,'-op:^s'});
              end
           end

           ss(ix,iy,iop)=getscalar(SL);
        end
     end
  end
  if vflag, fprintf(1,'\n'); end

  Ix=add2struct('-',kc,yidx,xops,isf);

  if kflag, wbstop, end
end

% -------------------------------------------------------------------- %

