function [HAM,Iout]=sweepPsi_tmpo(HAM,varargin)
% function [HAM,Iout]=sweepPsi_tmpo(HAM,Nkeep [,opts])
%
%   Complete one full time step tau based on a unitary MPO
%   that encodes the times step within 2nd order Trotter,
%   which must be setup prior to this routine using setupTrotter2.
%
%   This routine finds best variational fit for |Psi(dt)> := U(dt)*|Psi>
%   by sandwitching the MPO in mpo(k).trotter in between
%   <Psi(dt)| and |Psi>.
%
% Options
%
%    'ntau',..  number of tau steps to perform (100)
%    'nfit',..  (max) number of var. fitting sweeps for single time step (2)
%    'Nkeep',.. (max) number of variational states to keep (256)
%    'stol',..  tolerance, i.e., cutoff on singlular values (1E-12)
%    '--olap'   compute overlap of consecutive wave functions
%
% Wb,May31,19

% fitting sweeps appear to converge really fast
% (within two sweeps down to double precision in L=32 tight-binding chain!)

  getopt('INIT',varargin);
     ntau =getopt('ntau',100); ntau=ceil(ntau);
     nft  =getopt('nfit',2);

     Nkeep=getopt('Nkeep',512);
     stol =getopt('stol',1E-6);

     if getopt('--olap'), olap=1; 
     elseif getopt('--Olap'), olap=2; else olap=0; end

     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0;
     elseif olap>1, vflag=2; else vflag=1; end

  getopt('check_error');

  if ~exist('labindex')==5
       jid=1;
  else jid=labindex;
  end

  if ntau<=2, vflag=vflag+1; olap=2; end

  if isempty(getfield2(HAM.user,'trotter','info','dt',{[]})), wbdie(...
    'invalid usage (need to setup first [using setupTrotter2]'); end

  L=length(HAM);
  AF=QSpace(1,L); Sxt=[]; xft=[];

  e1=HAM.oez(1).op; oo={'Nkeep',Nkeep,'stol',stol};

  Itr=HAM.user.trotter.info;

  dt=Itr.dt;
  fop=Itr.fermionic; if fop, Zop=HAM.oez(2).op; end

  ops=getfield2(Itr,'ops',{[]}); nops=numel(ops);
  if nops
     if nops~=2 || ~iscell(ops) || ~isnumber(ops{2});
        wbdie('single operator needed for correlation function');
     end
  end

  Iout=add2struct('-','runtime=[]','tt=[]',Nkeep,stol,ntau,nft,Sxt,xft);

  if vflag>1 && olap
   % NB! all of the following converges quickly (within ~2 sweeps!)
   % where the first overlap <f|f'> becomes 1. up to 1E-12
   % and the second part <t|f'> describes the change of the wavefunction
   % within given time step; this also relates to the maximal change
   % in expectaction values later
     fprintf('%s |1-<f|f''>|  (1-<t|f''>)  | svd2tr_t  svd2tr_f %s\n',...
     repmat('-',1,29), repmat('-',1,6));
  end

  for itau=1:ntau; tic
    for ift=1:nft
       Il=load_trotter_data(HAM,1); AF(1)=Il.Af;
       if ift==1
          d=getqdir(Il.At); if any(d(1:3)<0)
          wbdie('invalid usage (current site not at k=1 !?)'); end
       end

       for k=2:L-1
          if vflag, s=sprintf('%g/%g >> trotter fit',ift,nft);
             if vflag>1,  fprintf(1,'\r %s: %g/%g ',s,k,L);
             end
          end

          if k>2 || ift==1
          Ir=load_trotter_data(HAM,k); end

          if olap, AF(k)=Ir.Af; end

          [HAM,Il,Ir,Ixt(k,1,ift),Ixf(k,1,ift)] = ...
          update_trotter_2site(HAM,Il,Ir,k,'>>',oo{:});

          Il=Ir;
       end

       Ir=load_trotter_data(HAM,L);

       if olap, AF(L)=Ir.Af; end
       if ift<nft, k1=3; else k1=2; end

       for k=L:-1:k1
          if vflag, s=sprintf('%g/%g << trotter fit',ift,nft);
             if vflag>1,  fprintf(1,'\r %s: %g/%g ',s,k,L);
             end
          end

          if k<L,
          Il=load_trotter_data(HAM,k-1); end

          [HAM,Il,Ir,Ixt(k,2,ift),Ixf(k,2,ift)] = ...
          update_trotter_2site(HAM,Il,Ir,k,'<<',oo{:});

          if ift==nft
             Q=Ir.Af;
               if fop, Q={Q,Zop,'-op:^s'}; end
               if k<L, Q={Q,XR(k+1)}; end
             XR(k)=contract(Ir.A0,'!1*',Q);
          end

          if olap && k>=3
             Q=AF(k); if k<L, Q={Q,XF}; end
             XF=contract(Ir.Af,'!1*',Q);
             if k==3, I1=load_trotter_data(HAM,1);
                Q=contract(AF(1),{AF(2),XF});
                Q=contract({Q,I1.Af,'*'},Il.Af,'*');
                xft(ift,1)=getscalar(Q);
             end
           if olap>1
             Q=Ir.At; if k<L, Q={Q,XT}; end
             XT=contract(Ir.Af,'!1*',Q);
             if k==3
                Q=contract(I1.At,{Il.At,XT});
                Q=contract({Q,I1.Af,'*'},Il.Af,'*');
                xft(ift,2)=getscalar(Q);
             end
           end
          end

          Ir=Il;
       end

       if vflag && olap, z=xft(ift,:)-1;

          if ift==nft && vflag<=1, fprintf(1,'\r%90s\r','');
          else
             a=[real(z); imag(z)];
             if olap>1
                  fprintf(1,'   %9.3g   %+10.3g %+.3gi',abs(z(1)),a(:,end));
             else fprintf(1,'   %+10.3g %+.3gi',a(:,end));
             end

             if vflag>1, nl='\n'; else nl='\r'; end
             fprintf(1,[' %10.2g %8.2g' nl], ...
                max([ Ixt(:,:,ift).svd2tr ]), ...
                max([ Ixf(:,:,ift).svd2tr ]));
          end

       elseif vflag>1, fprintf(1,'\r%90s\r',''); end

    end

    [Il,Ia]=check_normalization(Il, 'A0',      [1E-15,1E-12,1E-6]);
    [Il,Ib]=check_normalization(Il,{'At','Af'},[1E-15,1E-06,1E-2]);
    Iout.Inrm=Ib;

    if norm([Ia.wrn,Ib.wrn])==0 && vflag>1
       q=[Ia.nrm,Ib.nrm];
       wblog(' * ','normalization ok @ %.3g',max(abs(q-1)));
    end

    HAM=save_trotter_data(HAM,Il,k-1);

    Iout.runtime(itau,:)=[toc/nft nft];
    if olap
       Iout.xft(itau,:)=max(xft,[],1);
    end

  % ------------------------------------------------------------------ %
  % compute expectation value based on trotter.info.ops
  % NB! time steps need to be fully completed for this;
  % now, the first operator was already applied when calling
  % setupTrotter2()! therefore can evaluate simple
  % expectation values with A0 here!
    if nops
       for k=1:L, if k>1
          Il=load_trotter_data(HAM,k); end

          Af=Il.Af; if k>1, Af=contract(XL,Af); end
          Qf=Af; if k<L, Qf={Qf,XR(k+1)}; end

          for j=1:numel(ops{1})
             x=contract(Il.A0,'*',{Qf,ops{1}(j),'*','-op:^s:gop'});
             Sxt(k,j)=getscalar(x);
          end

          XL=contract(Il.A0,'!2*',Af);
       end
       Iout.Sxt(:,itau,:)=Sxt;
    end

    [HAM,ttot]=acceptTrotter2Step(HAM);

    Iout.tt(itau)=ttot;

  % ------------------------------------------------------------------ %
  % wblog_iteration

    if mod(itau,50)==1 && jid==1
       if olap
          f='%25s | %-16s | %-16s';
          s={'max(diff(Sxt)) @ rel ','   |d<f''|f>|','  d<f|t>'};
       else
          f='%-8s | %12s | %12s | %12s';
          s={'  d|Sxt|', ...
             'svd_tr2 At Af', '|d(sv)| At Af','dnrm: At Af'};
       end

       s=sprintf(['%15s  time | ' f],sprintf('dt=%g x %g',dt,ntau),s{:});
       if itau>1, fprintf(1,'\n%s\n\n',s);
       else
          l=repmat('-',1,84);
          i=find(s=='|'); i(find(s(i-1)~=' '))=[]; i(find(s(i+1)~=' '))=[];
          l(i)='+';
          fprintf(1,'%s\n%s\n%s\n',l,s,l);
       end
    end

    s={ sprintf('%g',itau), ''};
       if ttot<10, f='%7.3f';
       elseif ttot<100, f='%7.2f'; else f='%7.1f'; end
    s{2}=sprintf(f,ttot);

    s={sprintf('%s %4s %7s |',datestr(now,'HH:MM:SS'),s{:}),''};
    f={'%s',''}; if olap, f{2}=' %25s'; else f{2}=' %8.3s'; end

  % Iout.Ixt(:,itau)=Ixt(:,end); % too much data! // Wb,Aug15,19
    % WRN! Ixt contains useless inialized data for the case
    % when At carries no global Psi index! // IXT_INIT // Wb,Jul31,23
    % e.g., when applying a scalar operator, or when starting from
    % some inital scalar wave function that is not an eigenstate
    % cf. data Dimitris => use Ixf instead!
  % Iout.Ixt=Ixt(:,end); // Wb,Jul31,23
    Iout.Ixt=Ixf(:,end);

    if nops && itau>1
       q=Iout.Sxt(:,itau-1:itau,:); x=reshape(diff(q,[],2),[],1); a=abs(x);
       x=x(find(a==max(a),1)); a=[abs(x), max(abs(q(:,end)))];
       if a(1), a=a(1)/a(2); else a=a(1); end
       if olap
            s{2}=sprintf('%.3g %+.3gi @ %5.3f',real(x),imag(x),a);
       else s{2}=a; f{2}(end)='g';
       end
    end

    if olap
       x=abs(xft(end,1));
       s{3}=sprintf('%.4g @ %.3g',round(x),x-round(x)); f{3}=' | %16s';
       if olap>1
          x=xft(end,2);
          s{4}=sprintf('%.3g%+3gi',real(x),imag(x)); f{4}=' | %16s';
       end
    else
       s{3}=[ max([Ixt.svd2tr]), max([Ixf.svd2tr]) ]; 
       f{3}=' |%7.1g %6.1g'; Iout.svd2tr(itau,:)=s{3};

       dr=zeros(L,2,2); l=nft-1; n=[];
       n1=zeros(L,2);

       for i=1:L
       for j=1:2
          for r=1:2
             if r==1, q={ Ixt(i,j,l).svd, Ixt(i,j,l+1).svd };
             else     q={ Ixf(i,j,l).svd, Ixf(i,j,l+1).svd }; end

             for p=1:2, 
                if size(q{p},2)==2
                   q{p}=q{p}(:,1).*sqrt(q{p}(:,2));
                end
                n(p)=numel(q{p}); if n(p)
                n1(i,2*(p-1)+j,r)=norm(q{p}); end
             end
             if all(n)
                if diff(n)
                   if diff(n)>0
                        q{1}(end+1:n(2))=0;
                   else q{2}(end+1:n(1))=0; end
                end
                dr(i,j,r)=norm(q{1}-q{2});
             end
          end
       end
       end

       s{4}=max(reshape(dr,[],2),[],1);
       f{4}=' |%7.1g %6.1g'; Iout.dSVD(itau,:)=s{4};

       q=cell(1,2); for p=1:2
         q{p}=n1(:,3:4,p); i=find(q{p}); q{p}=max(abs(q{p}(i)-1)); end

       s{5}=[q{:}];
       f{5}=' |%7.1g %6.1g'; Iout.dnrm1(itau,:)=s{5};
    end

    fprintf(1,[f{:} '\n'],s{:});

  end

end

% -------------------------------------------------------------------- %
% updating neighboring pair of sites [k-1,k]

function [HAM,Il,Ir,Ixt,Ixf]=update_trotter_2site(HAM,Il,Ir,k,sdir,varargin)

  if     isequal(sdir,'>>'), qdir=+1;
  elseif isequal(sdir,'<<'), qdir=-1; else sdir, wbdie('invalid sdir'); end

  L=length(HAM); e1=HAM.oez(1).op;

  oo={'itag',Il.Af.info.itags{2},varargin{:}};

  ul=HAM.user.trotter.mpo(k-1); if k>2, Ql={Il.Xm,ul}; else Ql=ul; end
  Al=getIdentity(Il.Af,1,'-A:l',e1);
  Xl=contract(Al,'*',{Ql,Il.At},[1 3 2]); % LL',m(po)

  ur=HAM.user.trotter.mpo(k); if k<L, Qr={ur,Ir.Xm}; else Qr=ur; end
  Ar=getIdentity(Ir.Af,2,'-A:r',e1);
  Xr=contract(Ar,'*',{Ir.At,Qr}); % RR',m(po)

  if numel(Xl.Q)>4 || numel(Xr.Q)>4
     wbdie('unexpected rank (%g/%g)', numel(Xl.Q), numel(Xr.Q));
  end

  X2=contract(Xl,Xr);

  if numel(Il.At.Q)>3 || numel(Ir.At.Q)>3
     El=getIdentity(Il.At,1,'-A:l',e1); Yl=contract(El,'!3*',Il.At);
     Er=getIdentity(Ir.At,2,'-A:r',e1); Yr=contract(Er,'!3*',Ir.At);
     Y2=contract(Yl,Yr);

     if qdir>0
        [x,u,Ixt]=orthoQS(Y2,1,'>>',oo{:});
        Il.At=contract(El,u,[1 3 2]);
        Ir.At=contract(Er,x,[3 1 2]);
     else
        i=finditag(Y2,'r$',1); if i==3, Y2=permute(Y2,'132'); end

        [x,u,Ixt]=orthoQS(Y2,2,'>>',oo{:});
        Ir.At=contract(u,Er);
        Il.At=contract(El,x,[1 3 2]);
     end
  else Ixt=struct('svd2tr',0,'svd',[1]);
  end

  if qdir>0
     [x,u,Ixf]=orthoQS(X2,1,'>>',oo{:});
     Il.Af=contract(Al,u,[1 3 2]);
     Ir.Af=contract(Ar,x,[3 1 2]);

     HAM=save_trotter_data(HAM,Il,k-1);

     Ir.Xm=contract(Il.Af,'!2*',{Ql,Il.At},[1 3 2]); % LL',m(po)
  else
     i=finditag(X2,'r$',1); if i==3, X2=permute(X2,'132'); end

     [x,u,Ixf]=orthoQS(X2,2,'>>',oo{:});
     Ir.Af=contract(u,Ar);
     Il.Af=contract(Al,x,[1 3 2]);
     HAM=save_trotter_data(HAM,Ir,k);

     if k>2
     Il.Xm=contract(Ir.Af,'!1*',{Ir.At,Qr}); end % RR',m(po)
  end
end

% -------------------------------------------------------------------- %
% check/fix normalization of given wave function
% issuing warning/error for given thresholds weps
%    weps(1)   threshold for normalizing
%    weps(2)   threshold for warning
%    weps(3)   threshold for error

function [Ik,In]=check_normalization(Ik,ff,weps)

  if ~iscell(ff), ff={ff}; end;
  for i=1:numel(ff), f=ff{i}; Ak=getfield(Ik,f);

     nrm=normQS(Ak); e=abs(nrm-1); q=round(nrm); q=[q, nrm-q]; wrn=0;
     if any(e>weps)
        if e>weps(2)
           s=sprintf('%s not normalized %g @ %+.3g',f,q(1),q(2));
           if e>weps(3), wbdie('%s !?',s); wrn=9;
           else wblog('WRN',s); wrn=2; end
        else wrn=1;
        end
        Ik=setfield(Ik,f,(1/nrm)*Ak);
     end
     In(i)=add2struct('-',nrm,wrn);
  end
end

% -------------------------------------------------------------------- %

