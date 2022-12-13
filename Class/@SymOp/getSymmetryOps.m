function [FS,Io,FX]=getSymmetryOps(F,SOP,varargin)
% function [FS,I,FX]=getSymmetryOps(F,SOP [,FF])
%
%   generate proper irreducible operator multiplet
%   consistent with symmetry operations provided in SOP.
%
%   SOP is a structure array with one entry for every symmetry
%   where the order will be preserved. The structure must contain
%   the following field:
%
%     .info   info string (e.g. 'spin SU(2)')
%     .type   symmetry type (similar to QSpace.info.type; e.g. 'SU2')
%     .Sp     cell array of raising operators (e.g. single S_+ for SU2)
%     .Sz     cell array of z-operators to use (note that in SU3
%             for example, every Sz_i is not just equal to [Sp_i,Sp_i'])
%         
%   F is a single operator within the irop to start with.
%   It will be ascencded to max weight state (if it's not there yet).
%
%   NB! the returned IROP FS contains operators all of exactly
%   the same norm equal as the input operator F, i.e.
%   |FS{i}|^2 = |F(1)|^2 (using Frobenius norm).
%   No sign conventions are adopted.
%
%   FF is an optional larger pool of orthogonal (checked!)
%   yet not necessarily normalized operators which in the end
%   will be orthonormalized wrt. to the irop generated via F.
%   The resulting (reduced) operator set will be returned as FX.
%   If FX is not specified, log message is displayed, instead.
%
% Output:
%
%   FS  completed set of operators according to symmetry operations
%       provided.
%
% See also getSymStates.cc, getSymmetryStates_100706.m
% Wb,Jun29,10

% adpated from getSymmetryStates.m

   getopt('INIT',varargin);
      vflag=getopt('-v');
      Pflag=getopt('-P');
   FF=getopt('get_last',{});

   if nargin<2
      eval(['help ' mfilename]);
      if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
   end

   nS=numel(SOP); D=[]; dz=zeros(1,nS); cr=cell(1,nS);
   for i=1:nS
      cr{i}=checkSpSz(SOP(i).Sp,SOP(i).Sz,SOP(i).type);
      dz(i)=numel(SOP(i).Sz);
   end

   if isequal(FF,'-sort')
      n=numel(F); qz=cell(n,1);
      for i=1:n, qz{i}=get_symmetries_op(F(i).op,SOP); end
      [qz,is]=sortrows(chopd(2*cat(1,qz{:})));
      is=flipud(is); qz=flipud(qz);
      FS=F(is); if nargout>1, Io=add2struct('-',qz,is); end
      return
   end

   [FS,qz]=get_multiplet(F.op,SOP,vflag);
   qo=get_symmetries_op(FS{1},SOP);

   if numel(FF)>10
      wblog('TST','got %g F-ops',numel(FF));
   end

   if ~isempty(FF), ns=numel(FS); nf=numel(FF);
      if Pflag==0
         fac=zeros(size(FF)); kk=zeros(size(FF));
         for k=1:ns
            for i=1:nf
              [is,x,e]=sameup2fac(FF(i).op,FS{k});
              if is, fac(i)=x; kk(i)=k; i=i-1; break; end
            end
            if i==nf, error('Wb:ERR',...
           '\n   ERR failed to match input operators'); end
         end
         if nargout<3
            if any(kk(:)==0) || numel(FF)~=numel(FS)
                 wblog('ERR','failed to fully match input operator set'); 
            else wblog('ok.','all input operators matched.'); end
         else FX=FF(find(kk==0)); end
      else

         for i=1:numel(FF)
            f2(i)=norm(FF(i).op(:))^2;
         end

         f2=[std(f2), mean(f2)];
         if f2(1)/f2(2)<1E-8, f2=f2(2);
         else
            wblog('WRN',...
              'got input operators of differing norm (%.3g/%g) !?',f2);
            f2=0;
         end

         fac=zeros(nf,ns); ee=zeros(1,ns);
         for k=1:ns, Q=SymOp('',FS{k});
            for i=1:nf
               [Q,fac(i,k)]=project_ops(Q,FF(i));
               if fac(i,k) && norm(Q)<1E-12, break; end
            end
            ee(k)=norm(Q);
         end

         if norm(ee)>1E-12
            ol=zeros(nf,nf);
            for i=1:nf, for j=1:nf
                ol(i,j)=olap(FF(i).op,FF(j).op);
            end, end
            e=norm(ol-diag(diag(ol)),'fro'); if e>1E-12
               error('Wb:ERR',['\n   ERR input ' ... 
              'operator set FF is not orthogonal (@%.3g)!'],e);
            end

            error('Wb:ERR', ...
              '\n   ERR failed to fully match input operator set'); 
         end

         if nargout>2, kk=ones(size(FF)); FX=FF;
            for i=1:nf
               for k=1:ns
                  [FX(i),x]=project_ops(FX(i),FS{k});
                  if x && norm(FX(i))<1E-12, kk(i)=0; break; end
               end
            end
            FX=FX(find(kk)); nf=numel(FX); kk=ones(size(FX));

            for i=1:nf
               for k=1:i-1
                  [Q,x]=project_ops(FX(i),FX(k));
                  if x && norm(Q)<1E-12, kk(i)=0; break; end
               end
            end
            FX=FX(find(kk)); nf=numel(FX); ol=zeros(nf,nf);

            for i=1:nf, for j=1:nf
                ol(i,j)=olap(FX(i).op,FX(j).op);
            end, end

            if norm(ol-diag(diag(ol)),'fro')>1E-12
               n=size(ol,2); mark=zeros(1,n);
               for i=1:n, q=ol(:,i);
                  if norm(q)<1E-12
                     mark(i)=-1;
                  else q(i)=0;
                     if norm(q)<1E-12
                        mark(i)=+1;
                     end
                  end
               end
               F2=FX(find(mark>0)); i=find(mark==0);
               FX=FX(i);
               ol=ol(i,i);

             % WRN! eig() this introduces arbitrary signs in o,
             % and hence in FX => rather use qr() instead!
             % as this is equivalent to Gram-Schmitt, while the
             % sign in diag(r) needs to be kept positive! (see below)
             % [u,o]=eig(ol); o=diag(o);
               [u,r_]=qr(ol); r=diag(r_);

               i=find(abs(r)>1E-12); r=r(i); u=u(:,i);

               i=find(r<0);
               if ~isempty(i), u(:,i)=-u(:,i); r(i)=-r(i); end

               if ~isempty(ol)
                  for k=1:size(u,2)
                   % i=find(abs(u(:,k))>1E-3);
                   % [FQ,qz]=get_multiplet(FX(i(1)).op,SOP,vflag);
                   % if same_op_space(FX(i),FQ)
                   %    for j=1:numel(i)
                   %    end
                   % end

                     Q=FX(1).op*u(1,k);
                     for i=2:size(u,1), Q = Q + u(i,k) * FX(i).op; end
                     F2(end+1)=SymOp(sprintf('UFU''(%g@%.0f)',k,r(k)),Q);
                  end
               end
               FX=F2;
            end

            if f2
               for i=1:numel(FX)
                   FX(i).op=sqrt(f2/norm(FX(i).op(:))^2) * FX(i).op;
               end
            end
         end

      end
   end

   qz=mat2cell(qz,size(qz,1),dz); qq={}; QZ=[];
   Io=add2struct('-',qo,qq,qz,dz,QZ,'fac?','kk?','ol?');

   iz=1; eps_=1E-12;
   for i=1:nS, I=SOP(i);

      if ~isempty(I.qfac), n=numel(I.qfac);
         if n==1, qz{i}=qz{i}*I.qfac;
         elseif n>1, qz{i}=qz{i}*diag(I.qfac);
         end
      end

      t=I.type; if ~isempty(regexp(t,'^Z\d+$'))
         t='ZN'; N=str2num(I.type(2:end));
      end

      switch t

        case {'A','P','ZN'}
           if norm(diff(qz{i}))>eps_
              error('Wb:ERR','invalid Abelian z-labels'); end
           qq{end+1}=max(qz{i},[],1);

        case 'SU2'
           if norm(diff(uniquerows(sort(chopd(qz{i}))),2))>eps_
              error('Wb:ERR','invalid SU2 z-labels'); end
           qq{end+1}=max(qz{i},[],1);

        case {'SU3','SU4','SU5','SU6','SU7','SU8','SU9'}

           q3=qz{i}; r=str2num(I.type(3:end))-1;

           q=sortrows(fliplr(chopd(q3)));
           q=fliplr(q(end,:));

           if numel(q)~=r, error('Wb:ERR','invalid %s z-labels',I.type); end

           qq{end+1}=q;

        case {'Sp4','Sp6','Sp8','Sp10'}
           q3=qz{i}; r=str2num(I.type(3:end))/2;

           q=sortrows(fliplr(chopd(q3)));
           q=fliplr(q(end,:));

           if numel(q)~=r, error('Wb:ERR','invalid %s z-labels',I.type); end

           qq{end+1}=q;

        otherwise, error('Wb:ERR','invalid symmetry (%s)',I.type);
      end
      iz=iz+numel(I.Sz);

      if ~isempty(I.jmap), qq{end}=qq{end}*I.jmap; end

   end

   Io.qq=qq;

   for i=1:length(qq)
      qq{i}=repmat(qq{i},size(qz{i},1),1);
      if is_abelian_symmetry(SOP(i).type)
         qz{i}=zeros(size(qz{i},1),0);
      end
   end
   qq(2,:)=qz; Io.QZ=chopd(cat(2,qq{:}));

end

% ------------------------------------------------------------------- %
function G=get_disc_unitary(t,Z)

   if isequal(t,'P'), G=Z;
   elseif regexp(t,'^Z(\d+)$'), n=str2num(t(2:end)); 
      if ~isreal(Z) || norm(Z-diag(diag(Z)),'fro')
         error('Wb:ERR','\n   ERR invalid z-op !?'); end
      G=diag(exp((2i*pi/n)*diag(Z)));
   else
      error('Wb:ERR','\n   ERR got sym=%s !?',t);
   end
end

function z=get_disc_label(t,g)

   if isequal(t,'P'), z=g;
   elseif regexp(t,'^Z(\d+)$'), n=str2num(t(2:end));
      e=norm(abs(g)-1); if e>1E-12
         error('Wb:ERR','\n   ERR invalid g-values (e=%.3g) !?',e); end
      z=imag(log(g))*(n/(2*pi));
   else error('Wb:ERR','\n   ERR got sym=%s !?',t);
   end

end

% ------------------------------------------------------------------- %

function qq=get_symmetries_op(F,SOP)

   if nargin~=2, error('Wb:ERR','\n   ERR invalid usage'); end
   nS=numel(SOP); qq=[];

   for i=1:nS, I=SOP(i); Sz=I.Sz; m=numel(Sz);
      for j=1:m, isdisc=isequal(Sz(j).type,'disc');
         if isdisc
              G=get_disc_unitary(SOP(i).type,Sz(j).op);
              q=G*F*G';
         else q=comm(Sz(j).op,F); end

         [is,fac,e]=sameup2fac(q,F); if ~is, error('Wb:ERR',...
            'input op without well-defined symmetry labels'); end
         if isdisc, fac=get_disc_label(SOP(i).type,fac); end

         qq(end+1)=fac;
      end

      if 0
         if ~isempty(I.jmap)
            qq(end-m+1:end)=q(end-m+1:end)*I.jmap;
         elseif ~isempty(I.qfac)
            if ~isempty(I.qfac)
                x=I.qfac; if numel(x)>1, x=diag(x); end
                qq(end-m+1:end)=q(end-m+1:end)*x;
            end
         end
      end
   end

end

% ------------------------------------------------------------------- %

function [F,qz]=get_multiplet(F,SOP,vflag)

   if nargin<2 || nargin>3
      error('Wb:ERR','\n   ERR invalid usage'); end
   if nargin<3, vflag=0; end

   if ~isnumeric(F), error('Wb:ERR','\n   ERR invalid usage'); end

   Fnrm=norm(F,'fro'); Fin=F; F=F/Fnrm; 

   Sp={SOP.Sp}; for i=1:numel(Sp), Sp{i}=Sp{i}(:); end
   Sp=cat(1,Sp{:}); np=numel(Sp); found=1; nup=0;

   while found, found=0;
     for i=1:np
        X=comm(Sp(i).op,F); x=norm(X,'fro');
        if x>1E-12, F=X/x; found=1; nup=nup+1; end
     end
   end
   if nup && vflag, wblog('NB!',...
     'applied %g Sp ops to get MW seed (op)',nup);
   end

   F={F}; found=1;

   while found, found=0; nF=numel(F);
     for i=1:np, Sm=Sp(i).op';
        for k=1:nF
           X=comm(Sm,F{k}); if norm(X,'fro')<1E-12, continue; end

           for k2=1:numel(F)
              X=X-F{k2}*olap(F{k2},X);
           end
           x=norm(X,'fro'); if x<1E-12, continue; end
           F{end+1}=X/x; found=found+1;
        end
     end
   end

   for i=1:numel(F), F{i}=Fnrm*F{i}; end

   nF=numel(F); qz=cell(1,nF);
   for k=1:nF
      qz{k}=get_symmetries_op(F{k},SOP);
   end
   qz=cat(1,qz{:});

end

% ------------------------------------------------------------------- %

function F=get_MW_op(F,SOP)

   nS=numel(SOP); found=1; nup=0;

   while found, found=0;
     for i=1:nS, Sp=SOP(i).Sp; m=numel(Sp);
        for j=1:m
           q=comm(Sp{j},F); if norm(q,'fro')<1E-12, continue; end
           F=q; found=1; nup=nup+1;
        end
     end
   end

   if nup, wblog(' * ','applied %g Sp ops',nup); end

end

% ------------------------------------------------------------------- %

