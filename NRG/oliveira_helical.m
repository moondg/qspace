function [ff,gg,I2,I1]=oliveira_helical(ximp,Lambda,varargin)
% function [ff,gg,I2,I1]=oliveira_helical(ximp,Lambda [,opts])
%
%   Block-tridiagonalization for 2-impuritye setting within
%   single spinfull helical edge mode at impurity distance ximp
%   in units of a = pi*v/D. This routine uses units D=1 and a=1
%   which thus implies v=1/pi.
%
%   Note that for integer ximp, i.e., discrete lattice spacing
%   in units of a=1, the local f0's are exactly orthogonal.
%
% Options
%
%    'kloc',.. up to which f_k to transform symmetric to local basis
%              (default: -1, i.e., impurity only)
%
%    'Nmax',.. using energy discretization down to Lambda^(-Nmax)
%    'nz',..   zero-averaging (test for plain 1-band case only; 0) 
%    'z',..    z-shift of logarithmic discretization (0)
%    'sigma',. spectral broadening (together with non-zero nz only)
%
%    -q        quiet
%    -p        plot summary of block-tridiagonalization
%    'wfac',.. factor on width of line when plotting couplings
%              along ladder (only together with option -p; default: 5)
%    -k        stop with keyboard (wbstop)
%
% For testing purposes only
%
%   --block-trafo
%              block transformation on low-energy portion of chain
%              which shows with -p option, that chains can be transformed
%              to decoupled setting within numerical precision.
%
% adapted from oliveira.m
% Wb,Aug06,21

% actually not used
%    -v        more verbose
%
% UNIT_LATTICE_SPACING (see usage notes above)
% => lattice spacing a = pi*v/D = 1 

% example usage (e.g. Figure paper)
% [ff,gg,I2,I1]=oliveira_helical(1E6,2,'N',99,'-p'); % ,'--block-trafo'
% [ff,gg,I2,I1]=oliveira_helical(1E6,2,'N',99,'-p','btrafo',2E-7);

  vflag=1;
  getopt('init',varargin);
     if getopt({'nolog','-q'}), vflag=0;
     elseif getopt('-v'), vflag=2; end
     pflag=getopt('-p');
     kflag=getopt('-k');

     tflag=getopt('-t');

     btrafo=getopt('--block-trafo');
     if btrafo
          btrafo=1E-8;
     else btrafo=getopt('btrafo',0); end

     wfac=getopt('wfac',5);

     nz=getopt('nz',0); if nz<2, nz=0; end
     sigma=getopt('sigma',[]);
     kloc=getopt('kloc',-1);

     betaMin=getopt('betaMin',1E-16);
     Nmax=getopt('Nmax',[]);
     z=getopt('z',0);

     N=getopt('N',[]); if isempty(N)
     N=getopt('L',[]); end
  if isempty(N)
       N=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(Nmax)
     Nmax=ceil( abs( log(1E-4*betaMin) / log(Lambda) ));
  end
  if numel(Nmax)~=1 || ~isfinite(Nmax) || Nmax-round(Nmax)~=0 || Nmax<8
     Nmax, wberr('invalid Nmax');
  end

  if isempty(N), N=ceil(Nmax/2); end

  if numel(N)~=1 || ~isfinite(N) || N-round(N)~=0 || N<8
     N, wberr('invalid N');
  end
  if N>Nmax
     wblog('WRN','N too large; set N=%g -> %g (z=%g)',N,Nmax,z);
     N=Nmax;
  end

  s=sprintf('impurity distance x=%g',ximp);
  if ximp<=0, wbdie(['invalid ' s]); elseif vflag, wblog('<i>',s); end

  if nz
     zz=linspace(0,1,nz+1); zz(end)=[];
     a0=zeros(1E5,1); emin=1E-10; emax=10; n=size(a0,1);

     om=logspace(log10(emin), log10(emax), n)';
     dx=diff(log(om(1:2)));

     for i=1:nz
        [ek,vk]=get_star_1(Lambda,zz(i),Nmax,vflag);

        x=ek;
           x(find(x<emin))=emin; x=round((log(x)-log(emin))/dx);
           x(find(x<1))=1;
           x(find(x>n))=n;

        a0=a0 + sparse(x,1, abs(vk).^2,n,1);
     end
     a0=a0/nz;
  end

  [d0,d1]=get_star_1(Lambda,z,Nmax,vflag);

  H0=sparse(diag([d0;-d0]));    D=size(H0,1);
  u=[d1;d1]; u=u/norm(u);       UU=zeros(D); UU(:,1)=u; 

  nconv=0; iconv=0;
  alpha=zeros(D,1); beta=zeros(D,1);

  for i=1:D
     u=UU(:,i); v=H0*u; alpha(i)=u'*v;
     U=UU(:,1:i);
     v=v-U*(U'*v); v=v-U*(U'*v);

     beta(i)=norm(v); if i>=D, break; end
     if beta(i)<=betaMin, nconv=nconv+1; 
        if nconv==1, iconv=i; if vflag, wblog('<i>',...
          'Lanczos converged at iteration %d/%d @ %.3g',i,D,beta(i)); end
        end
        if ~beta(i), wblog(' * ','got beta(%g)=%g - stop',i,beta(i)); break
        elseif ~pflag, break, end
     elseif ~pflag && i>=N, break;
     end
     UU(:,i+1)=v/beta(i);
  end

  ilast=i;

  j=1:i-1; Htri=sparse(j,j+1, beta(1:i-1), i,i);
  j=1:i;   Htri=Htri+Htri'+sparse(j,j,alpha(1:i),i,i);

  q=U'*U; e=norm(q-eye(size(q)));
  if e>1E-8, wblog('WRN','U not untary (e=%.3g)',e); end

  if iconv, l=iconv; else l=ilast; end
  gg=alpha(1:l);
  ff=[ 1; beta(1:l-1)];

  I1=add2struct('-',N,Nmax,betaMin,d0,d1,U,alpha,beta,ilast,iconv,gg,ff);

if pflag
  Hx=U*Htri*U';

ah=smaxis(4,2,'tag',mfilename); addt2fig Wb
header('%M :: x=%g, \Lambda=%g, z=%g, N = %g / %g',ximp,Lambda,z,N,Nmax);

  mvaxis(ah(3:4,:),[0 0.01]);

  for i=1:4, mvaxis(ah(i,:),[0 0.02*i]); end

   p=[ get(ah(3,1),'Pos'); get(ah(3,2),'Pos') ];
   ah(5,1)=axes('Position', [ p(1), 0.04, p(2,1)+p(2,3)-p(1), 0.03 ]);

   yl=min(abs([ max(alpha(end-2:end,1)), max(beta(end-2:end,1)) ]));
   yl(2)=max(beta(:,1)); q=exp(0.1*diff(log(yl)));
   yl=yl.*[1/q q];

 setax(ah(1,1))
   i=find(alpha>0); semilogy(i, alpha(i),'+-','Disp','\alpha>0'); hold on
   i=find(alpha<0); semilogy(i,-alpha(i),'o-','Disp','\alpha<0');
   h=semilogy(beta,'o-'); set(h(1),'Disp','\beta');

   xtight(1.05); ylim(yl); sms(3);
   ymark(betaMin,'-','Color',[.7 .7 .7]);
   ymark(1E-16,'k:');

   if iconv, xmark(iconv,'Color',[.9 .9 .9]); end
   legdisp

   label('iteration step n','tridiagonal entries','one-impurity setting');

 setax(ah(1,2));
     mp(Hx,'-cm');

 setax(ah(2,2)); 
     mp(Hx-H0,'-cm');

 if nz
 setax(ah(3,1));
    if isempty(sigma), sigma=log(Lambda)/(nz/8); end
    eps=100*emin;
    [ox,ax,Is]=getSmoothSpec([om;-flip(om)], [a0; flip(a0)], ...
      'eps',eps,'nlog',512,'sigma',sigma);

    plot(ox,ax); xlim([-1.2 1.2]);
    label('\omega','A(\omega)',sprintf('n_z=%g',nz));

 setax(ah(3,2)); 

    plot(ox,ax); togglelogx; togglelogy
    xlim([emin/10, 10*emax]);
    ylim([1E-16, 10]);
    xmark([emin, emax]);
    xmark(eps,'r:');
 end
else ah=[]; end

 % --------------------------------------------- %
   [d0,d1]=get_star_helical(ximp,Lambda,z,Nmax,vflag);

   H0=diag([d0;d0;-d0;-d0]); 

   UH=(1/sqrt(2))*[1 1; 1 -1];

if 1
   H0=[ zeros(2) d1'; d1 H0 ];  D2=size(H0,1); D=D2/2;
   u=zeros(D2,2); u(1,1)=1; u(2,2)=1;

   EScale=Lambda.^(-(0:.5:Nmax)');
else
   u=d1;

   q=u'*u; u_=u;
   wblog(' * ','initial f0 overlap @ %.3g x (x=%.4g)',x*abs(q(2)),x);

   u=u*diag(1./sqrt(diag(u'*u)));
   e=norm(u'*u-eye(2),'fro'); if e>1E-14
     s=sprintf('got non-orthogonal f0 modes @ %.3g',e);
     if e>1E-8, wbdie(s); else wblog('WRN',s); end
   else wblog(' * ','got orthogonal f0 modes @ %.2g',e);
   end

   EScale=Lambda.^(-(0.5:.5:Nmax)');
end

   UU=zeros(D2,2,D); UU(:,:,1)=u;

   nconv=0; iconv=0;
   alpha=cell(D,1); beta=cell(D,1);

   for i=1:D
      u=UU(:,:,i); v=H0*u; alpha{i}=u'*v;
      U=reshape(UU(:,:,1:i),D2,[]);
      v=v-U*(U'*v); v=v-U*(U'*v); v_=v;

      v=v*diag(1./sqrt(diag(v'*v)));
      v=v*UH; v=v*diag(1./sqrt(diag(v'*v))); v=v*UH;

      beta{i}=v'*v_; if i>=D, break; end
      bi=norm(beta{i},'fro');

      if bi<=betaMin, nconv=nconv+1; 
         if nconv==1, iconv=i; if vflag, wblog(' * ',...
           'Lanczos converged at iteration ib=%d/%d @ %.3g',i,D,bi); end
         end
         if ~bi, wblog(' * ','got beta(%g)=%g - stop',i,bi); break
         elseif ~pflag, break, end
      end
      UU(:,:,i+1)=v;
   end

   r=sin(pi*ximp)/(pi*ximp); q=diag(sqrt([1+r, 1-r]));
   e=norm(beta{1}-q);
   if e>1E-12, wblog('WRN','got unexpected beta{1} @ e=%.3g',e); end

   ilast=i;
   EScale=EScale(1:ilast);

   if iconv, l=2*iconv; else l=ilast; end
   q=U(:,1:l)'*U(:,1:l); e=norm(q-eye(size(q)));
   if e>1E-8, wblog('WRN','U not untary (e=%.3g @ k=%g)',e,l/2); end

   HB=blkdiag(alpha{:});
   for k=1:ilast-1, i=2*(k-1)+(1:2);
      HB(i,i+2)=beta{k};
      HB(i+2,i)=beta{k}';
   end

   e=norm(H0 - U*HB*U');
   if e<1E-14, if vflag
        wblog(' * ',['block Lanczos consistent @ %.3g (D=%d)\n' ... 
       'truncating Wilson chain length to L=%d / %d ~ D/4'],e,D2,N,Nmax); end
   elseif e<1E-8
        wblog('WRN','block-diagonalization consistent @ %.3g',e);
   else wblog('ERR','block-diagonalization failed (e=%.3g)',e);
   end

   if 0
      k1=round(e2iter(1E-6,Lambda));
      for k=k1:ilast
         if mod(k,2), continue; end

         alpha{k}=UH*alpha{k}*UH;
         beta{k} =UH*beta{k} *UH;
      end
   end

   if btrafo
      k=round(e2iter(btrafo,Lambda)); i=2*k-1; ktrafo=k;
      [a,b]=block_trafo(HB(i:end,i:end));
      alpha(k:k+numel(a)-1)=a;
      beta( k:k+numel(b)-1)=b;
   else ktrafo=0;
   end

if pflag
 setax(ah(2,1));

   q=U'*U; 
   xd=0.5:0.5:D; i=1:(2*ilast);
   plot(xd(i), diag(q),'Disp','normalization'); hold on

   q=q-eye(size(q)); e=norm(q);
   if e>1E-8, wblog('WRN','U not untary (e=%.3g)',e); end

   plot(xd(i),(1/D)*sum(abs(q).^2,1),'Disp','non-orthogonality');
   legdisp({'W',[0 0]}); xlim([0 D]);
   ytight(1.1);

   yl=zeros(2,2,2);
   s={'11','21','12','22'};
   L=min(round(e2iter(1E-12,Lambda)),D);

   for p=1:2
      if p==1, 
           x3=reshape(cat(3,alpha{:}),4,[]).'; s1='\alpha';
      else x3=reshape(cat(3,beta{:} ),4,[]).'; s1='\beta';
      end

      setax(ah(3,p));
         h=plot(abs(x3)); mv2back(h(3:4));
         xlim([0,ilast+2]); togglelogy; yl(p,:,1)=ylim;
      setax(ah(4,p));
         h(:,2)=plot(x3./repmat(EScale,1,4)); mv2back(h(3:4,2));
         xlim([0,L+2]); yl(p,:,2)=ylim;

      for i=1:numel(s)
         set(h(i,:),'Disp',sprintf('%s_{%s}',s1,s{i}));
      end

      setax(ah(3,p)); 
      if 1 || p==1
         set(h(4,:),'Color',get(h(1),'Color'));
         set(h(3,:),'Color',get(h(2),'Color'));
             blurl(h(4),'LineW',4); blurl(h(4,2)); % ,'LineW',2);
             blurl(h(3),'LineW',4); blurl(h(3,2)); % ,'LineW',2);
         o={'pick',[3 4 1 2]};
      else o={}; end

      legdisp({'NE',[0.01 0.01]},o{:},'ysc',0.7,'-detach');

      if p==2, hold on
         h=plot(abs(diff(x3(:,[1 4]),[],2)),'Color',[.8 .8 .8],...
         'Disp',sprintf('|%s_{%s} - %s_{%s}|',s1,s{1},s1,s{4}));
         legdisp({'SW',[0 0]},'pick',5,'-detach');
      end
      ymark(1/ximp,'k:','istr',{'1/x',0.5});
      xmark(e2iter(1/ximp,Lambda),'k:')
   end

   for p=1:2
      q=[min(yl(:,1,p)), max(yl(:,2,p))];
      if p==1 
         q(1)=max(1E-40, q(1));
         set(ah(3,:),'YLim',q,'YTick',logspace(-40,0,5));
      else
         set(ah(4,:),'YLim',q);
      end
   end

   setax(ah(3,1)); label('iteration step n','diagonal block \alpha');
   title2('decoupled antisymmetric channel for n>%g ?',Nmax);
   setax(ah(3,2)); label('iteration step n','off-diagonal block \beta');

   for p=1:2, setax(ah(4,p));
      if p==1, s='\alpha';
      else     s='\beta'; end
      label('iteration step n',[s ' / \omega_n']);
      if ktrafo
         for i=1:2
            xmark(ktrafo,'istr',{'block-trafo',0.6,'Rotat',90});
            if i==1, setax(ah(3,p)); end
         end
      end
   end

 setax(ah(5,1));

   plot(1:L,repmat(1,1,L),'ko'); hold on
   plot(1:L,repmat(2,1,L),'ko'); sms(3);
   xlim([0 L+1]); axis off

   for k=[1 10:10:L]
      text(k,2.5,num2str(k),'HorizontalAl','center','FontSize',8);
   end

   q=reshape(cat(3,beta{1:L}),4,[])';
   q=max(abs(q(:,[2 3]) ./ repmat(EScale(1:L),1,2)));
   postext({'N',[0 2]},'rescaled Creutz-couplings <= %.3g',max(q));

   for k=1:L
      a=alpha{k}/EScale(k);
      e=norm(a-a'); if e>1E-4, wblog('WRN',...
       'alpha(%g) non-hermitian @ %.3g',k,e); end
      for i=1:2
         if abs(a(i,i))>1E-3
            if a(i,i)>0, o='bo'; else o='ro'; end
            plot(k,i,o,'MarkerSize',wfac*abs(a(i,i)));
         end
      end
      if abs(a(1,2))>1E-3
         if a(1,2)>0, o='k-'; else o='r-'; end
         plot([k k],[1 2],o,'LineW',wfac*abs(a(1,2)));
      end

      b=beta{k}/EScale(k);
      for i=1:2
         if abs(b(i,i))>1E-3
            if b(i,i)>0, o='k-'; else o='r-'; end
            plot([k-1, k],[i i],o,'LineW',wfac*abs(b(i,i)));
         end
         j=3-i;
         if abs(b(i,j))>1E-3
            if b(i,j)>0, o='k-'; else o='r-'; end
            plot([k-1, k],[i j],o,'LineW',wfac*abs(b(i,j)));
         end
      end
   end

   if ktrafo, xmark(ktrafo,'istr',{'block-trafo',0.6,'Rotat',90}); end
end

   I2=add2struct('-',N,Nmax,betaMin,...
      ximp,'v=1/pi',... % 'units=''D=1, v=1/pi'''
      d0,d1,U,HB,ilast,iconv,kloc,alpha,beta);

   Uh=diag([1 -1i])*UH;
   tau_x=[0 1; 1 0];

   for k=2:kloc+2, ak=alpha{k};
      e=norm(ak-ak(2)*tau_x);
      if e>1E-14, wblog('ERR','got alpha{%g} @ e=%.3g',k,e); end
      alpha{k}=Uh'*ak*Uh;
   end
   gg=reshape(cat(3,alpha{2:N+1}),[],N).';

   for k=1:kloc+2
      beta{k}=Uh'*beta{k}*Uh;
   end
   if ~isempty(k), beta{k}=beta{k}*Uh'; end
   ff=reshape(cat(3,beta{1:N}),[],N).';

   if kflag
      if btrafo
         k=ceil(e2iter(2E-12,Lambda));
         q=cat(3,alpha{k:k+2}); q=abs(diff(squeeze(q(2,1,:))));
         i=find(q==min(q)); k=k+i-1;

         n=8; i=2*[k, k+n-1]-1; i=i(1):i(2); d=numel(i);
         Hb=HB(i,i) / EScale(k);

         block_trafo(Hb,'-p');
      elseif tflag
         [om,is]=sort([d0;d0;-d0;-d0]); u2=d1(is,:);
         u2=u2*diag(1./sqrt(diag(u2'*u2)));
         a0=zeros(size(d1,1),4); l=0; S=cell(1,4);
         for j=1:2, for i=1:2, l=l+1;
            a0(:,l) = u2(:,i).*conj(u2(:,j));
            S{l}=sprintf('%g%g',i,j);
         end, end

         i=find(abs(om)>1E-16);

         if isempty(sigma), sigma=0.8*log(Lambda); end
         [ox,ax,Is]=getSmoothSpec(om(i),a0(i,:), ...
           'eps',min(d0)/10,'nlog',128,'sigma',sigma);

         gi=-pi*ax;
         gc=complex( KKreal(ox,gi), gi);
         G0=gc*kron(conj(Uh),Uh);

         t=[mfilename '-spec'];
         as=smaxis(2,2,'tag',t,'fpos',[530 220 525 485]);

        setax(as(1,1));
         plot(ox,ax); axis([-2 2, -.7 1]);
         ymark(0.5,'-','Color',[.8 .8 .8]);
         xmark([-1 1],'-','Color',[.8 .8 .8]);

        setax(as(1,2));
         h=loglog(abs(ox),abs(ax)); axis([0.1/ximp, 4, 1E-6 10]);
         for i=1:numel(h), set(h(i),'Disp',['A_{' S{i} '}']); end
         set(h(2),'LineSt','--'); hold on
         h=loglog(abs(ox),abs(mean(ax(:,[1 4]),2)),'Color',[.7 .7 .7],...
           'Disp','(A_{11}+A_{22})/2'); mv2back(h);
         legdisp({'SW',[0 0]}); 

        setax(as(2,1)); j=[1 4];
         h=plot(ox,-(1/pi)*imag(G0(:,j)));
         set(h(2),'LineSt','--'); xlim([-2 2]);
        setax(as(2,2)); j=[2 3];
         h=semilogx(abs(ox),abs(G0(:,j)));
         set(h(2),'LineSt','--'); xlim([1E-8, 10]);
      end
      wbstop
   end

end

% -------------------------------------------------------------------- %
function [alpha,beta]=block_trafo(Hb,varargin) 

   getopt('init',varargin);
      pflag=getopt('-p');
   getopt('check_error');

   if pflag, t=[mfilename '-block-trafo'];
      ah=smaxis(2,2,'tag',t,'fpos',[530 220 525 485]);
      setax(ah(1,1)); mp(Hb,'-cm');
   end

   if 0
      u=(1/sqrt(2))*[1 1; 1 -1];
      for l=1:2:n, i=2*l-1; i=i:i+1;
         Hb(i,:)=u*Hb(i,:);
         Hb(:,i)=  Hb(:,i)*u';
      end
   end

   d=size(Hb,1);

   for l=2:2:d-3, i=l:l+1;
      [u,s,v]=svd(Hb(i,i+2));
      if mod(l,4)==2
           u=blkdiag( u,-v);
      else u=blkdiag( u, v); end

      i=[i,i+2];

      Hb(i,:)=u*Hb(i,:);
      Hb(:,i)=  Hb(:,i)*u';
   end

   if pflag
      setax(ah(1,2)); mp(Hb,'-cm');
   end

   for l=4:4:d-1, i=l:l+1; j=[l+1,l];
       Hb(i,:)=Hb(j,:);
       Hb(:,i)=Hb(:,j);
   end

   if pflag
      setax(ah(2,2)); mp(Hb,'-cm');
   end

   if nargout, Hb_=Hb;
      n=floor((d-1)/2); alpha=cell(1,n); beta=cell(1,n-1);
      for i=1:n, j=2*i+[-1 0];
         alpha{i}=Hb(j,j);  Hb(j,  j)=0; if i>1
         beta{i-1}=Hb(j-2,j); Hb(j-2,j)=0; Hb(j,j-2)=0; end
      end

      e=norm(Hb);
      if e, wblog('Wrn','got residual off-diagonal blocks @ %.3g',e); end
   end
end

% -------------------------------------------------------------------- %
function [d0,d1]=get_star_1(Lambda,z,Nmax,vflag)

  z=mod(z,1);

  iz=(0:Nmax)'; if z
     if vflag, wblog(' * ','using z=%g (get_star_1)',z); end
     iz(2:end)=iz(2:end)-z;
  end

  Edisc=Lambda.^(-iz);

  d0=((1-1/Lambda)/log(Lambda)) * Edisc(1:end-1);
  d0(1)=z+(1-Lambda^(z-1))/log(Lambda);

  d1 = sqrt(abs(diff(Edisc)));
  d1 = d1/(sqrt(2)*norm(d1));

end

% -------------------------------------------------------------------- %
function [d0,d1]=get_star_helical(ximp,Lambda,z,Nmax,vflag,ah)

  z=mod(z,1);

  iz=(0:Nmax)'; if z
     if vflag, wblog(' * ','using z=%g (get_star_helical)',z); end
     iz(2:end)=iz(2:end)-z;
  end

  Edisc=Lambda.^(-iz);

  d0=((1-1/Lambda)/log(Lambda)) * Edisc(1:end-1);
  d0(1)=z+(1-Lambda^(z-1))/log(Lambda);

  dE=abs(diff(Edisc));
  Em=0.5*(Edisc(2:end)+Edisc(1:end-1));

  xDk = (ximp*pi/2)*dE;
  if ximp<=0, wbdie('invalid usage (got ximp=%g)',ximp); end
  d1=cell(1,2);

  for p=1:2
     q = (sin(xDk)./xDk) .* exp((1i*pi*ximp)*Em);
     if p==2, q=conj(q); end

     r=abs(q); 
     ph=0.5*angle(q);

     if any(r>1)
        wbdie('unexpected r in [%.3g %.3g]',min(r),max(r)); end

     sp=1+r; sm=1-r;
     if 1
        i=find(sm<1E-10); sm(i)=(1/6)*xDk(i).^2;
     end
     rp=sqrt(sp); rm=sqrt(sm);

     T=zeros(2,2,Nmax);
     for l=1:Nmax, c=cos(ph(l)); s=sin(ph(l));
        T(:,:,l)=sqrt(0.5*dE(l)) * [c -s; s c] * diag([rp(l), rm(l)]);
      % factor 0.5 due to FF1_NORM: Gamma/pi cancels, but 1/sqrt(2D) remains
      % actually 1/sqrt(2D * s0) with s0 ~=1 if ximp not an integer!
      % -> second column in T(:,:,l) diminishes for Em << 1/ximp
      % -> total normalization (sum^2 of both columns)
      %    still adds up to [ s0 2DGamma/pi ]^1/2
      %    which is 1 for integer ximp (enforced with FF1_ORTHO below!)
     end

     d1{p}=reshape(permute(T,[3 2 1]),[],2);
   % 1st column: coupling of 1st impurity to 2*Nmax bath sites  *)
   % 2nd column: coupling of 2nd impurity to 2*Nmax bath sites
   % 
   % *) 2 blocks of size Nmax each where
   %    1st block corresponds to symmetric bath channel, and
   %    2nd block corresponds to antisymmetric bath channel
  end

  d1=cat(1,d1{:});

  if abs(ximp)<1E12 && abs(mod(ximp,1))<1E-8
     e=norm(d1'*d1-eye(2));
     if e>1E-10, wblog('WRN',...
       'got unexpected non-orthogonal states @ e=%.3g',e); 
     end
  end

if nargin<6 || isempty(ah), return; end
setax(ah(2,1));

  for i=1:2
  for j=1:2, s=sprintf('T_{%g%g}',i,j);
     h=plot(squeeze(abs(T(i,j,:))),'-','Disp',s); hold on
     if i==2, set(h,'LineSt','--'); end
  end; end

  xlim([0,size(T,3)]); togglelogy; legdisp

  xmark(e2iter(1/ximp,Lambda)/2,'k:','istr',{'1/x'});

setax(ah(3,1));
  semilogy(abs(d1)); hold on
  xlim([0,size(d1,1)]);

end

% -------------------------------------------------------------------- %

