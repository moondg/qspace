function plot(HAM,varargin)
% function plot(HAM [,opts])
% Wb,Aug27,15

  getopt('init',varargin);
   % slightly randomize positions to double check
   % whether lines are lying on top of each other
     if getopt('-r'), rfac=0.05;
     else rfac=getopt('rfac',0); end
     mpflag=getopt('--map');
     vflag=getopt('-v');

     dx2=getopt('dx',0.025); nxb=0;

  getopt('check_error');

  N=numel(HAM.mpo);
  full_mpo=(isa(HAM.mpo,'QSpace') || isfield(HAM.mpo,'Q'));

  if numel(HAM)~=1, wbdie('invalid usage'); end
  i3=[];

  if isfield(HAM.info,'HH'), HH=HAM.info.HH;
     if full_mpo && isfield(HAM.info,'HS'), HS=HAM.info.HS; end
  elseif isfield(HAM.info,'HS')
     HS=HAM.info.HS; [HH,m]=convert_HS(HS,3);
     i3=find(m==3);
  else wbdie('invalid usage (missing field HAM.info.HH)'); end

  if ~N
     p=HAM.info.param;
     if isfield(p,'N') N=p.N;
     elseif isfield(p,'L') N=p.L;
        if isfield(p,'W'), N=N*p.W; end
     end
     wblog('WRN','HAM.mpo not yet initialized (using N=%g)',N);
  end

if mpflag
  XY=HAM.info.XY;

  [H24,I24]=uniquerows(sort(HH(:,[2 4]),2));

  m=1+numel(I24);
  if m<=2
       ah=smaxis(1,2,'tag',mfilename,'landscape');
  else ah=smaxis(2,ceil(m/2),'tag',mfilename,'landscape');
  end
else
  if full_mpo
       m=3; n=3;
  else m=2; n=2;
  end
  ah=smaxis(m,n,'tag',mfilename,'landscape');
end

  s=HAM.ops; [n,m]=size(s); s=reshape({s.info},[n m]); nm=n*m;
  if m==1
     for i=1:n, s{i}=sprintf('%g) %s',i,s{i}); end
  else
     for i=1:n
     for j=1:m, s{i,j}=sprintf('%g.%g) %s',j,i,s{i,j}); end, end
     s(end+1,:)={''};
  end

  if nm>8
       s=sprintf('\n%s',s{1:4},'...',s{end-1:end});
  else s=sprintf('\n%s',s{:});
  end

  header('%M :: %s with operator setup %s',HAM.info.istr,s); addt2fig Wb
  header({'SW',[0 0]},param2str(HAM.info.param));

setax(ah(1,1))

  if isfield(HAM.info,'XY')
     xy=HAM.info.XY;
  else xy=[]; end

  if isempty(xy)
     xy=[(1:N)',repmat(0.1,N,1)]; xy(1:2:end,2)=-0.1;
  end

  if full_mpo
       ah(1)=mergeax(ah(1:2,:));
  else ah(1)=mergeax(ah(1,:));
  end

  x=[min(xy(:,1)),max(xy(:,1))]; mflag=10;
  if diff(x)>2.6*mflag
     ah(1,1:2)=splitax(ah(1,1),[1 2],'dx',10);
     setax(ah(1,1));
  else mflag=0;
  end

  q=[];
  if isfield(HAM.info,'stype'), q=HAM.info.stype; end
  if isempty(q), q=ones(1,N); end

  [s,I,D]=uniquerows(q(:));
  for i=1:numel(D), c=getcolor(i);
     h=plot(xy(I{i},1),xy(I{i},2),'o','Color',c,'MarkerS',12); hold on
     if i>1
        set(h,'MarkerFaceC', 1-0.2*(1-c));
     end
  end

  for i=1:size(xy,1)
     text(xy(i,1),xy(i,2),0.1,sprintf('{%g}',i),...
     'Tag','sites:txt','FontSize',8,... %'Color',[.4 .4 .4],...
     'HorizontalAl','center','VerticalAl','middle');
  end

  i1=find(HH(:,1)==HH(:,3)); zc=get_zcolor(HH(i1,end));

  m={ 'o','*','^'}; %,'o'
  s=[ 12,  8,  8 ];
  z=[-.1,  0,  0 ];

  n=numel(m); n1=numel(i1);
  mm=repmat({''},n1,1); ms=repmat(6,n1,1); mz=repmat(0,n1,1); mi=repmat(0,n1,1);
  [q,I,D]=uniquerows(HH(i1,2));
  for i=1:size(q,1), l=mod(i-1,n)+1;
     mm(I{i})={m{l}};
     ms(I{i})=s(l);
     mz(I{i})=z(l);
     mi(I{i})=i;
  end

  q=zeros(1,n);

  for i=1:numel(i1), Hi=HH(i1(i),:);
     if Hi(end), k=Hi(1);
        c=1-0.5*(1-get_vcolor(zc(i))); 
        h=plot3(xy(k,1),xy(k,2),repmat(mz(i),numel(k),1),...
        mm{i},'Color',c,'MarkerFaceC',c,'MarkerS',ms(i));

        if mi(i)<=n && ~q(mi(i)), q(mi(i))=1;
           set(h,'Disp',sprintf('Hloc(%g)',Hi(2)));
        end
     end
  end

  i=1:size(xy,1);
  i=vkron(i,i);
  l=[ xy(i(:,1),1)-xy(i(:,2),1), xy(i(:,1),2)-xy(i(:,2),2) ];
  l=unique(sqrt(sum(l.^2,2)));
  if numel(l)>1 && l(2)>0.01, dx0=0.4*l(2); else dx0=0.1; end

  I2=find(HH(:,1)~=HH(:,3)); zc=get_zcolor(HH(I2,end),HH(I2,[2 4]));

  L=max(xy(:,1))-min(xy(:,1))+1;
  W=max(xy(:,2))-min(xy(:,2))+1;

  IJ=uniquerows(HH(:,[2 4 end])); if size(IJ,1)>6, IJ=[]; end
  nJ=size(IJ,1); bb=zeros(0,6);

  for i2_=1:numel(I2); i2=I2(i2_);
     k1=HH(i2,1); k2=HH(i2,3); c=get_vcolor(zc(i2_)); 
     x=xy([k1 k2],1); 
     y=xy([k1 k2],2); 

     if 1 && any(abs(diff(y))>2.1)
        j=find(y>0.6*W); y(j)=y(j)-W;
     end
     if 1 && any(abs(diff(x))>2.1)
        j=find(x>0.6*L); x(j)=x(j)-L;
     end

     l=norm([diff(x),diff(y)]);
     if l>dx0, q=0.5*(l-dx0)/l; else q=0.4; end

     xc=[mean(x); mean(y)];
     x=xc(1)+diff(x)*[-q q];
     y=xc(2)+diff(y)*[-q q];

     u=[-diff(y), diff(x)]; u=u/norm(u); if u(1)<0, u=-u; end
     q=[u, u*[x;y], 0 0 ]; e=abs(diff(q(3:4)));
     if e>1E-12, wbdie('invalid setup (e=%.3g)',e); end
     q(4:5)=sort([-u(2),u(1)]*([x;y]-q(3)*u'));

     i=find(sum((bb(:,1:3) - q(1:3)).^2,2)<1E-3);
     if ~isempty(i), b4=bb(i,4); b5=bb(i,5); e=0.01;
        j=find( q(4)>=b4 & q(4)<=b5 | ...
                q(5)>=b4 & q(5)<=b5 | q(4)-e<b4 & q(5)+e>b5);
        if ~isempty(j)
           dx=[-5:-1,1:5]';
           [~,~,Im]=matchIndex(unique(bb(i(j),6)),dx);
           dx=dx(Im.ix2);
           if ~isempty(dx)
              q(end)=dx(find(abs(dx)==min(abs(dx)),1)); nxb=nxb+1;
              dx=dx2*q(end);
              x=x+u(1)*dx;
              y=y+u(2)*dx;
           else
              if nxb<1E6, wblog('WRN',...
                'got more than 10 line segments on top of each other'); end
              nxb=nxb+1E6;
           end
        end
     end
     bb(end+1,:)=q;

     if rfac
        R=[diff(x); diff(y)]; R=R/norm(R); R(:,2)=[-R(2,1),R(1,1)];
        a=rfac*randn(1); U=R*[cos(a), -sin(a); sin(a) cos(a)]*R';
        x2=[xc,xc]+U*([x;y]-[xc,xc]);
        x=x2(1,:); y=x2(2,:);
     end

     h=plot(x,y,'Color',c);
     i=matchIndex(IJ,HH(i2,[2 4 end])); if ~isempty(i)
        set(h,'Disp',['J=' vec2str(IJ(i,end),'%.4g')]); IJ(i,:)=[];
     end
  end

  if vflag && nxb && nxb<1E6, wblog(' * ',...
    'got %g overlapping line segments (offset dx=%g)',nxb,dx2); 
  end

  xl=xtight(1.1); 
  ytight(1.2);

  if ~isempty(i3), H3=HS(i3); n3=numel(H3);
     scirc={ getsymb('\Rcirc'), getsymb('\Lcirc') };
     to={'HorizontalAl','center','VerticalAl','middle','FontSize',16,...
         'Tag','h3:txt'};

     for i=1:n3, q=H3(i);
        H3(i).H3aux=[q.iop(:,2)', q.hcpl];
        xi=xy(q.iop(:,1),:);
        if 1, y=xi(:,2);
           if any(abs(diff(y))>2.1), ym=max(y);
              if 0 && mean(y)/ym>.5
                   j=find(y<.5*ym); xi(j,2)=xi(j,2)+W;
              else j=find(y>.5*ym); xi(j,2)=xi(j,2)-W; end
           end
        end

        H3(i).xi=xi;
        H3(i).xc=mean(xi,1);
        H3(i).clw=(det(xi(2:end,:)-xi([1 1],:))>0);
     end

     [hx,Ix,Dx]=uniquerows(cat(1,H3.H3aux));

     ic=get(gca,'ColorOrderIndex');
     for p=1:numel(Ix), lo={'Color',getcolor(ic+p)};
        ii=Ix{p};
        for i=ii, q=H3(i);
           if 0
              xi=[q.xc; q.xi]; j=[1 2 1 3 1 4 1];
              h=plot(xi(j,1),xi(j,2),'Color',rand(1,3));
              mv2back(blurl(h,'cfac',.5,'LineW',1));
           end

           if i==ii(1)
              h=plot(q.xc(1),q.xc(2),'o',lo{:},'MarkerS',12);
              s={ sprintf('[%g,%g,%g] ',hx(p,1:3))
                   vec2str(hx(p,4),'%.4g')
              }; set(h,'Disp',[s{:}]);
           end
           if q.clw, j=1; else j=2; end
           text(q.xc(1),q.xc(2),scirc{j},lo{:},to{:});
        end
     end
  end

  if mflag
     x=[min(xy(:,1))-1.10, max(xy(:,1))+1.10]; m=mflag;
     ah(1,2)=repax(ah(1,2),ah(1,1)); xlim([x(2)-m,x(2)]);
     setax(ah(1,1));                 xlim([x(1),x(1)+m]);

     set(ah(1,2),'YTickLabel',{});
     for h=findall(ah(1,1),'type','text')'
         p=get(h,'Pos'); if p(1)>x(1)+m, set(h,'Visible','off'); end
     end
     for h=findall(ah(1,2),'type','text')'
         p=get(h,'Pos'); if p(1)<x(2)-m, set(h,'Visible','off'); end
     end
  end

  if nJ, legdisp({'NE',[0.044 0.010]},'-detach'); end % ,'-eraset'
  set(gca,'TickLength',[.010 .025]/10)

  label('X','Y','lattice structure');

if mpflag
   for k=1:numel(I24), setax(ah(k+1));
      q=sparse(HH(I24{k},1),HH(I24{k},3),HH(I24{k},end),N,N);
      mp(q,'-gca','cmap'); axis equal
      postext('NE',['HH ops = [' sprintf(' %g',H24(k,:)) ' ]']);
   end
   return
end

if full_mpo
setax(ah(3,1))
   hcpl=unique(HH(:,5));
   for i=1:numel(HS), q=HS(i);
      ic=find(hcpl==q.hcpl,1); c=getcolor(ic);
      k=q.iop(:,1);
      if numel(k)==1
           plot(k,i,'o','Color',c);
      else plot(k,repmat(i,size(k)),'Color',c);
      end
      hold on
   end
   [R,I,D]=uniquerows(abs(diff(HH(:,[1 3]),[],2)));
   n=numel(I); s=cell(1,n);
   for i=1:n, s{i}=sprintf('%g (%g)',R(i),D(i)); end
   label('k','index into info.HS',[
      sprintf('%g terms total with range {\\delta}k \\leq %g',sum(D),max(R)) ...
      10 strjoin(s,', ')]);
   sms(4);

setax(ah(3,2))
   plot(HAM.info.mpo.Nkeep,'o-'); sms(4);
   ylabel('MPO dimension')

   p=HAM.info.param;
   if isfield(p,'W')
      xmark(p.W:p.W:N+1,'-','Color',[.8 .8 .8])
      title2('W=%g / %g',p.W,W);
   end

setax(ah(3,3))
   [iop,Ib,Db]=uniquerows(HH(:,[2 4]));
   for i=1:size(iop,1)
      kk=HH(Ib{i},[1 3]);
      xx=mean(reshape(HAM.info.XY(kk(:),1),size(kk)),2);
      h=plot(xx,HH(Ib{i},5),'o-');hold on
      set(h,'Disp',sprintf('(%g,%g)',iop(i,:)));
   end
   xtight(1.1); sms(2); legdisp({'NE',[0.07 0]},'-detach','xsc',0.7); 
   label('(avg.) X','coupling strength',...
         'grouped by operator indices (i_{op1},i_{op2})');
   return

elseif ~isfield(HAM.info,'mpo') || ~isfield(HAM.info.mpo,'HM')

   return
end

setax(ah(2,1))

  HM=HAM.info.mpo.HM;
  stype=[HAM.mpo.stype];

  rfac2=2*rfac;

  for i=1:numel(HM), q=HM{i}; if ~isempty(q)
     if norm(q(:,1)-q(:,3))==0
        x=q(:,1); y=i;
        if rfac
           x=x+rfac *randn(size(x));
           y=y+rfac2*randn(size(y));
        end
        h=plot(x,y,'*','LineW',2); sms(h,4);
     else
        c=getcolor(stype(q(1)));
        x=q(:,[1 3])';
        y=linspace(i,i+1,size(q,1)+1); y=repmat(y(1:end-1),2,1);
        if rfac
           x=x+rfac *randn(size(x));
           y=y+rfac2*randn(size(y));
        end
        h=plot(x,y,'o-','Color',c); sms(h,2);
     end
     hold on
  end, end

  xtight(1.05); ytight(1.1,'y1',0.1);
  y=get(gca,'XTick'); set(gca,'YTick',[1,y(y>2)]);
  title('build sequence of pseudo MPO');
  label('site k','\rightarrow  step in mpo buildup  \rightarrow');

setax(ah(2,2));

  HM=cat(1,HM{:}); xx=mean(HM(:,[1 3]),2);
  [dx,Ix,Dx]=uniquerows([round(diff(HM(:,[1 3]),[],2)), HM(:,[2 4])]);

  nl=zeros(1,N);

  for i=1:numel(Ix), o='o-';
     l=unique(dx(i,2:end));
     t=sprintf('dx=%g',dx(i));

     if dx(i)==0, s='local'; o='*';
     elseif dx(i)==1, s='NN';
     elseif dx(i)==2, s='NNN';
     else s=t; end

     h=plot(xx(Ix{i}),HM(Ix{i},end),o); hold on % ,'Color',getcolor(dx(i))
     j=dx(i)+1; nl(j)=nl(j)+1;
     if numel(l)==1
          s=HAM.ops(l,1).info;
     else s=sprintf('%s/op%s',s,sprintf('%g',l)); end
     set(h(1),'Disp',s,'Tag',t);
  end

  for i=find(nl>6)
     h=findall(gca,'Tag',sprintf('dx=%g',i-1));
     if numel(h)>6
        set(h(3),'Disp','...');
        set(h(4:end-2),'Disp','');
     end
  end

  xtight(1.05); ytight(1.1); sms(3);
  legdisp({'E',[0 0]});

  title('energy / coupling strength vs. avg. site position');
  label('(k_1+k_2)/2','coupling / interaction strength');

end

% -------------------------------------------------------------------- %
% vary color within fluctations around mean value of coupling parameters

function z=get_zcolor(z,io2)
   if isempty(z), return; end
   if nargin>1, z=z+mean(io2,2); end

   q=[min(z),max(z)]; q(2)=diff(q);
   if q(2)
        z=(z-q(1))*(2/q(2))-1;
   else z(:)=0; end
end

function c=get_vcolor(z)

   c=0.5*([z, 4*z.^3-3*z, 2*z.^2-1 ]+1);
   if any(c<0 | c>1)
      i=find(c<0); c(i)=0;
      i=find(c>1); c(i)=1;
   end
   if norm(c-1)<0.01, c=[.9 .9 .9]; end
end

% -------------------------------------------------------------------- %

