function [Hq,Iq]=tuneHamQtot(H0,Qtot,varargin)
% function [Hq,Iq]=tuneHamQtot(H0,Qtot [,opts])
%    
%    Adjust `chemical potential' (dmu) in every symmetry label
%    such that Hq = H0+dmu*NQ has low-energy expectation value
%    equal or close to Qtot.
%
%     * here Qtot can be an arbitrary float
%     * symmetry labels with nan-entries in Qtot are ignored
%     * by default, also non-abelian symmetry labels with q=0
%       are ignored (use flag -f to also include them).
%
% Wb,Aug18,18

  getopt('init',varargin);
     pflag=getopt('-p');
     bfac =getopt('bfac',1.5);
     miter=getopt('miter',8);
     fflag=getopt('-f');
  Qop=getopt('get_last',[]);

  if numel(H0.data)==1 && isequal(H0.Q{1},0)
       oR={'-q'};
  else oR={};
  end

  nq=size(H0.Q{1},2); nd=numel(H0.data);
  if ~isempty(Qop)
     if ~isQSpace(Qop), error('Wb:ERR',...
       '\n   ERR invalid usage (invalid QSpace Qop)'); end
     if numel(Qop)~=nq, error('Wb:ERR',...
       '\n   ERR invalid usage (numel(Qop)=%g/%g)',numel(Qop),nq); end
     Q=H0.Q{1};
  else
     Qop=repmat(QSpace(getIdentityQS(H0,2)),1,nq);
     Q=Qop(1).Q{1};

     for i=1:nd, s=size(Qop(1).data{i});
         for k=1:nq
             Qop(k).data{i}=diag(repmat(Q(i,k),1,s(1)));
         end
     end
  end

  [ee,Ie]=eigQS(H0);
  EK=Ie.EK; nd=numel(EK.data);

  de=zeros(1,nd);
  l=0; for i=1:nd, l=max(l,numel(EK.data{i})); end;
  l1=min(8,max(2,ceil(0.6*l)));
  for i=1:nd
      q=EK.data{i}; l=numel(q); if l>=l1
      de(i)=mean(diff(q(1:l1))); end
  end
  de=de(find(de));
  de=min(de);

  beta=bfac/de;

  EK=QSpace(EK);

  R=getrhoQS(EK,beta);
  E=getIdentityQS(R,2);

  qr=zeros(2,nq);
  for i=1:nq
     x0(i)=getscalar(QSpace(contractQS(R,'*',Qop(i))));
     qr(:,i)=[min(Q(:,i)); max(Q(:,i))];
  end

  if pflag
    m=32; x1=linspace(-4,4,m); y1=zeros(m,nq);
    for j=1:nq
       for i=1:m
          R=getrhoQS(diag(EK)+x1(i)*Qop(j),beta,oR{:});
          y1(i,j)=getscalar(QSpace(contractQS(R,'*',Qop(j))));
       end
    end
  end

  if ~isequal(size(Qtot),[1,nq])
     error('Wb:ERR','\n   ERR invalid Qtot'); end

  if ~fflag
     r=getsym(H0,'-r'); s=r; s(find(r==0))=1; ns=numel(r);
     if sum(s)~=nq
        error('Wb:ERR','\n   ERR got nq=%g/%g !?',nq,sum(s)); end
     q=cell(1,ns); for i=1:ns, q{i}=repmat(r(i),1,s(i)); end
     q=[q{:}]; Qtot(find(q & Qtot==0))=nan;
  end

  y2=zeros(miter,nq); x2=zeros(miter,nq); y2(1,:)=x0;

  for j=1:nq
     w=1; if isnan(Qtot(j)), continue; end
     for i=2:miter
        dy=flipud(y2(max(1,i-2):i-1,j)-Qtot(j));
        if dy(1)==0
           x2(i:end,j)=x2(i-1,j);
           y2(i:end,j)=y2(i-1,j); continue;
        end

        if w==1
           if     all(dy>0), x2(i,j)=x2(i-1,j)+1;
           elseif all(dy<0), x2(i,j)=x2(i-1,j)-1;
           else w=2; end
        end
        if w==2
           if dy(1)*dy(2)<0
                x2(i,j)=mean(x2([i-1,i-2],j));
           else x2(i,j)=x2(i-1,j)+0.5*diff(x2(i-2:i-1,j));
           end
        end

        R=getrhoQS(diag(EK)+x2(i,j)*Qop(j),beta,oR{:});
        y2(i,j)=getscalar(QSpace(contractQS(R,'*',Qop(j))));
     end
  end

  Hq=H0; dmu=zeros(1,nq); qtot=zeros(1,nq);
  for j=1:nq, if isnan(Qtot(j)), continue; end
     q=abs(y2(:,j)-Qtot(j)); i=find(q==min(q),1);
     dmu(j)=x2(i,j);
     Hq=Hq+dmu(j)*Qop(j);
  end

  [~,I2]=eigQS(Hq);
  R=getrhoQS(I2.EK,beta);
  for j=1:nq
     qtot(j)=getscalar(QSpace(contractQS(R,'*',Qop(j))));
  end

  Iq=add2struct('-',qr,dmu,qtot,beta);

  if ~pflag, return, end
ah=smaxis(2,1,'tag',mfilename); addt2fig Wb
header('%M :: \beta=%.4g, m_{iter}=%g',beta,miter);

setax(ah(1,1))

  xl=[min(x1), max(x1); min(x2(:)), max(x2(:))];
  xl=[min(xl(:,1)), max(xl(:,2))];

  for i=1:nq
     lo={'Color',getcolor(i)};
     h=plot(x1,y1(:,i),lo{:}); if i==1, xlim(xl); hold on; end
     set(h,'Disp',sprintf('sym-label q_{%g}',i));
     if ~isnan(Qtot(i))
        plot(x2(:,i),y2(:,i),'*',lo{:});
        plot(xl,Qtot([i i]),':',lo{:});
     end
     ymark(qr(:,i),':',lo{:});
  end
  sms(3); xmark(0,'k:'); ytight(1.1);
  legdisp

setax(ah(2,1))

  for i=1:nq, if isnan(Qtot(i)), continue; end
      lo={'Color',getcolor(i)};
      plot(x2(:,i),y2(:,i)-Qtot(i),'o-',lo{:}); hold on
  end
  sms(3); togglelogy

  l=0;
  for i=1:nq, if isnan(Qtot(i)), continue; end
     q=[y2(end,i), Qtot(i)]; l=l+1;
     postext({'NW',[0 -0.1*(l-1)]},...
        'Qtot(%g) = %g @ %.2g  having d\mu=%+.3g\n',...
        i,q(1),abs(diff(q)/max(1,abs(mean(q)))),x2(end,i));
  end
  label('{\delta}\mu','{\delta}Qtot');

end

