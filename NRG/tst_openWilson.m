
% Wb,Nov11,16

if ~isvar('IDMA') || ~isvar('II')
  wsys='user'; nz=4;
     
      U=0.12; Gamma=0.01; epsd=-U/2;

    % setdef('Nkeep',256);
      setdef('Nkeep',512','Etrunc',8,'N',52);
      setupSIAM_SU2x2; locRho=1; % nostore=1;

    % use fixed Nkeep to obtain smoother data for omeba<<TK (!)
    % clear Etrunc => set Etrunc=0 if needed

      zz=get_zshifts(nz); cto lma
      mat=sprintf('%s_%s_L%02.0f_z%g',wbstamp,mfilename,10*Lambda,nz);
      if isbatch
         fprintf(1,'\n>> %s\n\n',mat);
      end

      clear T II IDMA
      TT=1E-6*Lambda.^(0:0.2:1); nT=numel(TT);

      TAG={'plain','zitko','open'};

  for iz=1:nz, z=zz(iz); oz={'z',z};
     banner('iz=%g/%g',iz,nz);

     for iw=1:3
        switch iw
          case 1, [ff,~,Ic]=getNRGcoupling(Gamma,Lambda,N,oz{:},'-w');
          case 2, [ff,~,Ic]=getNRGcoupling(Gamma,Lambda,N,oz{:});
          case 3, [ff,  Ic]=openWilson(Gamma,Lambda,N,oz{:}); % ,'--plot'
        end

        rnrg

        for it=1:nT, T=TT(it);
           rdma

           [Idma.g0,Idma.a0,Idma.g1]=... % see also rnrg_gg.m
           getConductance0(om,a0,T,pi*Gamma);

           IDMA(iz,iw,it)=add2struct(Idma,Ic,ff,ox,ax,om,a0);

           if iz==nz
              om=IDMA(1,iw,it).om;
              a0=cat(3,IDMA(:,iw,it).a0); a0=sum(a0,3)/size(a0,3);
              [ox,ax,ah,Idma.smo] = rsmoothSpec(om,a0,Idma,'RAW','eps',T);

              q=add2struct('-',ox,ax,om,a0);
              [q.g0,q.a0,q.g1]=getConductance0(om,a0,Idma.T,pi*Gamma);
              II(iw,it)=q;

              mfig('-f',['./' mat '_' TAG{iw} '.pdf']);
           end
        end
     end
  end

  if isbatch, save2(['./' mat '.mat']); end

end

ah=smaxis(2,2,'tag',mfilename,'dx',0.12,'dy',0.12); addt2fig Wb
if isset('Etrunc'), s=sprintf(', Etrunc=%g',Etrunc); else s=''; end
header('%M :: U=%g, Gamma=%g, \epsilon_d=-U/%g   [\Lambda=%g, Nkeep=%g%s]',...
  U,Gamma,-U/epsd,Lambda,Nkeep,s);

  gfac=pi*Gamma; it=1; yl=[];
  h3=ah([1 3 2]); % ah(1,1), ah(1,2), ah(2,1)

  for iw=1:3; setax(h3(iw))
     for iz=1:nz
        h=semilogx(abs(IDMA(iz,iw,it).ox),gfac*IDMA(iz,iw,it).ax(:,1)); hold on
        set(h,'Disp',sprintf('z=%g',zz(iz)))
     end
     blurl(findall(gca,'type','line'));
     h=semilogx(abs(II(iw,it).ox),gfac*II(iw,it).ax(:,1),'k'); hold on
     set(h,'Disp','z-avg');
     yl(iz,:)=ylim;

     xlim([1E-7 2]); set(gca,'XTick',logspace(-10,2,7));

     xmark(IDMA(1,1,it).T,'istr','T');
     ymark(1,'k--');

     if iw==1, legdisp('Location','NE','dx',[0 -0.06]); end % '-erase'

     label('|\omega|','A(\omega)',TAG{iw})
     poslabel(iw,'NE','dx',[-0.1 -.03]);
  end

  set(h3,'YLim',[0 max(yl(:,2))]);

setax(ah(2,2));

  for iw=1:3
     s=size(IDMA); yy=cat(1,IDMA(:,iw,:).g0); yy=reshape(yy(:,1),s([1 3]));
     h=plot(TT,yy); hold on
     if iw==3, set(h,'LineSt','--'); end
   % set(h(1),'Disp','z-shifted');
  end
  blurl(findall(gca,'type','line'));

  for iw=1:3
     if iw==1, o={'k--'};
     elseif iw==2, o={'-','Color',[0 .7 0]};
     elseif iw==3, o={'b--'};
     end

     yy=cat(1,II(iw,:).g0);
     h=plot(TT,yy(:,1),o{:},'Disp',['z-avg (' TAG{iw} ')']); hold on
  end

  legdisp('Location','E');
  label('T','conductance g [e^2/h]');

  poslabel(4,'NE','dx',[-0.1 -.03]);
