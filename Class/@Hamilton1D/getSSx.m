function SS=getSSx(HAM,Dk,varargin)
% function SS=getSSx(HAM, DK_or_IL [,opts])
%
% Adapted from / based on getSS.m
% Wb,Jul31,17

  pflag=0; ah=[];

  getopt('init',varargin);
     m  =getopt('m',2);
     dfull=getopt('--Dfull');
     kflag=getopt('-k');

     if getopt('-p'), pflag=1; else
        ah=getopt('ah',[]);
        if ~isempty(ah), pflag=2; end
     end
  getopt('check_error');

  if isempty(HAM) || ~isfield(HAM.info,'HSS')
     wbdie('invalid usage (missing field HAM.info.HSS)'); end

  [type,L,perBC]=lattice(HAM);

  if ~isequal(type,'ladder');
     wbdie('non-ladder not yet implemented');
  end

  if isfield(Dk,'Dk') && size(Dk,3)>1
     Dk=reshape({Dk.Dk},size(Dk));
     for k=1:size(Dk,3)
         q=cat(3,Dk{:,:,k});
         Dk{1,1,k}=max(permute(q(:,1,:),[3 1 2]),[],1);
     end
     Dk=cat(1,Dk{1,1,:});
  else
     s=size(Dk);
     if s(1)<s(2), Dk=Dk'; end
  end

  nsw=numel(HAM.info.HSS);
  if size(Dk,1)~=nsw, whos Dk
     wbdie('length mismatch (Dk: %g/%g) !?',size(Dk,1),nsw);
  end

  for isw=1:nsw
     [S1{isw},S2{isw}]=getSS(HAM,isw);
  end

  if dfull || size(Dk,2)<2
       Dk=Dk(:,end); xlb='1/D';
  else Dk=Dk(:,1); xlb='1/D^\ast'; end

  SS={ cat(2,S1{:}), permute(mean(cat(3,S2{:}),2),[1 3 2]) };
  Sl={'rung','leg'};

  nS=numel(SS);

  jf=max(1,min(floor(nsw/2),nsw-4)):nsw; jend=nsw+1;
  x=1./Dk(jf,1)'; x2=max(x);

  for k=1:nS, Sk=SS{k};
     if size(Sk,2)~=nsw
        wbdie('length mismatch (%g: %g/%g) !?',k,size(Sk,2),nsw);
     end
     for i=1:size(Sk,1)
        [p,s,mu]=polyfit(x,Sk(i,jf),m);
        Sk(i,jend)=polyval(p,0,s,mu);
     end
     SS{k}=Sk;
  end

  if ~pflag, return; end

if isempty(ah)
     ah=smaxis(nS,2,'tag',mfilename); header('%M :: m_p=%g',m); addt2fig Wb
else ah(end+1:2*nS)=ah(end);
end

setax(ah(1,1))

  x=[1./Dk(:,1); 0];
  l=size(SS{1},1); a=max(5/l,0.1); ic=round((0.5+[-a a])*l); ic=ic(1):ic(2);

  po={'FaceColor',getcolor(3),'EdgeColor','none'};

  for k=1:nS
     setax(ah(k,1));
       h=plot(x,SS{k}','o-'); sms(2); hold on
       label(xlb,['<S{\cdot}S>_{' Sl{k} '}']);
       xl=xlim; if xl(2)>2*x2, xl(2)=2*x2; xlim(xl); end
       ytight(1.2,'view'); xmark(x2,'k:');

       y=SS{k}(ic,end); y=[min(y),max(y)];
       patch(xl([1 2 2 1 1]),y([1 1 2 2 1]),repmat(-0.1,1,5),'r',po{:});

     setax(ah(k,2));
       h=plot(SS{k},'.-'); hold on
       i=1:min(jf)-1;
       blurl(h(i),'cfac',0.5,'LineW',1); set(h(i),'Marker','none');
       set(h(end),'Color','k','LineW',1.5);
       xl=xlim;
       patch(xl([1 2 2 1 1]),y([1 1 2 2 1]),repmat(-0.1,1,5),'r',po{:});
  end

  if kflag, keyboard; end

end

