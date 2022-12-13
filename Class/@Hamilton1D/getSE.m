function [se,Ie,IC]=getSE(HAM,Il,varargin)
% function [se,Ie,IC]=getSE(HAM,Il[,opts])
%
%    Get entanglement entropy (SE) from Il data
%
% Options
%
%   '--xr2'    extrapolate data to discarded weight r2->0 (default if nsw>1)
%   '--xD'     extrapolate data to number of states D->infty
%   '--xD*'    extrapolate data to number of multiplets D*->infty
%
%   'nsw',..   only use at most first nsw sweeps (this is also useful
%              when extrapolating to check convergence)
%
%   'pe',...   polynomial degree to use for extrapolation of entropy S (2)
%   'pc',...   polynomial degree to use for extrapolation of central charge c (2)
%   'p',...    polynomial degree to use for extrapolation (both, S and c)
%
% Wb,Aug09,17

% adapted from plotSS.m

  if isempty(HAM) || nargin<2 || ~isfield(Il,'se')
     wberr('invalid Il data'); end

  getopt('init',varargin);
     m=getopt('p',[]); 
     if ~isempty(m), pe=m; pc=m;
     else
        pe=getopt('pe',2);
        pc=getopt('pc',2);
     end

     sflag= getopt('-sym');

     if     getopt('--xr2'), xflag=1;
     elseif getopt('--xD' ), xflag=2;
     elseif getopt('--xD*'), xflag=3; else xflag=0; end

     nsw_  = getopt('nsw',0);

     ksplit=getopt('dk',[]);
     aflag =getopt('-a');
     pflag =getopt('--plot');
     p21   =getopt('p21',0);
     x1se  =getopt('x1se',-0.2);

  getopt('check_error');

  s=size(Il); s(end+1:3)=1; nsw=s(3);

  if xflag, if nsw<=1
     wberr('invalid usage (IL data required for extrapolation) !?'); end
  elseif nsw>1, xflag=1; end

  if nsw_>0 && nsw_<nsw
     nsw=nsw_; s(3)=nsw;
     Il=Il(:,:,1:nsw);
  end

  [type,L,perBC]=lattice(HAM);
  gotlad=~isempty(regexp(type,'ladder'));

  se=getdatafield(Il,'se');

  if nargout>1
     Ie=add2struct('-',type,gotlad,perBC,L);
  end

if ~xflag, return; end

  for i=nsw:-1:1
     q={Il(:,:,i).Dk}; for j=1:length(q), q{j}=q{j}(:,1); end
     q=cat(2,q{:});

     Dk(i,:)=max(q,[],2);
     r2(i)=mean([Il(:,:,i).r2]);
  end

  l=s(3)+1;
  nx=min(pe+2,s(3));

  kx=find(r2<1E-10);
  if numel(kx)<nx, kx=s(3)-nx+1:s(3); end

  Dk(l,:)=inf; r2(l)=0;

  if xflag==1
       xd=r2'; xlb='\delta\rho';
  elseif xflag>1 && size(Dk,2)>1
       xd=1./Dk; xlb='1/D^{\ast}';  
  else xd=1./Dk; xlb='1/D'; end

  xk=xd(kx);

  for i=1:s(1)
  for j=1:s(2)
      [p,su,mu]=polyfit(xk,permute(se(i,j,kx),[3 1 2]),pe);
      se(i,j,l)=polyval(p,0,su,mu);
  end
  end

if nargout<2 && ~pflag, return; end

  Ie.Ix=add2struct('-',xflag,Dk,r2,xlb,xd,pe,kx);

  Ie.se=se(:,:,1:l-1);
  se_=se(:,:,l);

  rs=symrank(HAM.oez(1).op);

  if sflag, se=0.5*(se+flipdim(se,1)); end

  if ~perBC
       cfac=1/6; if isempty(ksplit), ksplit=2*rs; end
  else cfac=1/3; if isempty(ksplit), ksplit=2; end
  end

  if ksplit
     ks=1:ksplit;
     if gotlad
        if aflag, ks=fliplr(ks); else ks=ks(end); end
     end
  else
     ks=1;
  end
  nks=numel(ks);

  x1c=1.00*x1se;

if pflag
ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
setax(ah(1,1))
end

  for k=nsw+1:-1:1
     q=struct('isw',k,'Dk',[],'se',[],'p',zeros(nks,pc+1),'c',nan(nks,2));
     q.se=se(:,:,k); q.xk=xd(k);

     x=cfac*log(sin(pi*(1:length(q.se))/L));

     if pflag
        plot(x,q.se,'-'); hold on
     end

     for i1=1:nks
        i=ks(i1):max(1,ksplit):length(x); j=find(x(i)>x1se);
        if numel(unique(x(i(j))))>pc
           p=polyfit(x(i(j)),q.se(i(j)),pc);
           q.p(i1,:)=p;
           q.c(i1,:)=[ p(pc), polyval(polyder(p),x1c) ];
        end
     end
     IC(k)=q;
  end

  if pflag
     xmark(x1se,'k:');
  end

  C2=permute(cat(3,IC.c),[3 1 2]); 
  s=size(C2); c2=zeros(s(2),s(3)); dc=c2;
  for i=1:s(2)
  for j=1:s(3)
     [p,s_,m_]=polyfit(xd(kx),C2(kx,i,j),pc);
     [c2(i,j),dc(i,j)]=polyval(p,0,s_,m_);
  end

  Ie.Ic=add2struct('-',cfac,ks,x1c,x1se,pc,ksplit);
  Ie.c2=c2;
  Ie.dc=dc;

  se=se_;

% NB! keep interleaved as this is required e.g. for the extraction
% of central charge where center is across central rung
% if 0 && gotlad
%    e=norm(diff(se,[],2))/norm(se); if e>1E-3
%      wblog('WRN','got se difference across legs @ %.3g !?',e); end
%
%    se=mean(se,2);
%    if 0
%     % se(:,1) entropy across full rungs
%     % se(:,2) entropy including half rungs
%       se=fliplr(reshape(se(1:end-1),2,[])');
%       Ie.L=Ie.L/2;
%    end
% end

end

