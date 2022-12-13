function varargout=getSS(HAM,varargin)
% function SS=getSS(HAM [,opts])
% function [s1,s2,I]=getSS(...)
%
%    Get <S.S> correlations from HAM.info.HSS{:}
%
% Options
%
%   if IL was specified as input:
%     '--xD'  extrapolate to x = 1/D    -> 0 (default)
%     '--xD*' extrapolate to x = 1/Deff -> 0
%     '--xr2' extrapolate to x = r2     -> 0
%     'nsw',..only use at most first nsw sweeps (this is also useful
%             when extrapolating to check convergence)
%     '-p'    plot data
%
% Wb,Jun22,17

% adapted from plotHSS.m

  varargin_=varargin;

  getopt('init',varargin);
    mp    = getopt('m',2);
    nsw_  = getopt('nsw',0);

    if     getopt('--xD' ), xflag=1;
    elseif getopt('--xD*'), xflag=2;
    elseif getopt('--xr2'), xflag=3; else xflag=0; end

  varargin=getopt('get_remaining'); narg=length(varargin);

  k=[]; wsys='';

  q=zeros(1,narg);
  for i=1:narg, q(i)=ischar(varargin{i}); end
  i=find(q); n=numel(i);
  if n==1, wsys=varargin{i}; varargin(i)=[]; narg=narg-1;
  elseif n, wberr('invalid usage');
  end

  if narg==1, k=varargin{1};
  elseif narg, wberr('invalid usage (narg=%g)',narg); end

  if ~isempty(regexp(wsys,'^-'))
     wberr('invalid wsys (%s) !?',wsys); end

  if isempty(HAM) || ~isfield(HAM.info,'HSS')
     wberr('invalid usage (missing field HAM.info.HSS)'); end

  [type,L,perBC]=lattice(HAM);

  if     isequal(type, 'ladder'), gotlad=1;
  elseif isequal(type,'iladder'), gotlad=2; else gotlad=0; end

  if gotlad>1 && ~perBC
     wblog('WRN','got type ''%s'' having open boundary !?',type);
  end

  nsw=numel(HAM.info.HSS);

  if nsw_>0 && nsw_<nsw
     nsw=nsw_;
  end

  if isnumeric(k)
     if isempty(k), k=nsw;
     else
        if nsw_, { nsw_, k}, wberr('invalid usage'); end
        if k<=0, k=nsw+k; end
     end
     if xflag
        if k==nsw, s='last'; else s=sprintf('%g/%g',k,nsw); end
        wblog('WRN','IL needed to extrapolate -> k=%s',s);
     end
     xflag=0;
  elseif isfield(k,'Dk')
     IL=k; k=[];
     if ~xflag, xflag=1; end
     if size(IL,3)<nsw
        wberr('invalid usage: size(IL,3) ~= numel(HSS) !?');
     end
  else k, wberr('invalid usage'); end

  if xflag
     if xflag<3
        if xflag==2
             xstr='D^\ast\rightarrow\infty';
        else xstr='D\rightarrow\infty';
        end
     else
        xstr='\delta\rho{\rightarrow}0';
     end

     xd=zeros(nsw,1);
     for k=nsw:-1:1
        [SS{k},H]=get_SS(HAM.info,k,L);

        if xflag<3, q=cat(3,IL(:,:,k).Dk);
           if xflag==2
                xd(k)=1/max(squeeze(q(end,1,:)));
           else xd(k)=1/max(squeeze(q(1,1,:)));
           end
        else
           xd(k)=mean([IL(:,:,k).r2]);
        end
     end

     k=max(1,nsw-(mp+2)):nsw;
     xk=xd(k);

     Q=permute(cat(4,SS{:}),[4 1:3]); s=size(Q);
     Q=reshape(Q,s(1),[]); SX=zeros(1,size(Q,2));

     for i=find(sum(abs(Q),1))
        [p,su,mu]=polyfit(xk,Q(k,i),mp);
        SX(i)=polyval(p,0,su,mu);
     end
     SS_=SS;
     SS=reshape(SX,s(2:end));
  else
     [SS,H]=get_SS(HAM.info,k,L);
     if k<nsw
          xstr=sprintf('sweep %g/%g',k,nsw);
     else xstr=sprintf('using last sweep %g',k);
     end
  end

  if nargout<=1, varargout={SS}; return; end

  lstr={};
  I=add2struct('-',lstr,SS,H,wsys,k,gotlad,perBC,xflag,xstr,mp);

  if isempty(wsys)
     if size(SS,3)>1
        wblog('WRN','got size(SS,3)=%g -> missing wsys? (add)',size(SS,3));
        SS=sum(SS,3);
     end

     [s1,s2,l,I.lstr]=getSS_12(SS,gotlad,perBC);

     varargout={s1,s2};
     if nargout>2
        if gotlad
           jj=2:10; jj(find(jj==l))=[];
           if iscell(SS), Sk=SS{end}; else Sk=SS; end
           for j=jj
              q=full(diag(Sk,j));
              if norm(q(5:end-5)), gotlad=0; end
           end
        end
        if xflag, I=add2struct(I,SS_); end
        varargout{3}=I;
     end
  elseif isequal(wsys,'HBLX')
     if gotlad || perBC
        wberr('unexpected wsys=%s (%g,%g)',wsys,gotlad,perBC); end

     [sl,sr,I.sx,I.lstr]=getSS_HBLX(SS);
     varargout={sl,sr,I};
  end

end

% -------------------------------------------------------------------- %
function [SS,H]=get_SS(Ih,k,L)

  SS=Ih.HSS{k}; nS=size(SS,3);

  HH=Ih.HH;
  if nS>1
     io=HH(:,[2 4]); io=[min(io); max(io)];
     [k,I,D]=uniquerows(HH(:,2)); H=zeros(L,L,nS);
     for k=1:nS
        Hk=HH(I{k},:);
        Hk=sparse(Hk(:,1),Hk(:,3),Hk(:,end),L,L);
           Hk=Hk+triu(Hk,1)';
           H(:,:,k)=Hk;
        Sk=SS(:,:,k); i=find(Sk);
        Sk(i)=Sk(i)./Hk(i);
        SS(:,:,k)=Sk;
     end
  else
     [k,I,D]=uniquerows(HH(:,[1 2]));
     if any(D>1) wblog('WRN',...
       'got multiple operators acting on same pair of sites'); end

     Hk=sparse(HH(:,1),HH(:,3),HH(:,end),L,L); Hk=Hk+Hk';
     i=find(SS); SS(i)=SS(i)./Hk(i);
     H=Hk;
  end

  s=size(SS);
  if ~isequal(s(1:2),size(Hk))
     wberr('got size mismatch L=%g/[%g %g] !?',L,s);
  end

  if ~isreal(SS)
     e=[norm(imag(SS(:)),'fro'), norm(real(SS(:)),'fro')];
     if e(1)/e(2)>1E12, wblog('WRN','got imag(HSS) @ %.3g %.3g!?',e); end
     SS=real(SS);
  end

end

% -------------------------------------------------------------------- %
function [s1,s2,l,lstr]=getSS_12(SS,gotlad,perBC)

  if gotlad==1, L=max(size(SS));
     if ~SS(L-2,L)
        q=[sum(SS(L,L-2:L-1)), SS(L-1,L)];
        e=norm(diff(q)/norm(q));
        if e<1E-6
           SS(L-2:L-1,L)=SS(L,L-2:L-1);
        else
           wblog('WRN','failed to fix SS for ladder with SS(L,L-2)=0 !?');
        end
     end
  end

  if nargout<=1, varargout={SS}; return; end

  s1=full(diag(SS,1)); lstr={};

  if gotlad && perBC, l=4; else l=2; end

  s2=full(diag(SS,l));
  if gotlad
     s1=s1(1:2:end); n=length(s2);

     lstr={ 'SS_{rung}', 'SS_{leg}' };

     if perBC
        sx=full(diag(SS,2));
        if s2(end)==0
           if s2(end-1)
              s2(end-1:end)=s2(end-1)/2;
           elseif sx(end-1)
              s2(end-1:end)=sx(end-1)/2;
           end
        end
        if all(sx([1 2])) && none(sx(3:end-2))
           s2=reshape([sx(1:2); s2],2,[])';
           s2=[flipud(s2(1:2:end,:)); s2(2:2:end,:)];
        end
     else
        if mod(n,l), n=n-mod(n,l); s2=s2(1:n); end
        s2=reshape(s2,l,[])';
     end
  elseif perBC
     s1=[0; flipud(s2(1:2:end,:)); s2(2:2:end,:)];
     s2=zeros(size(s1));
  end

end

% -------------------------------------------------------------------- %
function [sl,sr,sx,lstr]=getSS_HBLX(SS)

   q=size(SS,3); if q~=3
      wberr('unexpected size(SS,3)=%g/3 !?',q); end

   lstr={ '{\langle}SS_{leg}{\rangle}', 'SS_{rung}', 'SS_{X}' };

   q=SS(:,:,1); sr=diag(q); e(1)=norm(q-diag(sr));

   q=SS(:,:,2)/2; sl=diag(q,1);
   e(2)=norm(triu(q,2))+norm(tril(q,-2))+norm(diag(q));

   q=SS(:,:,3); sx=diag(q,1);
   e(3)=norm(triu(q,2))+norm(tril(q,-2))+norm(diag(q));

   if any(e)
      wblog('WRN','got unexpected HBLX data !? (@ %.3g)',max(e)); 
   end
end

% -------------------------------------------------------------------- %

