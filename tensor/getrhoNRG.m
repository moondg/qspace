function [R,Rf,I]=getrhoNRG(kk,varargin)
% function [R,Rf,I]=getrhoNRG(kk [,opts])
%
%    get reduced density matrices R in local state space (s)
%    from NRG/NRG data space for specified index set kk
%    which alsomay  be one of the following: 'last', [] = 'all'
%
% Options
%
%   'NRG',..  NRG data space ('./NRG/NRG')
%   'N',..    up to what length (default: full length)
%
% Usage #2: [RR,IR]=getrhoNRG(Inrg.HK[,opts])
%
%    Compute `reduced density matrices' ~ exp(-H/T) based on HK.
%
% Options for any usage
%
%   'T',..    effective temperature
%   '-v'      verbose flag
%
% Examples: [R,Rf,I]=getrhoNRG('all','NRG','NRG_2CK/NRG','-v');
%
% Wb,Jun01,08

% See also getrhoNRG_red.m
% outsourced from nrgtangle.m

  usage_2=isa(kk,'QSpace') || (isfield(kk,'Q') && isfield(kk,'data'));

  getopt('INIT',varargin);
     T     =getopt('T',[]);
     vflag =getopt('-v'); if ~vflag && getopt('-V'); vflag=2; end
     if ~usage_2
        nrg=getopt('NRG','./NRG/NRG');
        N  =getopt('N',[]);
     end
  getopt('check_error');

  if usage_2
    [R,Rf]=getrho_HK(kk,T,vflag); I=[];
    return
  end

  if isempty(findstr(pwd,'Data')), cto lma, end
  ff=dir2([nrg '_[0-9]*.mat']);

  if isempty(ff), wbdie('no files %s* found\n',nrg); end

  load([nrg '_info'],'Lambda','param');

  [x,i]=sortrows(strvcat(ff.name)); ff=ff(i);

  if isempty(N) || N>=length(ff), N=length(ff);
     load(ff(N).name); A=AT; H=HT;
     if ~isempty(QSpace(HK))
        wbdie('HK must be empty at last iteration !?');
     end
  else
     load(ff(N).name); A=AK; H=HK;
  end

  if isempty(kk), kk=1:N;
  elseif ischar(kk)
     switch kk
        case {'last','end'}, kk=N;
        case 'all',  kk=1:N;
        otherwise wbdie('invalid kk input');
     end
  elseif any(kk>N), wbdie('invalid kk-index set'); end

  if ~isempty(T)
       beta=( Lambda^(-N/2) * (Lambda+1)/2 )/T;
  else beta=100;
  end

  [R,I]=getrhoQS(H,beta); R=skipzeros(QSpace(R)); I.RN=R;
  X=R; se=nan(1,N);

% -------------------------------------------------------------------- %
% build reduced density matrix space

  R=QSpace(1,N); Rf=cell(1,N);
  kmin=max(1,min(kk));

  if vflag, inl(1); end

  P=PSet('k',N:-1:kmin);
  for ip=1:P.n, [p,pstr,tstr]=P(ip); structexp(p);

     if vflag
        if vflag>1, tstr=regexprep(tstr,'estimated.*finished: ','');
             fprintf(1,'   %s  %s: %s  (%s) \r',time('-t'),nrg,pstr,tstr);
        else fprintf(1,'   %s  %s_%02g/%g ...\r',time('-t'),nrg,k,N);
        end
     end
     q=contractQS(A,3,X,2);

     R(k)=contractQS(A,[2 3], q, [2 3]);
     X   =contractQS(A,[1 3], q, [1 3]);

     [Rf{k},i]=mpsFull2QS(R(k));
     se(k)=SEntropy(mpsFull2QS(X));

     if k==N, iN=i; elseif ~isequal(i,iN)
        wblog('WRN','local QSpace changes (%g)',k); end

     if k<2, break; end
     load(ff(k-1).name,'AK'); A=AK;
  end
  if vflag, inl(2); end

  if k==1, I.R0=X; end
  I.se=se;

  if length(kk)>1
       R=R(kk); Rf=Rf(kk);
  else R=R(kk); Rf=Rf{kk}; end

end

% -------------------------------------------------------------------- %
% get `reduced density matrices' based on Inrg.HK data
% e.g., may use this to track changes in energy flow diagram
% the points of strongest change are regurned in Iout.kpeak_ds
% Wb,Feb26,24

function [RR,Iout]=getrho_HK(HK,T,vflag)

   if isempty(T), beta=1; else beta=1/T; end
   if ~isfinite(beta) || beta<0, wbdie('invalid usage (T=%g)',T); end

   L=numel(HK); if isempty(QSpace(HK(L))), L=L-1; end
   RR=QSpace(1,L); 
   Iflag=nargout>1; if Iflag, se=zeros(1,L); end

   for k=1:L, Hk=QSpace(HK(k));
      if ~Hk, wbdie('got empty HK at k=%d/%d',k,L); end
      if k==1 && isdiag(Hk)<2
         [ee,Ie]=eigQS(Hk);
         Hk=QSpace(Ie.EK)-ee(1);
      end

      RR(k)=getrhoQS(Hk,beta); if ~Iflag, continue; end
      se(k)=SEntropy(RR(k));
   end

   if Iflag
      ds=abs(diff2(avgdata(se,2,'-l'),'len'));
      [xp,yp,Ip]=findpeak(1:L,ds,'--max'); n=numel(xp);

      Iout=add2struct('-',beta,se,ds);

      Iout.kpeak_ds=xp; if ~n, xp=nan; end
      Iout.k0=xp(end);
   end

end

% -------------------------------------------------------------------- %

