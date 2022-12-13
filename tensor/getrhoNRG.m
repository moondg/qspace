function [R,Rf,I]=getrhoNRG(kk,varargin)
% function [R,Rf,I]=getrhoNRG(kk [,opts])
%
%    get local reduced density matrix from NRG/NRG data space for given k-set.
%    kk may also be one of the following: last, all
%
% Options
%
%   'T',..    effective temperature
%   'NRG',..  NRG data space ('./NRG/NRG')
%   'N',..    up to what length (default: full length)
%   '-v'      verbose flag
%
% Examples: [R,Rf,I]=getrhoNRG('all','NRG','NRG_2CK/NRG','-v');
%
% See also getrhoNRG_red.m
% Wb,Jun01,08

% outsourced from nrgtangle.m

  getopt('INIT',varargin);
     T     =getopt('T',[]);
     vflag =getopt('-v'); if ~vflag && getopt('-V'); vflag=2; end
     nrg   =getopt('NRG','./NRG/NRG');
     N     =getopt('N',[]);
  getopt('check_error');

  if isempty(findstr(pwd,'Data')), cto lma, end
  ff=dir2([nrg '_[0-9]*.mat']);

  if isempty(ff)
  error('Wb:ERR','\n  ERR no files %s* found\n',nrg); end

  load([nrg '_info'],'Lambda','param');

  [x,i]=sortrows(strvcat(ff.name)); ff=ff(i);

  if isempty(N) || N>=length(ff), N=length(ff);
     load(ff(N).name); A=AT; H=HT;
     if ~isempty(QSpace(HK)) error('Wb:ERR',...
     'HK must be empty at last iteration !??'); end
  else
     load(ff(N).name); A=AK; H=HK;
  end

  if ischar(kk)
     switch kk
        case {'last','end'}, kk=N;
        case 'all',  kk=1:N;
        otherwise error('Wb:ERR','invalid k specs.');
     end
  elseif any(kk>N), error('Wb:ERR','invalid k-index set'); end

  if ~isempty(T)
       beta=( Lambda^(-N/2) * (Lambda+1)/2 )/T;
  else beta=100;
  end

  [R,I]=getrhoQS(H,beta); R=skipzeros(QSpace(R)); I.RN=R;
  X=R; se=nan(1,N);

% -------------------------------------------------------------------- %
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
     wblog('WRN','local QSpace changes (%g) !??',k); end

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

