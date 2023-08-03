function se = SEntropy(r,varargin)
% function se = SEntropy(sv [,opts])
%
%    Calculate entanglement entropy for given input vector
%    The input sv is considered singular values by default.
%
% Options
%
%    '-q'     quiet
%    '-n'     explicitely normalize to 1. (assumed otherwise)
%    '-rho'   expect eigennvalues of reduced density matrix
%             (by default, SVDs are assumed otherwise)
%
% Wb,Oct19,05  Wb,Aug05,06

  getopt('init',varargin);
     qflag  = getopt('-q');
     nflag  = getopt({'-n','norm'});
     gotsvd =~getopt({'-rho','rho'});
  getopt('check_error');

  if isempty(r), se=0; return; end

  [n,m]=size(r); d=[];
  if m==2 && all(r(:,2)==round(r(:,2)))
     d=r(:,2); r=r(:,1); if gotsvd, r=r.*conj(r); gotsvd=0; end
     nrm=d'*r;
     if nflag, r=r/nrm; nflag=0;
     else e=norm(1-nrm);
        if e>1E-8 || any(d~=round(d))
           wberr('got non-normalized (degenarate) svd spectrum (e=%g)',e);
        end
     end
  elseif ~isvector(r)
     e=(r-r').^2; e=sum(e(:));
     if e/max(abs(r(:)))>1E-14
        wblog('ERR','invalid usage (size %s)',vec2str(size(r),'sep','x'));
        se=nan; if ~isbatch, keyboard, end, return
     end
     r=0.5*eig(r+r'); gotsvd=0;
     if ~isreal(r), wblog('WRN','eigenvalues not real !??'); end
  end
  if ~qflag
     if min(r)<-1E-14, wblog('WRN',...
        'negative entries for shannon entropy ?? (%g)', min(r));
     if isempty(d), e=sum(r); if abs(e-1)>1E-10
        wblog('WRN','rho_i data not normalized (1%+.3g)',e-1); end
     end
  end

  r=r(:);
  if gotsvd, r=r.*conj(r); end
  if nflag, r=r/sum(r); end
  i=find(r<=0); r(i)=[];

  if ~isempty(d), d(i)=[];
       se = d' * (-r.*log(r));
  else se = -r'*log(r);
  end

end

