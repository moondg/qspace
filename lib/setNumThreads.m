function rval=setNumThreads(varargin)
% function setNumThreads(n [,opts])
%
%    (Re)set number of maximum active computational threads
%    where n can contain one or two numbers:
%    this sets feature('numthreads',n(1)), as well as
%    OMP_NUM_THREADS = MKL_NUM_THREADS = n(end) if n(end)>0;
%    it unsets the latter if n(end)<1.
%
% Wb,Apr05,16

% NB! leave QSP_NUM_THREADS out of this by default (see startup_numthreads.m)!
% rather focus on standard matlab environment
%  -> num_threads() in matlab environment, and
%  -> OMP_NUM_THREADS = MKL_NUM_THREADS

  vflag=1;
  getopt('init',varargin);
     if getopt('-q'), vflag=0;
     elseif getopt('-v'), vflag=vflag+1; end
     iflag=getopt('-i');
     fflag=getopt('-f');
     QSflag=getopt('--qs');
  n=getopt('get_last',[]);

  n0=num_threads();
  if ~nargin, rval=n0; return; end

  I=mlinfo;

  if nargin==1 && isequal(n,'-s')
     rval=I.nthreads; rval.nthreads=n0;
     return
  end

  if ischar(n), n=str2num(n); end
  if numel(n)>2, n, wbdie('invalid usage'); end

  if QSflag, qsp='QSP_NUM_THREADS';
     n2=n;
           if ~ischar(n2), n2=num2str(n2); end
           n2_=n2; if isempty(n2), n2_=''''''; end
     n_=getenv(qsp);
     if ~nargout || vflag>1, sout={'',''};
        if isequal(n_,n2), sout{2}=sprintf('already at %s = %s',qsp,n2);
        else
           if isempty(n_)
                sout{1}=sprintf('%s -> %s',qsp,n2_);
           else sout{1}=sprintf('%s=%s -> %s',qsp,n_,n2_); end
           setenv(qsp,n2);
        end
        fprintf(1,['\n   %s' char(27) '[38;5;8m' ... 
          '%s  (having numthreads=%d)' char(27) '[0m\n\n'],sout{:},n0);
     end
     if nargout, rval=n0; end
     return
  end

  nc=I.nthreads.cores;
  n1=max(n);

  if n(end)>=1
     if norm(diff(n)), s2=sprintf('%g/%g',n); else s2=num2str(n(1)); end
     sk=num2str(n(end));
  else
     s2=sprintf('%g/0',n(1));
     sk='';
  end

  q=getenv('MKL_NUM_THREADS');
  if isempty(q)
       s0=sprintf('%g/0',n0);
  else s0=sprintf('%g/%s',n0,q); end

  if n1>nc
     wblog('ERR',['\Nnum_threads exceeds number of cores (%s/%g)\n' ... 
     'keeping current numthreads=%s\N'],s2,nc,s0);
     return
  elseif ~fflag && n1==n0 && isequal(q,sk), if vflag, wblog(' * ',...
     'having num_threads %s (%g cores)',s2,nc); end
     return
  end

  sout=sprintf('setting num_threads %s -> %s (@ %g cores)',s0,s2,nc);

  num_threads(n(1));

  setenv('OMP_NUM_THREADS',sk);
  setenv('MKL_NUM_THREADS',sk);
  v='MKL_DOMAIN_NUM_THREADS';
  if ~isempty(sk) || ~isempty(getenv(v)), setenv(v,sk); end

  if iflag, fprintf(1,'%s\n',sout);
  elseif vflag, wblog('NB!','%s',sout);
  end

  if vflag>1, system(['set | egrep -v ''^BASH_EXEC'' | ' ...
    'egrep ''MKL|OMP[^A-Z]|THREADS|NSLOT''']);
  end

end

% -------------------------------------------------------------------- %
% see MEX/tst_numthreads.m // Wb,Mar30,16
% warning off MATLAB:maxNumCompThreads:Deprecated
% n=maxNumCompThreads(NTH);
% warning on MATLAB:maxNumCompThreads:Deprecated

function n=num_threads(varargin)

  n=feature('numthreads');
  if nargin, feature('numthreads',varargin{:}); end

end

% -------------------------------------------------------------------- %

