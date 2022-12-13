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

% NB! leave QSP_NUM_THREADS out of this (see startup_numethreads.m)!
% rather focus on standard matlab environment
%  -> num_threads() in matlab environment, and
%  -> OMP_NUM_THREADS = MKL_NUM_THREADS

  vflag=1;
  getopt('init',varargin);
     if getopt('-q'), vflag=0;
     elseif getopt('-v'), vflag=vflag+1; end
     iflag=getopt('-i');
     fflag=getopt('-f');
  n=getopt('get_last',[]);

  n0=num_threads();
  if ~nargin, rval=n0; return; end

  I=mlinfo;

  if nargin==1 && isequal(n,'-s')
     rval=I.nthreads; rval.nthreads=n0;
     return
  end

  if numel(n)>2, n, wberr('invalid usage'); end

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
function n=num_threads(varargin)

  n=feature('numthreads');
  if nargin, feature('numthreads',varargin{:}); end

end

% -------------------------------------------------------------------- %

