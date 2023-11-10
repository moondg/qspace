  
  i=which('rnrg','-all');
  if length(i)>1 , wblog('NB!','\Nrunning cpy of rnrg ...');
  disp(strvcat(i{:})); end

  setdef('wsys','SIAM');

  switch wsys % tag: WNRG

    case 1.0 % Krishna Murthy p.1016

      U=1E-3; epsd=-U/2; Gamma=U/(12.66*pi); B=0;
      Lambda=2.5; N=75; Nkeep=1024;

      A0=[]; setupSIAM;

    case 'SIAM'

      setdef('U',0.12); setdef('epsd',-U/2); SYM='A,A';
      setupSIAM_SU2x2;
      setdef('Etrunc',5,'Nkeep',600); % clear Nkeep

    case 'TICM' % Wb,Feb02,21
    % use_mem=1;
      setupTICM; locRho=1; nostore=1;
      setdef('Etrunc',5,'Nkeep',2048);

    case {'SIAM-SU2x2'} % Wb,Sep20,11

      setupSIAM_SU2x2; locRho=1; % nostore=1;
      setdef('Etrunc',5,'Nkeep',600); % clear Nkeep

    case {'Kondo-SU2x2','KondoJ-SU2x2'} % Wb,Sep19,11
      setdef('J',0.12);

      if isset('SYM'), s=[' [' SYM ']']; else s=''; end

      istr=sprintf('Kondo model (J=%g; using %s%s)',J,wsys,s);
      wblog('<i>','%s',istr);

      if ~exist('B','var'), B=0; end
      if ~exist('Nkeep','var'), Nkeep=256; end
      if ~exist('N','var'), N=51; end % double EPS reached at N=100
      if ~exist('Lambda','var'), Lambda=2; end
      
      locRho=1; % nostore=1;
      setupKondo_SU2x2

    case {'2CKondo'}
      istr='(2-channel) Kondo';
      wblog('<i>','%s',istr);

    % afac=[pi^2/2, 1] % see r2CK
    % see setupKondo_AM
    % NC=1;

      locRho=1; % store RHO for later analysis (e.g. using dmrho.m)
      setup2CKondo;

    case {'KondoJH'}
    % Wb,May19,08
      istr='two-channel Anderson model with JH (A.Rosch)';
      wblog('<i>','%s',istr);

    % locRho=1; store RHO for later analysis (e.g. using dmrho.m)
      setupAM_JJH

    case 'pg' % pseud-gap SIAM with symmetric standard parameter set

      U=0.12; epsd=-U/2; Gamma=0.01; if ~exist('B','var'), B=0; end
      Lambda=2; N=50; Nkeep=512;

      if ~exist('pg_r','var'), pg_r=0.3; end
      
      A0=[]; locRho=1;
      setupSIAM

      if length(op2)==3
      op1=[]; op2=op2(2:end); clear zflags; end

    case {'2HAM','2HKM'} % 2-imppurity helical Anderson model // Wb,Aug11,21
      if isequal(wsys,'2HAM'), q=1; else q=0; end

      setdef('Uflag',q,'use_mem',1,'locRho',1,'nostore',1);
    % setupHelical_old; NKEEP=4096;
      setupHelical

    % Lambda=4
    % - Etrunc=6 => stays below Nkeep ~ 2000
    % - Etrunc=8 => reaches Nkeep>4096 (~5000?) at iteration k=2 only)
      setdef('Etrunc',6,'Nkeep',4096);

    case 'pgk' % pseud-gap Kondo

      Lambda=2.0; N=70; Nkeep=256;

      if ~exist('B','var'), B=0; end
      if ~exist('J','var'), J=0.0666; end % TK=2E-4
      if ~exist('pg_r','var'), pg_r=0.3; end

      if ~exist('T','var'), T=1E-6; end
      eps=T*1E-2; emin=T*1E-4;

      locRho=1;
      setupKondo

    case {'SIAM2-off-diag'}
    % test setting for off-diagonal cross couplings between
    % the two flavors (`spins') along Wilson chain
    % therefore using Acharge only // Wb,Feb01,21

      SYM='Acharge'; B=1E-12; use_mem=1;
      setupSIAM_SU2x2; locRho=1; nostore=1;
      setdef('Etrunc',5,'Nkeep',600);

      setdef('ff_offdiag',[1E-4 0.05i]);

      ff_=ff; ff=reshape(ff,[],1);
      ff(:,2)=ff_offdiag(2)*ff; ff(1,2)=0;
      i=find(abs(ff(:,1))>abs(ff_offdiag(1))); ff(i,2)=ff_offdiag(1);

      ff=reshape([ff, fliplr(ff)],[N-1,2,2]); ff_={ff_,ff};
      ff=ff(:,:);

    case {'TICM'} % 2-imppurity chiral model // Wb,Feb06,21

      use_mem=1;
      setupTICM; locRho=1; nostore=1;
      setdef('Etrunc',5,'Nkeep',2048);

    case {'','user'}
    % assume parameters are already userdefined => do nothing

    otherwise % default / test parameter set

      if exist(wsys,'file'), run(wsys);
      else
         disp(wsys), error('Wb:ERR','Invalid wsys');
      end
  end % of WNRG

  if isset('CHECK_NKEEP')
     if ~exist('Ic','var') || ~isfield(Ic,'kx')
         error('Wb:ERR',['\n   ERR CHECK_NKEEP requires return value Ic ' ...
        'from getNRGcoupling()']);
     end
   % add safety margin of 1 iteration if allowed by CHECK_NKEEP
     kx=min(Ic.kx+1,CHECK_NKEEP);
     NKEEP=4^kx; if NKEEP>Nkeep
     NKEEP=repmat(NKEEP,1,kx); else clear NKEEP; end
  end

  if isempty(findstr(pwd,'/data/'))
  cto lma; end

  add2struct(param,wsys,'Etrunc?','ET1?','ETRUNC?','NKEEP?');

  if exist('T','var') && T>0 && (~exist('N','var') || isempty(N))
  N=ceil(-2*log(T/100)/log(Lambda)); N=N+mod(N,2); end

  if ~exist('calcflag','var'), calcflag=1; end
  if ~exist('plotflag','var'), plotflag=1; end
  onrg={}; % 'ionly' 

  if exist('user_fout','var') && ischar(user_fout), fout=user_fout;
  else
     if isempty(getenv('ML_DEBUG'))
          fout= './NRG/NRG';
     else fout= './NRG/TST';
     end

   % NB! make sure job output data is not overwritten by other jobs!
     s=getenv('JOB_ID'); if ~isempty(s), fout=[fout '.' s]; end
     s=getenv('SGE_TASK_ID');
       if ~isempty(s) && isempty(findstr(s,'undef'))
       fout=[fout '.' s]; end
     user_fout=fout; % rdma and other routines will look for it!
  end

if calcflag || isset('TST_RNRG')

  if exist('fout','var') && ~isempty(fout) && ~isset('use_mem')
       FOUT={'fout',fout};
  else FOUT={}; end

  if exist('Nkeep','var')
   % if Nkeep>1024 && isempty(FOUT)
   % FOUT={'fout',[getenv('LMA') '/NRG']}; end
     onrg=[ onrg, {'Nkeep',Nkeep} ];
     if Nkeep>9999 && isset('fflag'), onrg{end+1}='-f'; end % Wb,Nov07,23
  else onrg=[ onrg, {'Nkeep',256} ];
  end

  param.D=onrg{end}; % needed for nrg_header
  if exist('fout','var') && ~isempty(fout)
     if fout(1)=='/'
          param.nrgIO=[hostname '/' fout];
     else param.nrgIO=[hostname '/' pwd '/' fout]; end
  end

  if exist('NKEEP','var')
  onrg(end+1:end+2)={'NKEEP', NKEEP}; end

  if exist('ZFLAG') && ~isempty(ZFLAG)
   % required for particle-hole symmetry
     onrg(end+1:end+2)={'zflag',ZFLAG};
  end

  if isset('nrg_qflag'), onrg{end+1}='-q'; end

% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %

  if exist('deps','var') && deps>1E-6, wblog('WRN','deps=%g !??',deps); end
  if exist('db','var') && db>1E-4, wblog('WRN','db=%g !??',db); end

  setopts(onrg,'deps?','db?','dmax?','ETRUNC?')
  if exist('Etrunc','var') && Etrunc>0, setopts(onrg,Etrunc); end
  if isset('Estop')
     if Estop~=1, setopts(onrg,Estop);
     else setopts(onrg,'--Estop'); end
  end

  if exist('gg','var') && exist('FL','var') && ~isempty(gg)
     add2struct(param,gg,FL);
     onrg={ gg, FL, onrg{:} };
  end

  if isset('TST_RNRG')
     wblog('NB!','parameter setup mode - return')
     clear TST_RNRG; return
  end

% NB! FC is the preferred operator within rdma.m
% ==> keep FF as the preferred operator name for Wilson sites
% Wb,Sep28,13
  if     exist('FF','var') && ~isempty(FF), o={FF(:),Z,onrg{:}};
  elseif exist('FC','var') && ~isempty(FC), o={FC(:),Z,onrg{:}};
  else error('Wb:ERR',...
      '\n   ERR operators for Wilson chain undefined');
  end

% NRGWilsonQS => NRGWilsonQS
% NB! NRGWilsonQS/fdmNRG_QS also switches to col-major / LRs convention!
% Wb,Aug11,10

  if isempty(FOUT)
     [NRG,Inrg]=NRGWilsonQS(H0,A0,Lambda,ff,o{:});
  else    Inrg =NRGWilsonQS(H0,A0,Lambda,ff,o{:},FOUT{:});
  end

  if isset('Estop'), N=Inrg.N; param.N=N; end

  EE=Inrg.EE;
  E0=Inrg.E0;

end % of calcflag

% --------------------------------------------------------------------- %
  if plotflag, nrg_plot; end
% --------------------------------------------------------------------- %

