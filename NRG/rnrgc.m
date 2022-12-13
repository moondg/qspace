  
  i=which('rnrg','-all');
  if length(i)>1 , wblog('NB!','\Nrunning cpy of rnrg ...');
  disp(strvcat(i{:})); end

% if ~exist('wnrg','var'), wnrg='default'; end
  if ~exist('wnrg','var'), wnrg='default'; end

  switch wnrg % tag: WNRG

    case {'default'}
      setdef('Nkeep',256); % clear Nkeep
      setdef('N',50);

      if ~exist('U','var') || ~exist('Lambda','var')
         istr='default SIAM parameter set';
         wblog('<i>','Loading %s.',istr);

         U=0.12; epsd=-0.04; Gamma=0.01; if ~isvar('B'), B=0E-4; end
         Lambda=2;

       % setParam
      end

      setdef('T',1E-4);
      [FFLANCZOS,x,Ics]=oliveira_CS(Gamma,Lambda,'T',T);

    % wnrg=''; % take `otherwise' below next time
      setupSIAM;

    % NB! global offset important for actual time evolution
      H0 = H0 + ...
      complex(0,-Ics.gamma0) * length(ff) * QSpace(contractQS(A0,[1 3],A0,[1 3]));

      add2struct(param,T,Ics); clear Ics
      param.gamma0=param.Ics.gamma0;

      if length(op2)==3
      op1=[]; op2=FC; clear cflags zflags; end

    % locRho=1; store RHO for later analysis (e.g. using dmrho.m)

    case '' % do nothing

    % A0=[]; setupSIAM;
    % Nkeep=512; B=5*TK;

    otherwise % default / test parameter set
      disp(wnrg)
      error('Wb:ERR','Invalid wnrg');

  end % of WNRG

  if ~isempty(findstr(pwd,'/home/'))
  cto lma; end

  param.wnrg=wnrg;

  if exist('T','var') && T>0 && (~exist('N','var') || isempty(N))
  N=ceil(-2*log(T/100)/log(Lambda)); N=N+mod(N,2); end

  if ~exist('calcflag','var'), calcflag=1; end
  if ~exist('plotflag','var'), plotflag=1; end
  onrg={}; % 'ionly' 

  if exist('user_fout','var') && ischar(user_fout), fout=user_fout;
  else
     fout= './NRG/NRG';
   % NB! make sure job output data is not overwritten by other jobs!
     s=getenv('JOB_ID'); if ~isempty(s), fout=[fout '.' s]; end
     s=getenv('SGE_TASK_ID'); if ~isempty(s), fout=[fout '.' s]; end
     user_fout=fout; % rdma and other routines will look for it!
  end

if calcflag || exist('TST_RNRG','var') && TST_RNRG

  if exist('fout','var') && ~isempty(fout)
       FOUT={'fout',fout};
  else FOUT={}; end

  if exist('Nkeep','var')
     if Nkeep>1024 && isempty(FOUT)
     FOUT={'fout','/data/weichsel/Matlab/Data/NRG'}; end
     onrg = { onrg{:}, 'Nkeep', Nkeep };
  else
     onrg = { onrg{:}, 'Nkeep', 256 };
  end

  param.D=onrg{end}; % needed for nrg_header
  if exist('fout','var') && ~isempty(fout)
     if fout(1)=='/'
          param.nrgIO=[hostname '/' fout];
     else param.nrgIO=[hostname '/' pwd '/' fout]; end
  end

  setopts(onrg,'NKEEP?','stype?')

% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %

  if exist('gg','var') && exist('FL','var') && ~isempty(gg)
     add2struct(param,gg,FL);
     onrg={ gg, FL, onrg{:} };
  end

  if exist('TST_RNRG','var') && TST_RNRG
     wblog('NB!','parameter setup mode - return')
     clear TST_RNRG; return
  end

  if isempty(FOUT)
  [NRG,Inrg]=NRGWilsonQS_CS(H0,A0,Lambda,ff,FC,Z,onrg{:}        );  else
       Inrg =NRGWilsonQS_CS(H0,A0,Lambda,ff,FC,Z,onrg{:},FOUT{:});
  end

  EE=Inrg.EE;
  E0=Inrg.E0;

end % of calcflag

% --------------------------------------------------------------------- %
  if plotflag, nrg_plot; end
% --------------------------------------------------------------------- %

