function rkondo(varargin)
% Function rkondo([NC,Nkeep])

if nargin==1 && isfile2(varargin{1},'-mat')
   load2(varargin{1});
else

  getopt('init',varargin);
     kflag   = getopt('-k');
     tflag   = getopt('-t');
     noT     = getopt('-t0');
     chiflag = getopt('-chi');
     splitops= getopt('-split');
  varargin=getopt('get_remaining');

  onefile=1;

  Imain = struct(...
     'start',  datestr(now), ...
     'args',   {varargin}, ...
     'pwd',    pwd, ...
     'host',   hostname, ...
     'file',   mfilename, ...
     'JOB_ID', [getenv('JOB_ID') '.' getenv('TASK_ID')] ...
  );

% -------------------------------------------------------------------- %

  wsys='KondoAM';

  U=0; epsd=-U/2; Gamma=0.01; B=0; Lambda=2; J2=-0.008;

  [NC,i2,ftag]=dealargs(varargin{:});
  if isempty(NC), NC=1; end

  nostore=1;
  locRho=1;

  switch NC
    case 1
      Nkeep=512; N=50; nlog=512;
      J=2.5*J2/(NC+0.5);

      sigma=0.4; tfac=1/4; chiflag=1; clear nostore

      if 1
         Nkeep=1024; Lambda=1.7; N=90;
         TT=0.01; wflag=1; tfac=1;
         TT=-1; clear wflag
         N=50; Lambda=2; Nkeep=512; sigma=0.6; nlog=128;
      else
         Nkeep=512; Lambda=2; N=50;
         TT=-1;
      end

    case 2
    % Nkeep=6144; N=50; % frankie got stuck (requires > 4G)
    % Nkeep=5000; N=50;
    % J=-0.006; % T=1E-6;
    % Nkeep=6000; N=70; Lambda=1.7;
    %
      Nkeep=4096; N=50; chiflag=1;

      J=J2;

      TT=logspace(-8,-2,32); sigma=0.5; tfac=1/4; chiflag=1; Lambda=2.5;

    case 3
      Lambda=2.5; N=40; Nkeep=4500; splitops=1;
      J=2.5*J2/(NC+0.5);

      TT=logspace(-8,-4,10); sigma=0.6; tfac=1; Lambda=3;

    otherwise
      NC, wbdie('invalid NC'); 
  end

  if ~isempty(i2), Nkeep=i2; end

  if noT, TT=-1; end

  cto lma

  tstflag=1;

  i=getenv('JOB_ID');
  if isbatch && ~isempty(i)
     user_fout=['NRG/nrg' i];
  end

  o={'-f'};

  if isbatch
     tstflag=0;
     o{end+1}='-nodisp';
  end

  if (~exist('ftag','var')  || isempty(ftag)) && B~=0, ftag='b'; end

  if tstflag, mat=mfilename;
  else
     mat=sprintf('%s_rkondo%g_J%02.0f_G%02.0f_D%d',...
     wbtstamp, NC, -1000*J, 100*Gamma, Nkeep);
     if exist('ftag','var') && ~isempty(ftag), mat=[mat ftag]; end
  end
  fname=mat;

  fprintf(1,'\n   PWD: %s/%s\n   MAT: %s\n\n',hostname,pwd,mat);
  disp(add2struct('-',...
     'wsys?','NC?','Nkeep?','calcOps?','nostore?','locRho?',...
     'kflag?','chiflag?','splitops?' ...
  ));

if tflag, keyboard, return; end

  if 0 && tstflag
     rnrg; save([mat '.mat']);
     rdma; save([mat '.mat']);
     return;
  end

end

  if ~exist('Inrg','var') || exist('recalc','var') && recalc
     if NC>1, sysinfo, end
     rnrg; save([mat '.mat']); mpdf(o{:},[mat '_eflow']);
     if NC>1, sysinfo, end
  end

  if ~exist('tfac','var'), tfac=0.5; end

  if ~chiflag
     TZ=Lambda^(-(N-10)/2);
     TZ=linspace(TZ/3,TZ,3); nT=length(TZ);

     OP1=op1; op1=SZ; zflags=[0];
     OP2=op2; op2=SZ; cflags=[1];

     for it=1:nT, T=TZ(it); eps=tfac*T;
        fprintf(1,'\n>> magn. suszeptibility: it=%d/%d ...\n\n',it,nT);

        rdma

        TK = 1./((4/4)*sum(Idma.reA0(end-1:end)));

        ZDMA(it)=Idma;
        ZALL(it)=add2struct('-',om,a0,ox,ax,TK);

        save([fname '.mat']);
        if NC>1, sysinfo, end
     end

     swap(op1,OP1);
     swap(op2,OP2); clear zflags cflags

  end

  CFLAGS={}; ZFLAGS={};

  if splitops
     if isempty(op1), op1=op2; end
     m=numel(op2); em=ones(m,1);
     OPS={ mat2cell(op1(:),em,1), mat2cell(op2(:),em,1) };
     OP1=op1; OP2=op2;

     if exist('cflags','var'), CFLAGS=mat2cell(cflags,1,em); end
     if exist('zflags','var'), ZFLAGS=mat2cell(zflags,1,em); end
  else
     if isempty(op1), op1=op2; end
     OPS={ {op1}, {op2} };

     if exist('cflags','var'), CFLAGS={cflags}; end
     if exist('zflags','var'), ZFLAGS={zflags}; end
  end

  nops=size(OPS{1},1);

  if ~exist('TT','var') || isempty(TT), TT=-1; end; nT=length(TT);
  for it=1:nT

     T=TT(it); if T<0, clear T eps
     else
        if nT>1 && ~onefile, fname=[mat sprintf('_%02d',it)]; end
        eps=tfac*T;
     end

     for iop=1:nops
        if nT>1 || nops>1
        fprintf(1,'\n>> it=%d/%d  ops=%g/%g...\n\n',it,nT,iop,nops); end

        if ~isempty(CFLAGS), cflags=CFLAGS{iop}; end
        if ~isempty(ZFLAGS), zflags=ZFLAGS{iop}; end

        op1=OPS{1}{iop};
        op2=OPS{2}{iop};

        rdma, if ~exist('T','var'), T=Idma.T; end

        IDMA(it,iop)=Idma;
        DOPS(it,iop)=add2struct('-',om,a0,ox,ax,op1,op2);

        if iop>1, if ~isequal(om,DOPS(it,1).om) || ~isequal(ox,DOPS(it,1).ox)
            wblog('ERR','om or ox mismatch !??')
            whos om ox, DOPS(it,1)
        end, end

        if iop==nops
           wblog('<I>','\NCOLLECT SPECTRAL DATA\N');
           a0=cat(2,DOPS(it,:).a0);
           ax=cat(2,DOPS(it,:).ax);
        else
           save([fname '.mat']);
           continue
        end

        gc=getGC_JAM(ox, ax(:,1), ax(:,2));
        ax=[ ax(:,1:2), imag(gc), ax(:,3:end) ];

        [itp,iTp]=getDephasing(ox,ax(:,1), Gamma,T);
        [itc,iTc]=getDephasing(ox,gc,Gamma,T);

        if chiflag
           TK = 1./((4/4)*sum(Idma.reA0(end-1:end)));
        end

        IALL(it)=add2struct('-',om,a0,ox,ax,gc,...
        'it: [itp,itc]','iT: [iTp iTc]','TK?');

        if it==1 || ~onefile, mpdf(o{:},[fname sprintf('_a%g',iop-1)]); end

        save([fname '.mat']);
        if NC>1, sysinfo, end

     end
  end

  inl 1; disp(IALL(it).iT.')

% -------------------------------------------------------------------- %
  if isbatch % exist('user_fout','var')
     if kflag, o={}; else o={'rm',user_fout}; end
     sge_finish('fout',[mat '*'], o{:});
  end

end

