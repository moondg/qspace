% if isfield(param,'JJ')
%    banner('TEST rerun setup2CKondo'); % to save running cluster job ;)
%    setup2CKondo
% end

  if isset('Inrg') && isfield(Inrg,'Lambda') % isset('NRG')
     Lambda=Inrg.Lambda;
     if ~isempty(FOUT), NRG=FOUT{end}; end
  else
     if exist('user_fout','var') && ~isempty(user_fout)
     NRG=user_fout; else NRG='./NRG/NRG'; end

     fn=[NRG '_info.mat'];
     load2(fn,'Lambda','param','E0');
     if ~exist('FC','var') || ~exist('Z','var')
        q=load2(fn,'ops'); structexp(q.ops,'FC','Z');
     end
     if ~isempty(whos('-FILE',fn,'zflags')), load2(fn,'zflags'), end
     if ~isempty(whos('-FILE',fn,'cflags')), load2(fn,'cflags'), end
     if ~isempty(whos('-FILE',fn,'Gamma' )), load2(fn,'Gamma' ), end
  end

  if ~exist('Lambda','var')
  disp('Need Lambda for rescaling.'); return; end

% -------------------------------------------------------------------- %

  if exist('noDMA')~=1, noDMA=0; end
  if ~exist('plotflag','var'), plotflag=1; end

  if ~noDMA

    odma=setopts('-','-partial?','-PARTIAL?');

    if isset('dmaVersion')
         odma{end+1}=dmaVersion;
    else odma{end+1}='fDM'; end

    setopts(odma,'T?','NRho?','zflags?','cflags?','rhoNorm?','nlog?',...
      'mspec?','emin?','emax?',...
      '-calcOps?','-nostore?','-locRho?','-calcRho?','-keven?','-kodd?');

    if ~exist('Z0','var'),
    Z0=Z; end

    if 0
      if ~isvar('op1'), op1=[]; end
      if ~isvar('op2'), op2=FC; end
    else
      if ~isvar('op1'), op1=[FC FN]; end
      if ~isvar('op2'), op2=[FC FC]; end
    end

    if ischar(NRG) || ~isfield(NRG,'AK')
         [om,a0,Idma    ]=fdmNRG_QS(NRG,     op1,op2,Z0, odma{:});
    else [om,a0,Idma,NRG]=fdmNRG_QS(NRG,Inrg,op1,op2,Z0, odma{:});
    end
    a=a0; a(find(isnan(a)))=0; a=sum(a); a={ vec2str(a), max(abs(a-1)) };
    fprintf(1,'\nSum raw spectral data: %s (%.3g)\n',a{:}); 
    fprintf(1,'See structure Idma for more info.\n\n');

  end

% -------------------------------------------------------------------- %
  if plotflag && ~isempty(a0), dma_plot, end
