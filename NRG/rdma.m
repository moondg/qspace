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

    setopts(odma,...
      'T?','NRho?','rhoNorm?','nlog?','mspec?','emin?','emax?',...
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

    nops=max(numel(op1),numel(op2));

    if ~isset('use_mem')
       use_mem=(ischar(NRG) || ~isfield(NRG,'AK'));
    elseif ischar(NRG) || ~isfield(NRG,'AK')
       wbdie('invalid usage (got use_mem with invalid NRG variable)');
    end

    if ~isset('splitops') || ~nops
       setopts(odma,'cflags?','zflags?');
       if ~use_mem
            [om,a0,Idma    ]=fdmNRG_QS(NRG,     op1,op2,Z0, odma{:});
       else [om,a0,Idma,NRG]=fdmNRG_QS(NRG,Inrg,op1,op2,Z0, odma{:});
       end
    else
       if ~isset('nostore')
          wbdie('splitops should only be used with nostore');
       end

       q=splitops;
       if ~iscell(q)
          if isequal(q,1), q=matcell(1:nops);
          elseif sum(q)>nops || any(q<1),  wbdie('invalid splitops');
             q=mat2cell(1:sum(q),1,q);
          end
       end

       [i,p]=sort([q{:}]);
       if any(diff(i)<1) || any(i<1 | i>nops), wbdie('invalid splitops'); end
       if numel(i)~=nops || isequal(p,1:nops), p=[]; end

       n=cellfun(@(x) numel(x), q);

       Iops={ 'iop', 'n', 'zflags','cflags'
               q, matcell(n), [],  [] };
         if isset('cflags'), Iops{2,3}=mat2cell(cflags,1,n); end
         if isset('zflags'), Iops{2,4}=mat2cell(zflags,1,n); end
       Iops=struct(Iops{:});

       m=numel(Iops); clear om a0 Idma
       for k=1:m, i=Iops(k).iop;
          wblog('==>','using splitops (%d/%d) -> [%s ]',k,m,sprintf(' %d',i)); 

          o={[], op2(i), Z0};
          if ~isempty(op1), o{1}=op1(i); end
          q=Iops(k).cflags; if ~isempty(q), o=[o, {'cflags',q}]; end
          q=Iops(k).zflags; if ~isempty(q), o=[o, {'zflags',q}]; end

          if ~use_mem
               [om{k},a0{k},Idma(k)    ]=fdmNRG_QS(NRG,     o{:}, odma{:});
          else [om{k},a0{k},Idma(k),NRG]=fdmNRG_QS(NRG,Inrg,o{:}, odma{:});
          end

          if k>1
             if isequal(om{k},om{1}), om{k}=[];
             else wbdie('splitops got different om return data'); end
          end
       end

       om=om{1}; n=numel(p);
       a0=cat(2,a0{:}); if n, a0=a0(:,p); end

       Idma(1).finished=Idma(end).finished;

       q=[Idma.reA0];
          if n, q=reshape(q,2,[]); q=q(:,p); q=reshape(q,1,[]); end
          Idma(1).reA0=q;
       q=[Idma.a4];
          if n, q=reshape(q,4,[]); q=q(:,p); q=reshape(q,2,[]); end
          Idma(1).a4=q;
       Idma(1).symfac=[Idma.symfac];

	   pp=[Idma.paras];
	   for f={'B','C','cflags','zflags'}, f=f{1};
		  q=cell(1,m); for i=1:m, q{i}=getfield(pp(i),f); end
		  q=cat(2,q{:}); if ~isempty(q) && ~isempty(p), q=q(:,p); end
		  pp(1)=setfield(pp(1),f,q);
	   end
       Idma=Idma(1); Idma.paras=pp(1);

       Idma.Iops=Iops;
       Idma.splitops=splitops;
       Idma.p=p;

       clear pp Iops
    end

  % ================================================================ %

    a=a0; a(find(isnan(a)))=0; a=sum(a); a={ vec2str(a), max(abs(a-1)) };
    fprintf(1,'\nSum raw spectral data: %s (%.3g)\n',a{:}); 
    fprintf(1,'See structure Idma for more info.\n\n');

  end

% -------------------------------------------------------------------- %
  if plotflag && ~isempty(a0), dma_plot, end
