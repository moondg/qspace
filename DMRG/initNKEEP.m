
% outsourced from MEX/tst_Hamilton1D
% Wb,Apr08,16

  if isset('NKEEP_USER'), NKEEP=NKEEP_USER;
  else
   % WRN! nk1=2 is more likely to get stuck in local minima !?
   % Wb,Apr17,18
     setdef('nk1',4,'nk2',10,'Nspeed',2.8);
     if nk2<nk1, nk2=nk1; end

     if Nspeed>0 && Nspeed<1
        n__=nk1:Nspeed:nk2;
     elseif floor(Nspeed)==1
        n__=nk1:1/3:nk2; l__=findmel(n__,10*(Nspeed-1),-1);
        n__(2:2:l__-1)=[];
     elseif floor(Nspeed)==2
        if isset('MBOND') && MBOND>=4, n__=nk1:nk2;
        else
           n__=nk1:0.5:nk2; l__=findmel(n__,10*(Nspeed-2),-1);
           n__(2:2:l__-1)=[];
        end
     else
        Nspeed, error('Wb:ERR','\n   ERR invalid Nspeed');
     end

     NKEEP=round(2.^n__);

     q=whos('Nsw*');
     if ~isempty(q), NKEEP=matcell(NKEEP); s__={};
        for i__=1:numel(q)
           l__=str2num(q(i__).name(4:end));
           if ~isempty(l__)
              s__{end+1}=eval(q(i__).name);
              NKEEP{l__}=repmat(NKEEP{l__},1,s__{end});
              s__{end}=sprintf('Nsw(%g)=%g',l__,s__{end});
           end
        end
        if ~isempty(s__), wblog(' * ','using %s',strhcat(s__,'-s',', ')); end
        NKEEP=[NKEEP{:}];
     end
     if isset('Nfin')
        NKEEP=[NKEEP, Nfin];
     end
  end

  if isempty(NKEEP)
     error('Wb:ERR','\n   ERR got empty NKEEP (%g,%g)',nk1,nk2); 
  end

  if isset('NK1')
       Nkeep=NK1;
  else Nkeep=NKEEP(1); % setdef('Nkeep',NKEEP(1));
  end

  fprintf(1,'\n>> NKEEP = [%g;%s ] (%g sweeps)\n\n',...
  Nkeep, sprintf(' %g',NKEEP), numel(NKEEP));

  clear l__ i__ n__ s__

