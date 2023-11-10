function [Iout,Inrg]=specificHeatNRG(varargin)
% function [Iout,Inrg]=specificHeatNRG([Inrg|NRG,][opts])
%
% Options
%
%   'TT',...   set of temperatures (default: along Wilson shells)
%   'ff',...   uses these ff for analytic comparison (default: Inrg.ff)
%   '-q'       quiet mode
%   '-b'       include exact bath contribution (determined from ff)
%   '-p'       summarize results in plot
%   'ah',...   specify specific plot handles (at least 2x2 axes set)
%
% Wb,Jun10,09

  global param

  dflag=0;
  if nargin && (isstruct(varargin{1}) || isempty(varargin{1}))
     Inrg=varargin{1}; varargin=varargin(2:end);

     if ~isempty(Inrg) && isfield(Inrg,'IE') && isfield(Inrg.IE,'ED')
        dflag=1;
     end
  end

  getopt('init',varargin);
     TT    = getopt('TT',[]);
     ff    = getopt('ff',[]);
     qflag = getopt('-q');
     bflag = getopt('-b');
     pflag = getopt('-p');
     ah    = getopt('ah',[]); if ~isempty(ah), pflag=1; end

  if dflag, getopt('check_error');
     EK=Inrg.IE.EK;
     ED=Inrg.IE.ED; nrg='';
  else
     nrg=getopt('get_last',[getenv('LMA'), '/NRG/NRG']);

     if ~ischar(nrg), wbdie('invalid usage'); end
     [EK,ED,Inrg]=getEDdata(nrg,'-m','-d','-K');
     Inrg.IE=add2struct('-',EK,ED);
     Inrg.nrg=nrg;
  end

  if (bflag  || pflag) && isempty(ff)
     ff=Inrg.ff(:,1).*Inrg.EScale(2:end)';
  end

  N=Inrg.N; Lambda=Inrg.Lambda;
  if isempty(TT), TT=Lambda.^(-(N/2:-.25:-5)); end

  gotM=any(ED(:,3)>1);

  nT=numel(TT); C=nan(size(TT)); S=nan(size(TT)); dk=1;
  for it=1:nT, T=TT(it); if ~qflag
     fprintf(1,'\r   temperature %3g/%g ...  \r',it,nT); end

     k=ED(:,1); E=ED(:,2); xx=E/T;

   % xx=xx*0.988; % shifts S to towards 0, yet shifts C away from 0!
   % NB! slight shift towards lower energies shifts entropy
   % ==> comparison to exact result for bath has artifacts
   % ==> systematic shift of ED to higher energies within FDM
   %     due to missing discarded space from earlier iterations(!)
   % hint: compare within NRG (ie. one extra run with the dot
   % eliminated, e.g. by having the ground state at the impurity
   % split-off by a large energy -- Wb,Aug12,11
     zz=exp(-xx + (Inrg.N-k-dk)*log(Inrg.dloc));
     if gotM, zz=zz.*ED(:,3); end

     Z=sum(zz); zz=zz*(1/Z);
     E1=zz'*E;
     E2=zz'*(E.^2);
     C(it)=(E2-E1^2)/T^2;
     S(it)=E1/T + log(Z);

  end
  if ~qflag, fprintf(1,'\r%50s\r',''); end

  Iout=add2struct('-','T=TT',C,S);

if isempty(ff), return; end

  nb=log(Inrg.dloc)/log(2);

  [Cb,Sb,Eb]=specificHeatTB(ff(2:end),TT,'nb',nb);

  add2struct(Iout,Cb,Sb,Eb,ff,nb);

if ~pflag, return; end
if isempty(ah)
   ah=smaxis(2,2,'tag',mfilename);
   header('%M :: %s',nrg); nrg_header; addt2fig wb
elseif any(size(ah)<[2 2])
   wbdie('invalid ah (at least 2x2 required)');
end

setax(ah(1,1))
  bs=sprintf('exact bath^{\\otimes%g}',nb);

  semilogx(TT,C,'o-','Disp','FDM'); hold on
  semilogx(TT,Cb,'r','Disp',bs); hold on

  xtight; sms(2); label('T','total specific heat C');
  legdisp('Location','SouthWest');

setax(ah(1,2))

  semilogx(TT,S,'o-','Disp','FDM'); hold on
  semilogx(TT,Sb,'r','Disp',bs); hold on

  xtight; sms(2); label('T','total entropy S');
  legdisp('Location','SouthEast');

setax(ah(2,1))

  semilogx(TT,C-Cb,'o-'); xtight; sms(2); hold on
  label('T','relative C (impurity contribution)');

  if isfield(Inrg.param,'U')
       q=[Inrg.param.U, Inrg.param.Gamma, TKondo(Inrg.param)];
  else q=nan(1,3); end
  xmark(q);

setax(ah(2,2))

  loglog(TT,exp(S-Sb),'o-','Disp','FDM'); xtight; sms(2); hold on
  label('T','relative S (impurity contribution)');

  ytight(1.1); set(gca,'YTick',1:10); grid on

  legdisp('Location','SouthEast');

end

