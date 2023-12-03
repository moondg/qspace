function varargout=getEDdata(varargin)
% function [[EK,]ED,Inrg]=getEDdata([NRG ,opts])
%
%    get energy data of discarded state space from specified NRG run
%    (default: 'NRG/NRG'). All energies are converted to physical
%    energies w.r.t. the ground state energy of the last iteraion,
%    i.e. all NRG scaling and shifting are undone.
%
%    The returned cell array ED{k} contains one matrix [Ed,deg]
%    for every Wilson iteration k.
%
% To be used for calculating specific heat, entropy, and similar.
%
% Options
%
%    -K   also return combined data for each iteration (i.e. K+D) as first argument
%    -Q   include symmetry labels
%    -m   merge ED data into a single 3-colmun matrix [k,Ed,deg]
%         of Wilson iteration k, discarded energies Ed and
%         their respective degeneracies
%    -d   get also local dimensions along Wilson chain from AK
%         (returned as Inrg.d)
%    -q   quiet mode.
%
% Wb,Jul17,11

% adapted from specificHeatNRG.m

  getopt('INIT',varargin);
     merge=getopt('-m');
     Kflag=getopt('-K');
     dflag=getopt('-d');
     qflag=getopt('-q');
     Qflag=getopt('-Q');
  NRG=getopt('get_last','NRG/NRG');

  if Kflag && nargout>3 || ~Kflag && nargout>2
     wbdie('invalid usage (too many output arguments requested)'); 
  end

  if ischar(NRG)
     f='_info.mat'; m=[NRG f];
     if ~exist(m,'file'), NRG, wbdie('invalid NRG data'); end
     if ~qflag
        nstr=sprintf('   %s: loading %s',mfilename,repHome(NRG));
        fprintf(1,'\r%s%s ...  ',nstr,f);
     end
     Inrg=load(m);
  else
     if ~iscell(NRG) || numel(NRG)~=2 || ...
        ~isfield(NRG{1},'AK') || ~isfield(NRG{2},'EScale')
        wberr('invalid NRG data specification');
     end
     Inrg=NRG{2};
  end

  ES=Inrg.EScale;
  E0=Inrg.E0.*ES;

% dE=cumsum(E0); Eref=min(dE); dE=dE-Eref;
% NB! this looses energy resolution towards the end of the system
% NB! with Eref=sum(E0), above cumsum() is equivalent to performing
% (-1)*cumsum() wrt. to the end of the system! => numerically stable.
% Wb,Jul19,13

  dE=-fliplr(cumsum(fliplr(E0)));
  dE=[ dE(2:end), 0 ];

  if min(dE)<0, wblog('WRN','got min(dE)=%.3g !?',min(dE)); end
  if abs(Inrg.phE0-sum(E0))/Inrg.phE0>1E-8
     wbdie('ground state energy discrepancy (%.3g)',sum(E0)-Inrg.phE0);
  end

  if dflag, vars={'AK','HD'}; else vars={'HD'}; end
  if Kflag, vars{end+1}='HK'; end

  HD=QSpace; k=0; 
  while 1
     if ischar(NRG)
        f=sprintf('_%02d.mat',k); m=[NRG f];
        if ~exist(m,'file'), break; end
        if ~qflag, fprintf(1,'\r%s%s ...  ',nstr,f); end
        load(m,vars{:}); k=k+1;
     else
      % bug-fix: Matan Lotem [11/2023]
        if k<numel(NRG{1}), k=k+1; else break; end
        structexp(NRG{1}(k));
     end

     gotHD=~isempty(HD.data);

     if k==1
        if gotHD, wbdie('got truncation at 0th NRG iteration !?'); end
        if Kflag && isdiag(QSpace(HK),'-d')<2
           [~,q]=eigQS(HK); AK=q.AK; HK=q.EK;
        end
     elseif k==2
        if Kflag && isdiag(QSpace(HK),'-d')<2
           wbdie('got non-diagonal HK at 1st NRG iteration');
        end
     end

     if gotHD, HD_=HD.data;
      % bug-fix: loaded HD by now QSpace // Matan Lotem [11/2023]
        ED(k).E = ES(k)*[HD_{:}]' + dE(k);          % physical energies
        ED(k).deg=getzdim(QSpace(HD),2,'-p','-x');  % degeneracies

        if Qflag
           nd=numel(HD.data); qq=cell(1,nd); Q=HD.Q{1};
           for i=1:nd
              qq{i}=repmat(Q(i,:),length(HD.data{i}),1);
           end
           ED(k).Q=cat(1,qq{:});
        end
     else HD_={};
     end

   % consider only iterations where truncation occured
   % eg. allow impurity to have altered state space
     if dflag && gotHD && ~isempty(AK) && ~isempty(AK.data)
        q=getDimQS(AK); d(k)=q(end);
     end

     if Kflag
        if ~isempty(HK.data), HK_=HK.data; else HK_={}; end

      % bug-fix: loaded HK by now QSpace // Matan Lotem [11/2023]
        EK(k).E = ES(k)*[HK_{:}, HD_{:}]' + dE(k);
        EK(k).deg=[
            getzdim(QSpace(HK),2,'-p','-x')
            getzdim(QSpace(HD),2,'-p','-x') ];
        if Qflag
           nd=numel(HK.data); qq=cell(1,nd); if nd, Q=HK.Q{1}; end
           for i=1:nd
              qq{i}=repmat(Q(i,:),length(HK.data{i}),1);
           end
           qq=cat(1,qq{:}); if gotHD
           EK(k).Q=[qq; ED(k).Q ]; else EK(k).Q=qq; end
        end

        [EK(k).E,is]=sort(EK(k).E);
        EK(k).deg=EK(k).deg(is); if Qflag
        EK(k).Q  =EK(k).Q(is,:); end
     end
  end
  if ~qflag, fprintf(1,'\r%70s\r',''); end

  if dflag, 
     Inrg.d=d;
     d=unique(d(find(d))); if numel(d)~=1, d
       wbdie('failed to determine local d'); end
     Inrg.dloc=d;
  end

  if k<6, wblog('WRN','Wilson chain appears to short'); end
  N=k;

  if merge
     for k=1:N
        if isempty(ED(k).E), continue; end
        q={ repmat(k-1,numel(ED(k).E),1), ED(k).E, ED(k).deg };
        if Qflag, q{end+1}=ED(k).Q; end
        ED(k).dd=cat(2,q{:});
     end
     ED=cat(1,ED.dd);
  end

  if ~nargout, varargout={ED};
  elseif Kflag
       varargout={EK,ED,Inrg};
  else varargout={   ED,Inrg};
  end

end

