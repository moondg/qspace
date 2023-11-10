function [wd,TT]=getWeightsNRG(varargin)
% function [wd,TT]=getWeightsNRG([NRG,]TT)
%
%    Get FDM weight distribution for given NRG run
%    (default: $LMA/NRG/NRG.mat).
%
% Wb,Jun20,13

  getopt('init',varargin);
     vflag=getopt('-v');
     TT=getopt('T',-1);
  args=getopt('get_remaining'); nargs=length(args);

  if nargs
     if ischar(args{1})
        NRG=args{1}; args(1)=[];
     elseif iscell(args{1}) && numel(args{1})==2 && ...
        isfield(args{1}{1},'AK') && isfield(args{1}{2},'EScale')
        NRG=args{1}; Inrg=NRG{2}; args(1)=[];
     else NRG=[getenv('LMA') '/NRG/NRG'];
     end
  end

  if numel(args)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if ischar(NRG)
     f=sprintf('%s_info.mat',NRG);
     if ~exist(f,'file')
        wbdie('invalid NRG data (inaccessible info file) !??');
     end
     Inrg=load(f);
  end

  if vflag, o={}; else o={'-q'}; end
  [ED,IN]=getEDdata(NRG,'-d',o{:}); N=numel(ED); nT=numel(TT);

  dloc=IN.dloc;
  if dloc<2 || dloc~=round(dloc)
     wbdie('invalid d_loc = %g !??',dloc);
  end

  for k1=1:N
     if ~isempty(ED(k1).E), break; end
  end

  E0=min(cat(1,ED.E)); if E0~=0
     wblog('WRN','got small energy E0=%.3g !??',E0);
     for k=k1:N, ED(k).E = ED(k).E - E0; end
  end

  for it=1:nT, T=TT(it);
     if T<0, if T==-1, T=-10; end
        T=Inrg.Lambda^(-(N+T)/2);
        TT(it)=T;
     end
     rw=zeros(N,1);

     if T==0, rw(end)=1;
     else
        for k=k1:N, q=ED(k); if ~isempty(q.E)
           rw(k) = q.deg' * exp(-q.E/T + log(dloc)*(N-k));
        end, end
        rw=rw/sum(rw);
     end

     wd(:,it)=rw;
  end

end

