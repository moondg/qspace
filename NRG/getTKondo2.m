function [TK,Itk]=getTKondo2(Simp,varargin)
% function [TK,Itk]=getTKondo2(Simp, [,NRG,TT])
%
%    get (improved) Kondo temperature by using a combination
%    of the dynamical spin-spin susceptibility as well as
%    the mixed correlator <Simp|Stot>.
%
% Options
%
%  'T',..   set of temperatures to calculate spin susceptibility
%           (default: single T towards the end of the Wilson chain)
%  'iS',..  for the construction of Stot, the position of the
%           SU(2) spin symmetry in QSpace.info.qtype is required.
%           if iS is not specified, this routine tries to
%           determine iS from the structure of Simp.
%  '-dyn'   also calculate dynamical susceptibility,
%           and return improved estimate for TK.
%
% Wb,Jul04,13

  getopt('init',varargin);
     vflag=getopt('-v');
     TT=getopt('T',[]);
     iS=getopt('iS',{});
     mt=getopt('mt',4);
     dyn=getopt('-dyn');
  varargin=getopt('get_remaining');

  opts={};
  if vflag, opts{end+1}='-v'; end
  if ~isempty(iS), opts=[opts {'iS',iS}]; end

  if numel(varargin) && ischar(varargin{1})
       NRG=varargin{1}; varargin=varargin(2:end);
  else NRG=[getenv('LMA') '/NRG/NRG']; end

  if numel(varargin) || ~ischar(NRG)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  f=[ NRG, '_info.mat' ]; if ~exist(f,'file')
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  Inrg=load(f);
  Lambda=Inrg.Lambda;

  if isempty(TT)
     TT=(100*Inrg.EScale(end))*Lambda.^((0:prod(mt)-1)/mt(end));
  end

  [ch2,Ich]=getSpinSuscept2(Simp,NRG,'T',TT,opts{:});

  e=std(ch2)/mean(ch2);
  if numel(ch2)>2 && mod(mt,2)==0, l=mt/2;
     q=ch2-mean(ch2); e(end+1)=norm(q(1:l)+q(l+1:end))/mean(ch2);
  end
  Ich.err=e;

  if e(end)>0.05, wblog('WRN',...
    'got significant fluctuation on chi0 (%.3g)',e(end));
  end

  TK=0.25/mean(ch2);

  if dyn
    [om,a0,Ich.Idma] = fdmNRG_QS(...
        NRG,[],Simp,[],'calcRho','nostore','cflags',1);
    TK0=0.25/sum(Ich.Idma.reA0);
    TK2=TK;
    TK=1/(1/TK0-2*(1/TK0-1/TK2));
    add2struct(Ich,TK0,TK2,TK); Ich.dTKrel=(TK-TK2)/TK2;
  end

  if nargout>1
     Itk=add2struct(Ich,ch2);
  end

end

