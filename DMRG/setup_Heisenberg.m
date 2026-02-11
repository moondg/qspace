function [HAM,Ix]=setup_Heisenberg(varargin)
% function [HAM,Ix]=setup_Heisenberg([L,][opts])
%
%    Setup of plain spin-S Heisenberg chain (default: S=1/2) 
%    with (anisotropic) nearest-neighbor couplings [J [,Jz]].
%    This setup uses either SU(2) [Jz=J] or U1 [Jz~=J] symmetry.
%
%    For more elaborate Hamiltonian parameter differentiation
%    see also @Hamilton1D/setup_Heisenberg
%    e.g. as called with runDMRG with wsys='HAM::Heisenberg'.
%
% Wb,Jan16,25

  if nargin && isnumber(varargin{1})
       L=varargin{1}; varargin(1)=[];
  else L=[]; end

  getopt('init',varargin);
     J    = getopt('J',   1);
     S    = getopt('S', 1/2);

     L    = getopt('L', L); if isempty(L), L=32; end
     perBC= getopt('-perBC');

     Nspeed=getopt('Nspeed',[]); % see initNKEEP() for usage
     nk1  =getopt('nk1',3);      % -> Nkeep ∈ [2^nk1, 2^nk2]
     nk2  =getopt('nk2', max(nk1+1,floor(log(L)/log(2))));
	   %  L-range   nk2
	   %  [ 16, 32[    3
	   %  [ 32, 64[    4
	   %  [ 64,128[    5
	   %  [128,256[    6

     use_mem = getopt('--mem');
     odir = getopt('odir','');
     ftag = getopt('ftag','');

     if     getopt('-v'), vflag=2; o={'-V'};
     elseif getopt('-q'), vflag=0; o={};
     else                 vflag=1; o={'-v'}; end

     tflag= getopt('-t');
     kflag= getopt('-k');

     use_mpo = getopt('--mpo',{'not_specified'});

  getopt('check_error');

  HAM=struct(Hamilton1D); % calls setup_empty()

  param=add2struct(L,S,J,perBC);

  if isempty(L), wbdie('length L not specified'); end

  n=numel(J); if n==1, dJz=0;
  elseif n==2, dJz=diff(J); if ~dJz && J(1)>0, J=J(1); end
  else J, wbdie('invalid usage'); end

  if norm(J)<1e-3, J
   % try to keep energies of order 1.
     wbdie('invalid usage (got tiny couplings)'); 
  end

  q=2*S; if q<0 || q>10 || mod(q,1), S, wbdie('invalid spin S'); end
  istr={ ['spin-' wbrat(S)], '' };

  U1sym=0; if vflag, o={'-v'}; else o={}; end
  if dJz || J(1)<0
     U1sym=1; o{end+1}='-A'; % Abelian U(1)
     if numel(J)==1, J=[J J]; end
  end

  if U1sym
       istr{2}=sprintf('J=[%g %g]',J);
  else istr{2}=sprintf('J=%g',J); end
  istr=['spin-' istr{1} ' Heisenberg chain with ',istr{2}];

% central definition of local state space of site
% described then by the spin operator(s) S
  [S,IS]=getLocalSpace('Spin',S,o{:});

  initNKEEP
  Ix=add2struct('-',NKEEP,Nspeed,Nkeep,param);

  HAM.info.istr=istr;
  HAM.info.IS=IS;
  HAM.info.param=param;

  HAM.oez={{IS.E,'local identity operator (E)'}}; % init_ops

  if U1sym
       Y=[sum(S(2:3)), S(1)]; % [ +/-, z ]
  else Y=S; end

  nJ=[numel(J), numel(Y)];
  if diff(nJ), wbdie('unexpected setup (nJ=%d/%d)',n);
  else nJ=nJ(1); end

  j2=sqrt(abs(J));
  is=find(J<0);

  Q=QSpace; for i=1:nJ, Q=Q+j2(1,i)*Y(i); end
  HAM.ops={{Q,'Heisenberg-coupling','~hconj'}};       % init_ops

  if ~isempty(is), j2(is)=-j2(is);
     Q=QSpace; for i=1:nJ, Q=Q+j2(i)*Y(i); end
     HAM.ops{1}{2} = [ HAM.ops{1}{2} '-1' ];
     HAM.ops{2,1}={Q,'Heisenberg-coupling-2','~hconj'}; % init_ops
  end

% coordinates for plot(HAM)
  XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
  HAM.info.XY=XY;

  HH=zeros(L-1,5); l=1;
  i2=numel(HAM.ops); stype=[];

  if perBC
     wblog('NB!','using periodic BC (interleaved setup)');
     for i=1:L
        HH(l,:)=[ [i-1, 1  ], [i+1, i2 ], 1]; l=l+1;
     end
     HH(1)=1; HH(end,3)=L;
  else
     wblog(' * ','using plain open BC');
     for i=1:L-1
        HH(l,:)=[ [i, 1  ], [i+1, i2 ], 1]; l=l+1;
     end
  end

  if isempty(ftag), ftag='Heisenberg'; end
  if ~isempty(odir)
       fout = [odir '/' ftag];
  else fout = [ 'DMRG/' ftag]; end

  oH=setopts('-',fout,{'--mpo',{0}},HH,'stype?');

  [HAM]=Hamilton1D(HAM,oH{:});

  if kflag, keyboard, end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

