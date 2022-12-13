function [HAM,IS]=Hamilton1D(varargin)
% function [HAM]=Hamilton1D(type, L [,opts])
%
%    setup class for one-dimensional DMRG Hamiltonians
%    in a traditional setting in a quasi-MPO structure.
%    The class already also applies for the non-abelian
%    setting.
%
% Wb,Apr06,14 ; Wb,Aug25,15

% Examples:
% [HAM]=Hamilton1D('Heisenberg',{12,'sym','SU4','qloc',[0 2 0],'-perBC'},'-mat'); 
% [HAM]=Hamilton1D('tightbinding',{12,'-perBC'},'-mat'); 

  if ~nargin
   % e.g. called at the very beginning of each private/setup_*.m routines
     HAM=class(setup_empty(),'Hamilton1D'); IS=[];
     return
  end

  if nargin==1, H=varargin{1};
     if isa(H,'Hamilton1D'), HAM=H; return
     elseif isstruct(H)
        H_=setup_empty(); [H,e,s]=addfields(H,H_);
        if e<0, wblog('WRN','%s() altered input structure:\n%s',mfilename,s);
        elseif e>0, wberr('invalid input structure:\n%s',s); end

        if isfield(H.mpo,'idx') && ~isfield(H.mpo,'iop')
           for i=1:numel(H.mpo), H.mpo(i).iop=H.mpo(i).idx; end
           H.mpo=rmfield(H.mpo,'idx');
        end
        if ~isfield(H.mpo,'user'), H.mpo(1).user=[]; end

        HAM=class(H,'Hamilton1D');
        return
     end
     helpthis; wberr('invalid usage');
  end

  if nargin<2 || ~ischar(varargin{1}), helpthis
     if nargin || nargout, wberr('invalid usage'), end, return
  end

  getopt('init',varargin(3:end));
     tflag=getopt('-t');
     mflag=getopt('-mat'); store=''; fout=''; if ~mflag
     store=getopt('store', store);  if isempty(store)
     fout=getopt('fout',[]); end; end
  getopt('check_error');

% -------------------------------------------------------------------- %
  switch varargin{1}
    case 'Heisenberg',   [HAM]=setup_Heisenberg(varargin{2}{:});
    case 'AKLT',         [HAM]=setup_AKLT(varargin{2}{:});
    case 'HsbgLadder',   [HAM]=setup_HeisenbergLadder(varargin{2}{:});
    case 'HBTriLadder',  [HAM]=setup_HeisenbergTriLadder(varargin{2}{:});
    case 'HBLadderX',    [HAM]=setup_HeisenbergLadderX(varargin{2}{:});
    case 'HBLadderP2',   [HAM]=setup_HeisenbergLadderP2(varargin{2}{:});
    case 'Hubbard',      [HAM]=setup_Hubbard(varargin{2}{:});
    case 'TLHubbard',    [HAM]=setup_TLHubbard(varargin{2}{:});
    case 'KondoNecklace',[HAM]=setup_KondoNecklace(varargin{2}{:});
    case 'HBL_BaIrO',    [HAM]=setup_Heisenberg_BaIrO(varargin{2}{:});
    case 'SLH_chiral',   [HAM]=setup_HeisenbergSL_chiral(varargin{2}{:});
    case 'Kitaev',       [HAM]=setup_Kitaev(varargin{2}{:});
    case 'tightbinding', [HAM]=setup_tightbinding(varargin{2}{:});
    case 'tb_ladder',    [HAM]=setup_tb_ladder(varargin{2}{:});
    case 'tb_testferm',  [HAM]=setup_tb_testferm(varargin{2}{:});
    case 'empty',        [HAM]=setup_empty(varargin{2}{:});
    otherwise wberr('invalid system ''%s''',varargin{1});
  end
  if isfield(HAM.info,'XY')
     [Q,I,D]=uniquerows(round(100*HAM.info.XY)); i=find(D>1);
     if ~isempty(i)
        HAM.info.XY([I{i}],:)
        wberr('got overlapping XY data !?');
     end
  end

  if ~ischar(store), wberr('invalid store specification'); end

  HAM.mat=[];

  if mflag
     if ~isfield(HAM,'store') || isempty(HAM.store) || ...
        ~ischar(HAM.store), wberr(...
        'model setup must specify default HAM.store (global)');
     end

     HAM.mat=[getenv('LMA') '/DMRG/' HAM.store];
     HAM.store=[];
  elseif ~isempty(store)
     if ~isempty(findstr(store,'/')) 
          HAM.mat=store;
     else HAM.store=store; end
  elseif ~isempty(fout)
     HAM.mat=fout;
     HAM.store='';
  else
     if ~xor(isempty(HAM.store),isempty(HAM.mat))
     wberr('invalid output setting'); end
  end

  if nargout>1
     IS=getfields(HAM.info.IS,'istr','sym','E','SOP');
     IS.SOP=getfields(IS.SOP,'info','type');
     IS.RC_STORE=getenv('RC_STORE');
  end

  if tflag, o={'-t'}; else o={}; end
  e=consistency_check(HAM,o{:});
  HAM=class(HAM,'Hamilton1D');

  if tflag
       wblog('TST','skipping storage setup'); 
  else setupStorage(HAM); end

end

% -------------------------------------------------------------------- %
