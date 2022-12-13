function ah=nrg_header(varargin)
% function ah=nrg_header([opts])
%
%    Display content of global param structure as left footer
%    If param.istr is present, it is shown as left header.
%    This also calls the routine addfinfi (add file info)
%    which lookad at the `mat' string variable if present,
%    as well as Imain (main info structure generated in jobs)
%    The extracted mat/host file info is shown as right footer.
%
%    Global variables consulted: param savemat
%    Variables read from caller workspace if they exist:
%    Inrg Imain
%
% Options
%
%   '-x'         regex exclude pattern of fields in param to ignore
%   '--pval',..  replace (or set) values in local copy of global param
%                prior to displaying it (* x 2 cell array)
%   '--regex',.. list of regexprep pattern to apply to param string
%                prior to showing it in left footer (* x 2 cell array)
%
% Wb,Jul04,07

% -------------------------------------------------------------------- %
% Examples
%
%    nrg_header('-x','\<(z|ALambda)\>');
%
%    nrg_header('-x','gg|FL',...
%        '--regex',{'ximp','x_{imp}'; 'TK1','T_K^{(1)}'}, ...
%        '--pval',{'J',J});
%
% -------------------------------------------------------------------- %

  global param savemat

% setglobal was reported as missing in MCC // Wb,Mar31,22
% -> repeat to make sure setglobal is actually called
  setglobal param

  evalin('caller','setglobal savemat');

  getopt('init',varargin);
     dy =getopt('dy',0);
     xp_=getopt('-x','');
     rx_= getopt('--regex',{});
     px_= getopt('--pval', {});
  Inrg=getopt('get_last',[]);

  try   checkParam(0);
  catch l
     wblog('WRN','checkParam: %s',l.message);
  end

  p=param;
  if ~isempty(px_)
     if ~iscell(px_) || size(px_,2)~=2
        wbdie('invalid ''--pval'' (..x 2 cell array required)'); end
     for i=1:size(px_,1)
        if ~ischar(px_{i,1})
           wbdie(['invalid ''--pval'' (first column in cell array ' ... 
             'must be field names, i.e., strings)']); end
        p=setfield(p,px_{i,1},px_{i,2});
     end
  end
  structexp(p)

  if 0 && isfield(p,'Lambda') && ...
     sum(isfield(p,{'J','JH','Gamma'}))==1 && ~isfield(p,'U3')

   % -------------------------------------------------------------- %
     global gTK; nTK='T_K';
     if ~isempty(gTK)
        if iscell(gTK)
           TK=gTK{2}; nTK=strrep(gTK{1},'\','\\');
        else TK=gTK; end
     else TK=TKondo; end

     sTK=sprintf([nTK '=%s'], num2tex(TK,'%.3E'));

     if ~exist('D','var')
        if exist('hbdim','var'), D=hbdim; end
     end

     if ~exist('L','var')
        if exist('Psi','var'), L=(length(Psi)-2)/2; end
        if exist('N','var'), L=N; end
     end

     if exist('Psi','var'), ss='DMRG'; else ss='NRG'; end

     if ~exist('D','var') && exist('Nkeep','var'), D=Nkeep; end

     lam={sprintf('\\Lambda=%g',Lambda)};
     if isfield(p,'z') && p.z, lam{end+1}=sprintf('z=%.3g',z);
     elseif isfield(p,'nz') && p.nz, lam{end+1}=sprintf('n_z=%g',nz);
     end

     if isfield(p,'ALambda') && ~isempty(p.ALambda)
        if ischar(p.ALambda), lam{end+1}='A_\Lambda';
        else lam{end+1}=sprintf('A_\Lambda=%g',ALambda); end
     end
     if numel(lam)>1, lam=[ lam{1},' (',strhcat(lam(2:end),'-s',', '),')'];
     else lam=lam{1}; end

     ss = {
        sprintf('%s  N=%g, %s, D=%d',ss,L,lam,D);
        ''
        sTK
     };

     if exist('U','var') && exist('Gamma','var')
        s={sprintf('\nSIAM: U=%g', U)};
        if exist('epsd','var')
        s{end+1}=sprintf(', \\epsilon_d=%s', vec2str(epsd, 'sep',', ')); end
        s{end+1}=sprintf(', \\Gamma=%s', mat2str2(Gamma,'fmt','%.5g','rowsep','; '));
        if exist('J','var'), s{end+1}=sprintf(', J=%g', J); end
        ss{2}=cat(2,s{:});
     elseif exist('J','var') || exist('fJ','var')
        if ~exist('J','var'), J=fJ;
        elseif exist('fJ','var') && J~=fJ
        wblog('ERR','J inconsistency (%g,%g)',J,fJ); end

        ss{2}=sprintf('; Kondo: J=%g',J);
     end

     if exist('Bfield','var') && ~exist('B','var'), B=Bfield; end
     if exist('B','var')
        if TK~=0 && ~isnan(TK)
             ss{2}=[ss{2}, sprintf(', B=%.3g TK', B/TK)];
        else ss{2}=[ss{2}, sprintf(', B=%.3g', B)]; end
     end
     ss=sprintf('%s %s (%s)', ss{:});

  elseif isfield(p,'Lambda') && any(isfield(p,{'J','JH','Gamma'}))

     xpat='nrgIO|wsys|ff|mat|sym|istr'; if ~isempty(xp_), xpat=[xpat,'|',xp_]; end
     ss=param2str(p,'--tex','-x',xpat);

     if isfield(p,'T'), ss=regexprep(ss,'T=[\deE\.+-]+',...
        ['T=',regexprep(num2tex(T,'%.2E'),'\\','\\\\')]); end

  elseif isfield(p,'tLR')
     xpat='IO|mat|sym'; if ~isempty(xp_), xpat=[xpat,'|',xp_]; end
     ss=param2str(p,'-x',xpat);

     ss=regexprep(ss,'([^\\])Lambda','$1\\Lambda');
     ss=regexprep(ss,'([^\\])Gamma','$1\\Gamma');
     ss=regexprep(ss,'([^\\])delta','$1\\delta');
     ss=strrep(ss,'Vbias','V_{b}');
     ss=strrep(ss,'tLR','t_{LR}');
     ss=strrep(ss,'tdot','t_{dot}');
  elseif ~isempty(p)
     xpat='IO|mat|sym|istr|POP'; if ~isempty(xp_), xpat=[xpat,'|',xp_]; end
     ss=param2str(p,'-tex','-x',xpat);
  else ss=''; end

  rx = {
     '(ALambda|A\\Lambda).*energies''', 'A_{\\Lambda}'
     'Nkeep',       'N_K'
     'Etrunc',      'E_{tr}'
     '\<([Nn])z\>', '$1_z'
     'sym=', ''
  };
  if ~isempty(rx_)
     if ~iscell(rx_) || size(rx_,2)~=2
        wbdie('invalid ''--regex'' (..x 2 cell array required)'); end
     rx=[rx; rx_];
  end

  for i=1:size(rx,1)
     ss=regexprep(ss,rx{i,1},rx{i,2});
  end

  ss=regexprep(ss,'([a-zA-Z])*=''[a-zA-Z]*(=|\\leq|\\geq)([^'']*)''','$1$2$3');

  getbase Imain
  if isempty(Inrg), getbase Inrg; end

  if isempty(Imain) && isempty(Inrg)
     q0='---xvar-tmp---';

     setuser(0,'xvar',q0);
     evalin('caller','if exist(''Imain'',''var''), setuser(0,''xvar'',Imain); end');
     q=getuser(0,'xvar'); if ~isequal(q,q0)
        Imain=q;
        setuser(0,'xvar',q0);
     end

     evalin('caller','if exist(''Inrg'',''var''), setuser(0,''xvar'',Inrg); end');
     q=getuser(0,'xvar','-rm'); if ~isequal(q,q0)
        Inrg=q;
     end

  end

  if isfield(Inrg,'NK') && isfield(p,'Etrunc')
     q=max(Inrg.NK,[],1); q={ int2str2(q(1)), int2str2(q(2)) };
     ss=regexprep(ss,'[, ]*D=[^,]*','');
     ss=regexprep(ss,'([, ]*)N(k|_K)[^,]*',sprintf('$1N_K\\\\leq%s (%s)',q{:}));
  end

  if isfield(Imain,'finished')
     if ~isempty(findstr(Imain.finished,'not finished'))
        if exist('Idma','var') && isfield(Idma,'finished')
           Imain.finished=Idma.finished;
        elseif exist('Inrg','var') && isfield(Inrg,'finished')
           Imain.finished=Inrg.finished;
        end
     end
  end

  addt2fig wb

  fh=addfinfo;

  if isfield(p,'istr') && ~isempty(p.istr)
     s=dbstack; if numel(s)>1
          s1=[ strrep(s(2).name,'_','\_') ': '];
     else s1=''; end
     s=p.istr;
     s=regexprep(s,'one-channel *','');
     [hh,ah]=header([s1 s]); set(hh,'FontSize',10);
  end

  [hh,ah]=header('fleft',ss);
  set(hh,'FontSize',8, 'VerticalAlign','bottom');

  if dy
    if ~isempty(fh) && isempty(get(hh,'UserData'))
       set(hh,'Pos', get(hh,'Pos')+[0 dy 0], 'UserData',1);
    end
  end

  if isfield(p,'sym') && ~isempty(p.sym)
     s=p.sym;
     s=regexprep(s,'SU2(charge|spin)','SU2_{$1}');
     s=regexprep(s,'A(charge|spin)','A_{$1}');

     h=header('hleft');
     set(h,'String',[ get(h,'String') ' using ' s]);
  end

  if ~nargout, clear hh ah; end

end

