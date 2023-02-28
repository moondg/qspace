function mfig(varargin)
% mfig(fname [,opts])
%
%     create pdf-file with given name from current figure.
%     (default name if fname is not specified: figure.pdf)
%
% Options
%
%    'ppi',..  pixels per inch (90; larger value = smaller fig)
%              note that the actual figure size determines output figure size.
%    '-eps'    make .eps instead of .pdf
%    '-ps'     make .ps instead of .pdf
%    '-epg'    same as eps but in grayscale
%    '-png'    make .png instead of .pdf
%    '-all'    make .pdf and .eps
%    '-q'      quiet mode
%    '-f'      force overwriting of preexisting pdf-file with same name.
%    '~h'      hide header
%
% Wb (C) 2005 / Jan26,05 / Feb06,07 / Wb,Aug02,16

  ppi=90;
  if ismac, ppi=70; end

  wfig={'pdf'};
  getopt ('INIT',varargin);
     ppi=getopt('ppi',ppi);
     force=getopt('-f');
     tflag=getopt('-t');

     for s={'-eps','-ps','-epg','-all','-png'}, s=s{1};
        if getopt(s), s=s(2:end);
           if isequal(s,'all'), wfig={'eps','pdf'}; end
           wfig={s}; break
        end
     end

     if getopt('-Q'); vflag=0;
     elseif getopt({'-q','nodisp'}), vflag=1;
     elseif getopt('-v'); vflag=3;
     elseif getopt('-V'); vflag=4;
     else vflag=2;
     end

     hflag=getopt('~h');

  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg && numel(varargin{1})==1 && isfig(varargin{1})
       fh=varargin{1}; varargin(1)=[]; narg=numel(varargin);
  else fh=gcf; end

  if narg
     if narg~=1 || isempty(varargin{1}) || varargin{1}(1)=='-'
        wberr('invalid usage'); end
     fname=varargin{1};
  else
     fname=getfigname(fh);
  end
  oo={ 'PaperOrientation'
       'PaperPosition'
       'PaperPositionMode'
       'PaperUnits'
       'PaperSize'
       'PaperType'
       'Position'
  };

  f0=get(fh,oo);

  fp=f0{end};
  fs=fp(3:4)/ppi;

  if isbatch && isequal(get(0,'ScreenSize'),[1 1 1 1])
     fs=1.5*fs;
  end

  set(fh,'PaperUnits','inch', ...
        'PaperOrientation','portrait',...
        'PaperSize',fs,'PaperPosition',[0 0 fs] ...
   );

  if hflag
     h=findall(groot,'type','axes','tag','frame');;
     if ~isempty(h)
        set(h,'visible','off')
     end
  end

% wblog('TST','ppi=%g', ppi);
% fprintf(1,'\n   %18s: %s\n   %18s: %s\n   %18s: %s\n\n',...
%   'screensize',vec2str(get(0,'ScreenSize')), ...
%   'figpos', vec2str(fp), 'figsize',vec2str(fs));
% for i=1:numel(oo)
%    x={ get(fh,oo{i}), f0{i} };
%    if isnumeric(x{1})
%       x={ sprintf('[%s]',vec2str(x{1})), sprintf('[%s]',vec2str(x{2})) };
%    end
%    fprintf(1,'   %18s: %-26s %s\n',oo{i},x{:});
% end

% -------------------------------------------------------------------- %
  if vflag
     s=cell(1,2);
     if tflag, s{1}='TEST: save'; else s{1}='Saving'; end

     i=find(fname=='/',1);
     if isempty(i), fname=['./' fname]; end

     f=[ regexprep(repHome(fname),'\.(pdf|png|eps)$','') '.' strjoin(wfig,'|') ];
     if fh==gcf, s{2}=[s{2} ' current']; end

     fprintf(1,'\n   %s%s figure (%g,%s) to file\n   %s\n', ...
     s{1:2}, get(fh,'Number'), get(fh,'Renderer'), f);
  end

  for i=1:numel(wfig)
      mfig_1(fh,fname,ppi,wfig{i},force,vflag,tflag);
  end
  if vflag, fprintf(1,'\n'); end

% -------------------------------------------------------------------- %
  set(fh,...
    'PaperOrientation',f0{1}, 'PaperType', f0{6},...
    'PaperUnits',f0{4}, 'PaperPosition', f0{2}   ...
  );

  set(fh,'Position',fp);

  if hflag && ~isempty(h), set(h,'visible','on'); end

end

% -------------------------------------------------------------------- %

function mfig_1(fh,fname,ppi,fext,force,vflag,tflag)

  eps=0;
  switch fext
    case 'eps', ffmt='epsc2';           eps=1;
    case 'epg', ffmt='eps'; fext='eps'; eps=1;
    case  'ps', ffmt='epsc2';           eps=1;
    case 'pdf', ffmt=fext;
    case 'png', ffmt=fext;
    otherwise fext, wberr('invalid fext');
  end

  if isempty(regexpi(fname,['.' fext '$'])), fname = [ fname '.' fext ]; end
  if isempty(fileparts(fname)), fname=['./' fname]; end

  if ~force 
     [i,fname]=fexist(fname,tflag);
     if i, % disp(' huh? ')
     return, end
  end

  set(fh,'FileName',fname);

  if vflag>1 || tflag
     fprintf(1,'   saveas(fh,''%s'',''%s'');\n',fname,ffmt);
     if tflag, return; end
  end

% wblog('TST','saveas(%g,%s,%s) ...',fh,fname,ffmt);
% THIS CAN TERMINATE MATLAB INSTANTANEOUSLY WITHOUT ANY FURTHER MESSAGE
% use direct command, instead [eval: matlab-bug] // Wb,Feb18,11
% CMD=sprintf('saveas(fh, ''%s'', ''%s'');', fname, ffmt);
% eval(CMD);
  saveas(fh,fname,ffmt);

  if vflag>3
     if eps, eval(['! gv --watch ' fname ' &']);
     else
        eval(['! vi.pl ' fname ' &']);
     end
  end

end

% -------------------------------------------------------------------- %

