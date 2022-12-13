function mpdf(varargin)
% Function mpdf([fname,OPTS])
% create pdf file from current figure (default fname='figure.pdf')
%
% Options
%
%   '-f'           force overwrite if pdf already exists
%   '-q','-nodisp' do not show pdf in default previewer
%   '-v'           verbose
%
% Wb,Jan26,05

  getopt ('init', varargin);
     nodisp = getopt ({'-q','-nodisp'});
     force  = getopt ('-f');
     verbose= getopt ('-v');
  varargin=getopt('get_remaining');

  if numel(varargin)
     fname=varargin{1};
     if isempty(findstr(fname,'.pdf')), fname = [ fname '.pdf' ]; end
  else
     fname=getfigname(gcf); % get(gcf,'Name');
     fname=regexprep(fname,'.*\/\/','');
  end

% ----------------------------------------------------------------------
  f=gcf;

  if verbose
  fprintf(1,['\n   Saving current figure (%g,%s) to \n' ...
    '     pwd : %s\n' ...
    '     file: %s\n\n'], f, get(gcf,'Renderer'), pwd, fname);
  end

  if isempty(fileparts(fname)), fname=['./' fname]; end

  if exist(fname) && ~force
     s='exists.\n   Overwrite? <[1]|0> ';
     if verbose
          inp = input (['   File ' s]);
     else inp = input ([sprintf('\n   File `%s'' ',fname), s]); inl(1)
     end
     if isempty(inp); inp=1; end
     if inp~=1; disp ' '; return; end
  end

  CMD = sprintf('saveas(gcf, ''%s'', ''pdf'');', fname);
  if verbose, fprintf('   %s\n\n',CMD); end

  if isequal(get(0,'ScreenSize'),[1 1 1 1])
     set(gcf,'PaperPosition',[0.2 0.2 7.86 11.29 ])
  else
     set(f,'PaperPositionMode','auto')
     set(f,'PaperType',get(0,'DefaultFigurePaperType'))
  end

  eval(CMD);

  if ~nodisp
     p='kpdf';
     eval(['! ' p ' ' fname ' &']);
  end

end

