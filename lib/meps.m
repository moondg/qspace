function meps(fname)
% meps([fname='figure.eps'])
%
%     make (e)ps file from current plot window
%     default output file: figure.eps
%
% Wb (C) 2001  (Apr22,01)

% when called like 'meps fname' then all the input args are strings!

  if ~nargin, fname=getfigname(gcf); end
  if isempty(strfind(fname,'.eps')), fname=[fname '.eps']; end

  cdir = pwd;

% ----------------------------------------------------------------------
% cd /home/weichsel/public_html/ML

  if 0
     d = dir('./');
     ii = find(cat(2,d.isdir)); d = d(ii);

     d2 = [];
     for i=1:length(d)
         if d(i).name(1)~='.'
             if strncmpi(fname, d(i).name, length(d(i).name))
                d2 = d(i).name;
             end
         end
     end
     if ~isempty(d2), eval(['cd ' d2]), end
  end
  disp([ 10 pwd 10])

  if isempty(fileparts(fname)), fname=['./' fname]; end
  if exist(fname)
     inp = input (['File ' fname ' exists. Overwrite? <[1]|0> ']);
     if isempty(inp); inp=1; end
     if inp~=1; disp ' '; return; end
  end

  fprintf(1,['\n   Print current figure (%g,%s) to ' fname ...
  '\n   pwd : %s\n\n'], gcf, get(gcf,'Renderer'), pwd);

  CMD = sprintf( ...
  'print(gcf, ''-depsc2'', ''-tiff'', ''%s'')', fname);

  set (gcf, 'PaperPositionMode','auto')

  h=findall(groot,'type','axes','tag','frame');;
  if ~isempty(h)
     set(h,'visible','off')
  end

  disp([CMD 10])
  eval(CMD);

  eval(['! gv ' fname ' &']);

  eval(['cd ' cdir]);

  if ~isempty(h), set(h,'visible','on'); end

end

