function cmap0w(varargin)
% function cmap0w([opts])
%
%    set colormap ensuring that z=0 is shown as white
%    positive values in blue, and negative values in red
%    using symmetric colormap for neg and pos values.
%
% Options:
%
%     'dim',..   length of colormap (256)
%     'cfac',..  factor on clim (1) where cfac<1 increases the contrast
%     'alpha',.. linear scaling for alpha=1 (otherwise: [0..1] .^ alpha)
%     'ah',..    axis handle(s) to look at for surface handles
%     'sh',..    explicitely specified surface handle
%     'all'      set all handles to the same clim value (scale all the same)
%
% Wb, Oct12,02

   getopt ('init',varargin);
      ndim = getopt('dim', 256);
      alpha= getopt('alpha',1.);
      cfac = getopt('cfac',1.);
      sall = getopt('all'    );
      ah   = getopt('ah',[]);
      sh   = getopt('sh',[]);

      cp   = getopt('cp',[0 0 1]);
      c0   = getopt('c0',[1 1 1]);
      cn   = getopt('cn',[1 0 0]);

   a=getopt('get_last',[]);
   od={}; 

   if ~isempty(a)
      if isempty(ah) && all(isaxis(a)), ah=a;
      else od={'-data',a}; end
   end

   if cfac<=0
      eval(['help ' mfilename]); return
      wbdie('invalid cfac = %g', cfac);
   end

   if ~isempty(sh), ah=get(sh,'parent'); end

   if isempty(ah)
      if sall
           ah=findall(gcf,'Type','axes')';
      else ah=gca; end
   end

   if isempty(sh), sh = [
      findall(ah,'Type','surface')
      findall(ah,'Type','patch')
   ]'; end

   ndim = ndim - rem(ndim,2);

   dc = 1/(ndim/2);
   xx = (0:dc:1) .^ alpha;

   cm = [ ...
      ones(ndim/2+1,1) * c0 + flipud(xx'       ) * (cn-c0); ...
      ones(ndim/2  ,1) * c0 +       (xx(2:end)') * (cp-c0)];
   colormap(ah,cm);

   if sall
       if isempty(sh)
          disp('ERR - no suface handles found.')
          eval(['help ' mfilename]); return
       end

       [a,b,cl]=getDataRange_(od{:},sh);

       if ~isempty(cl)
       set(h,'clim',cl*cfac); end
   else
       for h=ah
           sh = [ findall(h,'Type','surface'); findall(h,'Type','patch')]';

           [a,b,cl]=getDataRange_(od{:},sh);

           if ~isempty(cl)
           set(h,'clim',cl*cfac); end
       end
   end

end

% -------------------------------------------------------------------- %
function [a,b,cl]=getDataRange_(varargin)

 if nargin>1 && isequal(varargin{1},'-data')
    a=0; b=0;
    for i=2:length(varargin);
        dd = varargin{i}(:);
        a = min ([a, min(dd)]);
        b = max ([b, max(dd)]);
    end
 else
    a=0; b=0; sh=varargin{1};
    for h2=sh
        dd = get(h2,'CData');
        a = min ([a, min(min(dd))]);
        b = max ([b, max(max(dd))]);
    end
 end

 ab = max(abs(double([a b])));

 if ab~=0
       cl=[-ab, +ab];
  else cl=[]; end

end

% -------------------------------------------------------------------- %

