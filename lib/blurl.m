function [lh,Ib]=blurl(varargin)
% function blurl([lh, opts])
%
%    blur line by making it brighter with linewidth 2
%
% Options
% 
%    lh       set of line handles (if not provide, then all lines
%             within current axis set are taken)
%   'cfac',.. factor of how much to blur (0=no blur; 1=white) (0.79)
%             (may also be specified as plain number right after lh)
%             -1<=cfac<0 tunes the color to darker with -1 being black.
%
%   'cmin',.. min of color to keep [relative to sqrt(white)=sqrt(3)]
% 
%   --undo    undo prior blurring by restoring earlier formats
%   --line    blur line only (i.e. keep marker as it is)
%   --bg      move lines to background
%
%   Any remaining options are directly applied to line handles (lh).
%
% Wb,Oct31,06

% Wb,Feb24,21 : permitting cfac<0 (make lines darker)

  if ~nargin
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if isempty(varargin{1}), lh=varargin{1}; return; end

  bflag=0;

  if ishandle(varargin{1}), lh=varargin{1}; varargin(1)=[];
  elseif isequal(varargin{1},'-bar'), bflag=1;
       lh=findall(gca,'type','hggroup');
  else lh=findall(gca,'type','line'); end

  if nargin>1 && isnumeric(varargin{1})
       cfac=varargin{1}; varargin(1)=[];
  else cfac=[]; end

  getopt('init',varargin);
     if isempty(cfac),       cfac=getopt('cfac',cfac);
        if isempty(cfac),    cfac=getopt('fac', cfac);
           if isempty(cfac), cfac=0.79; end
     end, end

     cmin =getopt('cmin',0.);
     toback=getopt('--bg');

     undo_blur=getopt('--undo');
     lflag=getopt('--line' );

  opts=getopt('get_remaining'); nopts=numel(opts);

  if bflag
     for i=1:length(lh), ch=get(lh(i),'Children');
        c=get(lh(i),'Color'); set(ch,'Color',c,'ZData',[],'Visible','on');
        if undo_blur, continue; end

        if cfac<=1 && cfac>=-1
           if cfac>=0
                c = (1-cfac) * (c-1) + 1;
           else c = (1+cfac) *  c;    end

           for j=2
              x=get(ch(j),'XData');
              z=repmat(-0.01,size(x)); z(find(isnan(x)))=nan;
              set(ch(j),'ZData',z,'Color',c);
           end
        else
           set(ch(2),'Visible','off');
        end
     end

   % errorbars will still overlap with data points even though put to
   % background by setting z<0 above.
   % however, replotting the data without error bars still will find
   % itself partly covered with errorbars !-93_()*@!&$_
   % Wb,May07,09

   % if ~undo_blur, hold on
   %    for i=1:length(lh), lo=lineopts(lh(i));
   %       x=get(lh(i),'XData');
   %       y=get(lh(i),'YData'); z=repmat(0.1,size(x));
   %       plot3(x,y,z,lo{:},'tag','blurl_bar');
   %    end
   % end

     if ~nargout, clear lh; end
     return
  end

  if undo_blur
     for h=reshape(lh,1,[])
        u=get(h,'UserData');
        if isstruct(u) && isfield(u,'w') && isfield(u,'col')
        set(h,'LineWidth', u.w, 'Color', u.col); end
     end
  else
     cmin=sqrt(3)-cmin;

     lw=[]; if nopts
        i=findstrc(opts,'LineW','-i','-x');
        if ~isempty(i), lw=opts{i(1)+1}; end
     end
     if isempty(lw)
        if lflag, lw=1;
        else lw=max(1,min(2, 1 + (cfac-0.5)/(0.79-0.5)));
        end
     end

     for h=reshape(lh,1,[]), c=get(h,'Color');
        if lflag
           set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c);
        end
        u=struct('w',get(h,'LineWidth'),'col',c);
        set(h,'UserData',u);

        if cfac>=0
            if cfac<= 1, c2 = (1-cfac) * (c-1) + 1; else c2=[1 1 1]; end
        elseif cfac>=-1, c2 = (1+cfac) *  c;        else c2=[0 0 0]; end

        if cfac<=0 || norm(c2)<cmin, c=c2; end

        set(h,'Color',c,'LineW',lw);
        if ~lflag && ~isequal(h,'Marker','none')
        set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c); end
     end
     if nargout>1
        Ib=add2struct('-',cfac,'LineW=lw',cmin,toback,'lines_only=lflag');
     end
  end

  if nopts, set(lh,opts{:}); end

  if toback, mv2back(lh); end
  if ~nargout, clear lh; end

end

