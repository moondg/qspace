function hm=allmarks(varargin)
% function hm=allmarks([-opts])
%
% Options
%
%   '-fgz'   move all marks to the foreground (by setting zdata)
%
% Wb,Mar21,19

  if nargin && isaxis(varargin{1})
       ah=varargin{1}; varargin(1)=[];
  else ah=gca; end

  if isempty(varargin)
     helpthis, if nargout
     wbdie('invalid usage'), end, return
  end

  hm=find_all_marks(ah);

  getopt('INIT',varargin);
     if getopt('--fgz'), mv_to_frontZ(ah,hm);
     elseif getopt('--rm'), delete([hm{:}]); end
  getopt('check_error');

end

% -------------------------------------------------------------------- %

function hm=find_all_marks(ah)
   na=numel(ah); hm=cell(na,2);
   for r=1:2
      if r==1, t='line'; else t='text'; end
      for l=1:na
         lh=findall(ah(l),'type',t);
         tt=get(lh,'tag');

         I=regexp(tt,'^\w?mark$');
         for i=1:numel(I), if isempty(I{i}), I{i}=0; end, end
         I=[I{:}]; i=find(I); hm{l,r}=lh(i)';
      end
   end
end

% -------------------------------------------------------------------- %
% move *mark data to the front by assigning ZData

function mv_to_frontZ(ah,hm)
   for l=1:numel(ah);
      zmax=max([ get(ah(l),'ZLim'), get(ah(l),'CLim') ]);
      if zmax<0, zmax=0; end

      for h=hm{l,1}
         set(h,'ZData',repmat(zmax,size(get(h,'XData'))));
      end

      for h=hm{l,2}
         set(h,'Units','data'); p=get(h,'Position'); p(3)=zmax;
         set(h,'Position',p);
      end
   end
end

% -------------------------------------------------------------------- %

