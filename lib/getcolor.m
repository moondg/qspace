function c=getcolor(varargin)
% Function: getcolor(i,[cmap])
%
%    get i-th color of color order of current axis set
%    taken modulo if i exceeds colormap range
%    (i may also be a vector).
%
% Options
%
%    'cm',..  explicitly specify user defined color order
%    --cm     use default colormap instead of colororder
%    --flip   flip color order upside down
%    --plot   plot color order
%
% See also recolor.m Wb,Apr20,07
% Wb,Jan18,08

  if nargin==1, q=varargin{1};
     if numel(q)==3 && all(q>=0) && all(q<=1), c=q; return; end
     if ischar(q) && numel(q)==1
        if ~isempty(find('brmcyk'==q)), c=q; return; end
        if isequal(q,'g'), c=[0 0.5 0]; return; end
     elseif all(ishandle(q) & ~isnumeric(q))
        c=get(q,'Color'); return
     end
  end

  getopt('init',varargin);
      flip =getopt('--flip');
      cmap =getopt('--cm'  );
      cm   =getopt('cm', []);
      pflag=getopt('--plot');
  varargin=getopt('get_remaining');
  narg=length(varargin);

  if pflag
     plot_ColorOrder(gca);
     if ~narg, return; end
  end

  if ~narg || narg>2
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'); end
     return
  end

  ic=varargin{1};
  varargin(1)=[]; narg=length(varargin);

  if ischar(ic), c=ic;
     if ~isempty( regexp(c,'^gray[0-9][0-9]$') )
        c=repmat(0.01*str2num(c(5:end)),1,3);
     else
        switch c
           case 'brown',     c = [ 0.88  0.68  0.3  ];

           case 'blue',      c = [ 0.1   0.5   1    ];
           case 'lightblue', c = [ 0.5   0.75  1    ];

           case 'orange',    c = [ 1     0.5   0    ];
           case 'orange.in', c = [ 0.88  0.68  0.34 ];

           case 'g',         c = [ 0 .5 0 ];
           case 'y2',        c = [ .85 .85 0 ];

        otherwise
           disp(c)
           wbdie('unknown user color specification');
        end
     end

  else
     if numel(ic)==3 && all(ic<=1) && all(ic>=0)
        c=ic; return
     end

     if isempty(cm)
        if cmap, cm=colormap;
        else cm=get(gca,'ColorOrder'); end
     else
        if ~isnumeric(cm) || size(cm,2)~=3
        wbdie('invalid input color map'); end
     end

     nc=size(cm,1); ic=ic(:);

     if flip
          j=mod(nc-ic,nc)+1;
     else j=mod(ic-1, nc)+1; end

     if isequal(j,round(j)), c=cm(j,:);
     else
        j_=floor(j); x=repmat(j-j_,1,3);
        j =ceil(j);  j(find(j>nc))=1;

        c=cm(j_,:).*(1-x) + cm(j,:).*x;
     end
  end

if ~narg, return; end

  il=varargin{1}(:);

  s=size(c);
  c=mat2cell(c,ones(s(1),1),s(2));

  lo={'-','--','-.',':'}; nl=numel(lo); il=mod(il-1, nl)+1;

  if numel(il)==1
     c(:,2)={il};
  elseif numel(il)==numel(ic)
     c(:,2)=mat2cell(il,ones(s(1),1),1);
  else
     wbdie('invalid specification of line style');
  end

  for i=1:s(1), c{i,2}=lo{c{i,2}}; end

  if s(1)==1
     c={'Color',c{1},...
    'MarkerEdgeColor',c{1},'MarkerFaceColor',c{1},'LineSt',c{2}};
  end

end

% -------------------------------------------------------------------- %
function plot_ColorOrder(h0)

   cc=get(h0,'ColorOrder');
   k=get(h0,'ColorOrderIndex');

   ah=smaxis(1,1,'tag',mfilename); setax(ah(1,1))

   x=[.1 .9]; x1=x(1)-0.05*diff(x);

   for i=1:size(cc,1), y=i+[-.4 .4];
       h=patch(x([1 2 2 1 1]), y([1 1 2 2 1]),cc(i,:),'EdgeColor',cc(i,:));
       if i==k
          set(h,'EdgeColor','k','LineW',4);
       end
       text(x1,i,sprintf('%g',i),'FontSize',14);
   end
   axis off ij
end

% -------------------------------------------------------------------- %

