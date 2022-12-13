function cm=whot(varargin)
% function cm=whot(n [,opts])
%
%    like colormap hot, but start from white, to hot, to black
%    (fitted to the colormap used by Isaac Castillo)
%
% Wb,Aug28,19

  getopt('init',varargin);
     pnflag=getopt('-pn');
  n=getopt('get_last',64);

  cm=zeros(n,3); xx=linspace(0,1,n);

  cm(:,1)=get_weight(xx,[1/3],...
  [2.96 0],[0.02 0.98]);

  cm(:,2)=get_weight(xx,[1/3, 0.36, 2/3],...
  [0.46 0 0],[0 0.051],[1.94 -0.65],[1.07 -0.07]);

  cm(:,3)=get_weight(xx,[0.3, 2/3, 0.73],...
  [1.10 0 0],[0.36 -0.20 0.13],[0.22 0.01],[3.05 -2.05]);

  if pnflag
       cm=[fliplr(cm);flipud(cm)]; cm=cm(1:2:end,:);
  else cm=flipud(cm);
  end
end

% -------------------------------------------------------------------- %

function c=get_weight(xx,xk,varargin)

   c=zeros(numel(xx),1); xk=[-inf,xk,inf];

   for k=1:numel(xk)-1
      i=find(xx>=xk(k) & xx<xk(k+1));
      if ~isempty(i)
         c(i)=polyval(varargin{k},xx(i));
      end
   end

end

% -------------------------------------------------------------------- %

