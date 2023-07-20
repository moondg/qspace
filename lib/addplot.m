function [hh]=addplot(xd0,yd0,varargin)
% Function addplot(xd,yd,varargin)
%
%    add data to plot while checking linear/log-scale plot
%    for log-scale, negative data is shown in dashed format.
%
% Options
%
%  'cfac',.. light / blurred color instead of dashed (cfac>0 and <=1)
%            cfac=0 -> just dashed
%            cfac=1 -> data disappears
%
% See also loglog2.m
% Wb,Apr10,12

% adapted from loglog2.m

  getopt('init',varargin);
     cfac=getopt('cfac',0);
  varargin=getopt('get_remaining');

  if ~isscalar(cfac) || ~isnumeric(cfac) || cfac<0 || cfac>1
  wbdie('invalid usage'); end

  n=length(xd0); s=size(yd0);
  if find(size(xd0)>1)>1, wbdie('invalid x-data'); end
  if numel(s)>2 || all(s~=n), wbdie('invalid y-data'); end
  if s(1)~=n, yd0=yd0.'; end

  xlin=isequal(get(gca,'XScale'),'linear');
     if xlin, ii={ 1:n };
     else ii={ find(xd0>0), find(xd0<0) }; end
  ylin=isequal(get(gca,'YScale'),'linear');
  hold on

  for i=1:numel(ii)
     xd=nan(size(xd0)); xd(ii{i})=xd0(ii{i}); if i>1, xd=-xd; end

  for k=1:size(yd0,2), yd=yd0(:,k);
     if ylin, jj={ 1:n };
     else jj={ find(yd>0), find(yd<0) }; end
     for j=1:numel(jj), yd=yd0(:,k);
        yd=nan(n,1); yd(jj{j})=yd0(jj{j},k); if j>1, yd=-yd; end

        hh=plot(xd,yd,varargin{:},'tag','addplot'); HH{i,j}=hh;
        if i>1 || j>1
           if ~cfac
              set(hh,'LineSt','--','tag','addplot::neg');
           else
              if cfac<1
                 for h=reshape(hh,1,[])
                    c=get(h,'Color'); c=(1-cfac)*(c-1) + 1;
                    set(h,'Color',c,'LineW',2);
                 end
              else set(hh,'Visible','off'); end
           end
        end
     end

  end,end

  hh=cat(1,HH{:});
  mv2front(HH{1});

  if ~nargout, clear hh; end

end

