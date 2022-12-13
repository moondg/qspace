function untagf(varargin)
% Function untagf(varargin)
%
%    unset tag of current figure
%
% Options
%
%   '-a'  unset tags from all figures
%   '-q'  quite mode
%
% Wb,Dec11,06

  getopt('init',varargin);
     aflag=getopt('-a'  );
     qflag=getopt('-q'  );
  varargin=getopt('get_remaining');

  if aflag,                  h=findall(groot,'type','figure');
  elseif ~isempty(varargin), h=cat(1,varargin{:});
  else                       h=gcf;
  end

  k=0;
  for i=1:length(h)
     t=get(h(i),'tag'); if isempty(t), continue; end
     set(h(i),'tag',''); k=k+1;
     set(h(i),'Name', ['[' get(h(i),'name') ']']);
  end

  if qflag || ~aflag, return; end

  if k==0
     fprintf(1,'   no figure tags set.\n');
  else
     if k~=1, s='s'; else s=''; end
     fprintf(1,'   untagged %d figure%s.\n',k,s);
  end

end

