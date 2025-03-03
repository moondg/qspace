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

  vflag=1;
  getopt('init',varargin);
     aflag=getopt('-a');
     if getopt('-q'), vflag=0; end
  args=getopt('get_remaining'); nargs=numel(args);

  if aflag
     if nargs, wbdie('invalid usage'); end
     h=findall(groot,'type','figure');
  elseif nargs
       h=[args{:}];
  else h=gcf;
  end

  k=0;
  for i=1:numel(h), t=get(h(i),'tag');
     if ~isempty(t)
        set(h(i),'tag',''); k=k+1;
        set(h(i),'Name', ['[' get(h(i),'name') ']']);
     end
  end

  if vflag
     if aflag
        if k, if k==1, s=''; else s='s'; end
             fprintf(1,'\n   untagged %d figure%s\n\n',k,s);
        else fprintf(1,'\n   (no figure tags set)\n\n');
        end
     elseif ~nargs
        s=sprintf('current Fig. %g',double(gcf));
        if k
             fprintf(1,'\n   untagged %s\n\n',s);
        else fprintf(1,'\n   %s already untagged\n\n',s);
        end
     end
  end

end

