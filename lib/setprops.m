function setprops(hh,varargin)
% function setprops(hh [,opts,'prop1',val1,...])
%
%    set options of given object handles
%
% Options
%
%   '-s'      save old value userdata 
%   '-reset'  return to saved values (prop list must be empty)
%
% Wb,Apr16,11

  getopt('init',varargin);
     sflag=getopt('-s');
     reset=getopt('-reset');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if reset
     if narg, helpthis, wbdie('invalid usage'); end

     for i=1:numel(hh), h=hh(i);
        v=getuser(h,'-rm','props'); if isempty(v), return; end
        o=fieldnames(v); n=numel(o); o(2,:)={[]};
        for j=1:n, o{2,j}=getfield(v,o{1,j}); end
        set(h,o{:});
     end
     return
  end

  e=0; if mod(narg,2), e=2; end
  if ~e, for i=1:2:narg,
     if ~ischar(varargin{i}), e=2; break; end; end
  end
  if e
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if isempty(hh), return; end

  if sflag
     o=varargin(1:2:end);
     for i=1:numel(hh), h=hh(i); o(2,:)=get(h,o(1,:));
        v=getuser(h,'props');
        if isempty(v)
             v=struct(o{:});
        else v=setfields(v,o{:}); end
        setuser(h,'props',v);
     end
  end

  set(hh,varargin{:});

end

