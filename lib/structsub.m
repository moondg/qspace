function Iout=structsub(S,varargin)
% function S=structsub(S,fields)
% Wb,Dec11,20

  if ~nargin || ~isstruct(S) || numel(S)~=1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     vflag=getopt('-v');
     xflag=getopt('-x');
  ff=getopt('get_remaining');

  if numel(ff)==1 && iscell(ff{1}), ff=ff{1}; end
  nf=numel(ff);
  for i=1:nf
     if ~ischar(ff{i}), wbdie('invalid usage'); end
  end

  Iout=struct;
  if xflag, f0=fieldnames(S);
     for i=1:nf
        for j=find(~cellfun(@isempty,regexp(f0,ff{i})))'
           Iout=setfield(Iout,f0{j},getfield(S,f0{j}));
        end
     end
  else
     for i=1:nf
        if ~isfield(S,ff{i}), wbdie('invalid field ''%s''',ff{i}); end
        Iout=setfield(Iout,ff{i},getfield(S,ff{i}));
     end
  end

end

