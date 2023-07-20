function rmuser(varargin)
% Function rmuser([hh,] tag)
%
%    Remove field tag of userdata of current axes (or hh)
%    to value specified.
%
%    See also getuser(<field>,'-rm'); % also removes field.
%
% Wb,Apr17,09

  if nargin && isnumeric(varargin{1})
     ah=varargin{1}; varargin(1)=[];
     narg=length(varargin);
  else
     ah=gca; narg=nargin;
  end

  if narg<1
     if ~nargout, eval(['help ' mfilename]); end
     return
  end

  ah=reshape(ah,[],1); n=length(ah);
  for i=1:n, h=ah(i);

     u=get(h,'UserData');
     if ~isstruct(u), wbdie('UserData not of type structure.'); end

     for i=1:narg, f=varargin{i};
        if ~ischar(f), wbdie('invalid usage'); end
        if isfield(u,f), u=rmfield(u,f);
        else wblog('WRN','no field %s in UserData',f); end
     end
     set(h,'UserData',u);
  end

end

