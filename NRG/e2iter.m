function [n,istr]=e2iter(E,Lambda,z)
% Usage: [n,istr]=e2iter(E,Lambda,z)

  if nargin<1 || nargin>2
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if nargin<2, z=0;
     global param
     if isfield(param,'Lambda')
        Lambda=param.Lambda;
        if isfield(param,'z') && ~isempty(param.z), z=param.z; end
     else
        getcaller Lambda z
        if isempty(z), z=0; end
        if isempty(Lambda)
        error('Wb:ERR','failed to obtain Lambda from caller'); end
     end
  elseif nargin<3, z=0; end

  a=(Lambda^(z-1))*(Lambda-1)/log(Lambda);
  n=1-2*log(E/a)/log(Lambda);

  if nargout>1
     if z==0
        istr='\omega_n = [(1-\Lambda^{-1})/log(\Lambda)] \Lambda^{-n/2}';
     else
        istr='\omega_n = [\Lambda^{z-1}(\Lambda-1)/log(\Lambda)] \Lambda^{-n/2}';
     end
  end

end

