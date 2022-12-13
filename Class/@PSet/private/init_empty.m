function P=init_empty(vars,data)
% function P=init_empty([vars,data])
% Wb,Mar13,19

  if nargin<1, vars={}; data={};
  elseif nargin<2, data={}; end
  n=numel(vars);

  if nargin
     if ~iscell(vars) || ~iscell(data), error('Wb:ERR',...
     '\n   ERR invalid usage (cell array required)');
     elseif n~=numel(data), error('Wb:ERR',...
     '\n   ERR invalid usage (len=%d/%d)',n,numel(data));
     end
  end

  P=struct(...
      'vars',{vars}, 'data',{data}, ...
      'r',0, ... % `rank'
      'n',0,  ...
      's',[], ...
      'i',0,  ...
      't',[], ...
      'info',[]...
  );

  if n
     s=zeros(1,n); for i=1:n, s(i)=numel(data{i}); end
     P.s=s; P.r=n;
  end

end

