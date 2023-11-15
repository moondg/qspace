function [varargout]=getEQdata(HK,varargin)
% function [E,Q,D]=getEQdata(HK [,opts])
% Wb,Sep22,12

% outsourced from $NRG/nrgphase.m

  if ~nargin || isdiag(HK,'-d')~=2
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  getopt('init',varargin);
     iq=getopt('iq',[]);
  getopt('check_error');

  qq=HK.Q{1}; qq=mat2cell(qq,ones(size(qq,1),1),size(qq,2));
  for i=1:length(qq)
     qq{i}=repmat(qq{i},length(HK.data{i}),1);
  end
  qq=cat(1,qq{:}); if isempty(iq), iq=1:size(qq,2); end

  varargout={
     cat(2,HK.data{:})'
     qq(:,iq)
     getzdim(QSpace(HK),2,'-p','-x')
     qq
  };

  [varargout{1},is]=sort(varargout{1});
  for i=2:numel(varargout)
      varargout{i}=varargout{i}(is,:);
  end

  if nargout<=1
     varargout = { cat(2,varargout{1:3}); };
  elseif nargout==2
     varargout = { cat(2,varargout{1:3}), varargout{end}};
  end

end

