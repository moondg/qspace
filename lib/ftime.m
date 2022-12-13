function ts=ftime(varargin)
% Function: timestamp = ftime(file_name);
% Wb,Mar03,06

  if nargin==0
  eval(['help ' mfilename]); return, end

  ts=zeros(nargin,1);

  for i=1:nargin
      f=dir(varargin{i});
      if length(f)~=1
         ts(i)=-1; continue
      end

      ts(i)=datenum(f.date);
  end

return

