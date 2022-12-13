function b=bsize(varargin)
% function bsize(varargin)
%
%    return size in bytes for each input argument
%
% Wb,Mar03,22

  b=zeros(1,nargin); sflag=[];
  for i=1:nargin
      if ischar(varargin{i}) && isequal(varargin{i},'-s')
         sflag(end+1)=i;
      else
         q=varargin{i}; I=whos('q');
         b(i)=I.bytes;
      end
  end
  if isempty(sflag), return; end
  if numel(sflag)~=1, wblog('WRN','%g -s flags encountered',numel(sflag)); end

  b(sflag)=[]; b=matcell(b);
  for i=1:numel(b), b{i}=num2str2(b{i},'--bytes'); end
  if numel(b)==1, b=b{1}; end

end

