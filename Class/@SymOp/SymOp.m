function S=SymOp(varargin)
% function S=SymOp([type,] istr, Sop)
% Wb,Nov17,11

  if nargin==0 || nargin==1 && isempty(varargin{1})
     S=struct('istr',{},'op',{},'hc',{},'type',{});
     S=class(S,'SymOp'); return
  end

  v1=varargin{1};

  if nargin==1
     if isa(v1,'SymOp'), S=v1;
     elseif isnumeric(v1) && numel(v1)<=3
        s=v1; if numel(s)==1, s=[1 s]; end
        q=repmat({[]},s);
        S=struct('istr',q,'op',q,'hc',q,'type',q);
        S=class(S,'SymOp');
     end
     return
  end

  v2=varargin{2};

  if nargin==2 && ...
     isnumeric(v1) && numel(v1)==1 && isnumeric(v2) && numel(v2)==1 
     q=repmat({[]},v1,v2);
     S=struct('istr',q,'op',q,'hc',q,'type',q);
     S=class(S,'SymOp');
     return
  end

  if nargin && ischar(v1) && numel(v1)==1, t=lower(v1);
     if isempty(find(t=='+-z'))
        wbdie('invalid type=%s !?',v1);
     end
     varargin=varargin(2:end);
  else t=''; end

  narg=numel(varargin); e=0;
  if narg<2 || narg>3, e=1;
  else
     if isempty(varargin{1}) && isempty(varargin{2}), hc={};
     else hc=[];
        if ~ischar(varargin{1}) || ~isnumeric(varargin{2}), e=2; end
     end
     if ~e && narg>2
        if ~isempty(t), e=3;
        elseif ~isequal(varargin{3},'-disc'), e=4;
        else t='disc'; end
     end
  end

  if e
     eval(['help ' mfilename]), if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  S=struct('istr',varargin{1},'op',varargin{2},'hc',hc,'type',t);
  S=class(S,'SymOp');

end

