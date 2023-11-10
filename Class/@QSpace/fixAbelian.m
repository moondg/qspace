function varargout=fixAbelian(varargin)
% function fixAbelian(A,B,...)
% Options
%
%    -xop3    reduce "scalar" rank-3 operators to regular rank-2.
%
%    -addCGR  add (trivial) CGR space to all-abelian operators.
%    -op3     together with -addCGR: make operator of rank 3 if rank 2.
%
% Wb,Feb20,13

  isopt=zeros(1,nargin);
  for ia=1:nargin
     if ischar(varargin{ia}), isopt(ia)=1; end
  end

  i=find(isopt); IA=1:nargin;
  if ~isempty(i)
     opts=varargin(i); IA(i)=[]; nargs=numel(IA);
     getopt('init',opts);
        xop3=getopt('-xop3');
     getopt('check_error');
  else xop3=0; nargs=nargin; end

  if nargout && nargs~=nargout
     wbdie('invalid usage (must match input to output variables)');
  end

  mark=zeros(size(IA));

  for ia=1:nargs, A=varargin{IA(ia)};
     if ~isempty(A), isa=isAbelian(A); isa=isa(:);
     if any(isa)
        for k=1:numel(A)
           if isa(k) && xop3
              d=getDimQS(A(k));
              if length(d)==3 && d(3)==1, A(k).Q(3)=[]; end
              if numel(A(k).info.itags)==3, A(k).info.itags(3)=[]; end
           end
           if isa(k)>1, A(k).info={}; end
        end
        varargin{IA(ia)}=A; mark(ia)=1;
     end, end
  end

  if nargout, varargout=varargin(IA);
  else
     for ia=find(mark), k=IA(ia); n=inputname(k);
        if isempty(n)
           if numel(IA)==1, varargout=varargin(IA);
           else wbdie(...
             'invalid usage (failed to access name of input variable)');
           end
        else
           assignin('caller',n,varargin{k});
        end
     end
  end

end

