function disp_qdim(A,varargin)
% function disp_qdim(A,varargin)
%   example: disp_qdim(A0,2,IS.SOP)
% Wb,Jan26,12

  if nargin>1 && (isequal(varargin{1},'op') || isnumeric(varargin{1}))
       dim=varargin{1}; varargin(1)=[];
  else dim='op'; end

  getopt('init',varargin);
     split=getopt('split',[]);
     P    =getopt('P',[]);
  SOP=getopt('get_last',[]);

  if nargout
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if ~isempty(SOP)
     for i=1:numel(SOP)
        dz(i)=numel(SOP(i).Sz);
     end
     if ~isempty(split) && ~isequal(split,dz), split,dz
        error('Wb:ERR','\n   ERR inconsistent split/SOP');
     else split=dz; end

     sym={ SOP.info }; if ~isempty(P), sym=sym(P); end
     fprintf(1,'\n');
     for i=1:numel(sym), fprintf(1,'   %2d. %s\n',i,sym{i}); end
     fprintf(1,'\n');
  end

  [Q,dd,dc]=getQDimQS(A,dim);

  ns=size(dc,2);
  if isempty(split)
     if ns~=size(Q,2), error('Wb:ERR',...
       '\n   ERR got rank>1 symmetry => must specify split'); end
     split=ones(1,ns);
  end

  Q=mat2cell(Q,size(Q,1),split);

  if ~isempty(P)
     if numel(P)~=ns, error('Wb:ERR',...
       '\n   ERR invalid permutation (expected length %d)',ns); end
     Q=Q(P); dc=dc(:,P);
  end

  Q=cat(2,Q{:});
  [Q,is]=sortrows(Q); dd=dd(is); dc=dc(is,:);

  l=[14 % length of 'symmetry labels' in header
     3*sum(split)+2*(numel(split)-1) ];

  s=''; if l(1)>l(2), s=repmat(' ',1,-diff(l)); end
  L=repmat('-',1,72); L(1:3)=' ';
  L([9, max(l)+[14 23]])='+';

  ss=cell(1,ns);
  for i=1:ns, ss{i}=['; ' repmat(' %2g',1,split(i)) ]; end
  ss=cat(2,ss{:}); ss=[ '%6d. |' ss(3:end), ' > ' s ' |%6g  |' repmat(' %3g',1,ns) '\n'];

  fprintf(1,['\n%s\n    ##  |  %-' num2str(max(l)) 's | %s| %s\n%s\n'], ...
    L,'symmetry labels','D(data)','dd(CGS)',L);

  fprintf(ss,[ 1:size(Q,1); Q'; dd; dc' ]);

  fprintf(1,'%s\n   Total space dimension: %d\n\n',L,...
    sum(prod([dd',dc],2)));

end

