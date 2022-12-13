function Q=contract (A, B, ma, mb, permQ);
% Function - contract
% Usage: Q = contract (A, B, ma, mb, permQ);
%
%   general routine to contract certain
%   set of indizes of two tensors A and B
%
%   Usage 2: generalized trace: contract(A, [i1 i2; i3 i4; ...])
%   traces out i1 with i2, i3 with i4, ...
%
% AWb, F.Verstraete

% check generalized trace: contract(A, [i1,i2], [i3,i4], ...)
  if nargin==2
     if ~isequal(size(A),size(B))
     eval(['help ' mfilename]); error('Wb:ERR',''); end

     Q=reshape(A,1,[])*reshape(B,[],1);
     return
  end

  if nargin<4, eval(['help ' mfilename]); return; end

  if iscell(ma), ra=ma{2}; ma=ma{1}; else ra=max(ma); end
  if iscell(mb), rb=mb{2}; mb=mb{1}; else rb=max(mb); end

  sa=size(A); sa(end+1:ra)=1; ra=length(sa);
  sb=size(B); sb(end+1:rb)=1; rb=length(sb);

  Ia=1:ra; Ia(ma)=[];
  Ib=1:rb; Ib(mb)=[];

  AA = permute(A,[Ia ma]);
  BB = permute(B,[Ib mb]);

  if ~isequal(sa(ma),sb(mb))
     error('Wb:ERR','misfit in contracted index range [%s; %s]',...
     vec2str(sa(ma)), vec2str(sb(mb)));
  end

  Q = reshape(AA,[prod(sa(Ia)) prod(sa(ma))]) * ...
      reshape(BB,[prod(sb(Ib)) prod(sb(mb))]).';

  snew=[sa(Ia) sb(Ib)]; snew(end+1:2)=1;
  Q = reshape(Q,snew);

  if nargin==5, if ~isempty(permQ), Q=permute(Q,permQ); end; end

end

