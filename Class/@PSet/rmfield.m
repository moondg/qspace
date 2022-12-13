function P=rmfield(P,f)
% Function i=rmfield(P)
% Wb,Jun16,09

  if ~iscell(f), f={f}; end

  nf=length(f); mark=repmat(-1,1,nf);

  for k=1:nf
     for i=1:length(P.vars)
     if isequal(P.vars{i},f{k}), mark(k)=i; break; end
  end

  if any(mark<0), f{find(mark<0)}
  error('Wb:ERR','field(s) not found'); end

  i=1:P.r; i(mark)=[];

  P.vars=P.vars(i);
  P.data=P.data(i);

  P.t=[]; P.i=0; P=update_size(P);

end

