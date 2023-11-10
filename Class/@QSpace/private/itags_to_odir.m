function odir=itags_to_odir(t,lflag)
% function odir=itags_to_odir(itags [,opts])
%
%    Default behavior (i.e., With no options specified):
%
%    This routine returns the location of up to one out-index
%    (if more than one out-index is encountered, error is issued).
%    E.g. for standard rank-3 A-tensors, the value of odir
%    indicates the direction of the orthonormalization
%    with reference to the index position in t:
%
%        1  got RL orthonormalized (i.e. 1st index has '*', i.e. is out)
%        2  got LR orthonormalized (i.e. 2nd index has '*', i.e. is out)
%        0  for current A-tensor   (i.e. got no out index)
%       -1  if not yet initialized (got empty itags)
%
% Options
%
%  -l   long format (show direction for each index, returning
%       +1 for ingoing
%       -1 for outgoing indices (having trailing conjugate flag *)
%
% Wb,Jan15,15

% formerly associated with Hamilton1D/private // Wb,May15,17

  if nargin>1
     if ~isnumeric(lflag)
        if ~isequal(lflag,'-l'), wbdie('invalid usage'); end
        lflag=-1;
     end
  else lflag=0; end

  if isempty(t)
     if lflag, odir=[];
     else odir=-1; end
     return
  end

  if ischar(t)
     t=strread(t,'%s','whitespace',' ,;|\n\r\t')';
  end

  nt=numel(t);
  if any(lflag)
     if isequal(lflag,-1), I=1:nt;
     else I=lflag;
        if any(I<1 | I>nt), wbdie(...
          'invalid usage (index out of rank:%s / %d)',sprintf(' %d',I),nt);
        end
     end
     n=numel(I); odir=ones(1,n);
     for i=1:n
        if ~isempty(regexp(t{I(i)},'\*$')),  odir(i)=-1; end
     end
     return
  end

  odir=findstrc(t,'\*$');
  nout=numel(odir);

  if nt==3

     if ~nout, odir=0;
     elseif nout>1, t
        wbdie('got unexpected itags (got %g out-indizes)',nout);
     end

  elseif nt==4

     if ~nout || nout>2 || odir(end)~=4, t
        wbdie('got unexpected itags (got %g out-indizes!?)',nout);
     end
     if nout==1, odir=0;
     else odir=odir(1); end

  else t
     wbdie('got unexpected itags of rank %g',nt);
  end

end

