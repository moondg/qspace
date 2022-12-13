function odir=itags_to_odir(t)
% function odir=itags_to_odir(itags)
%
%    specifies location of up to one out-index (if more than
%    one out-index is encountered, error is issued).
%    E.g. for standard rank-3 A-tensors, the value of odir
%    indicates the direction of the orthonormalization
%    with reference to the index position in t:
%
%        1  got RL orthonormalized (i.e. 1st index has '*', i.e. is out)
%        2  got LR orthonormalized (i.e. 2nd index has '*', i.e. is out)
%        0  for current A-tensor   (i.e. got no out index)
%       -1  if not yet initialized (got empty itags)
%
% Wb,Jan15,15

% formerly associated with Hamilton1D/private // Wb,May15,17

  if isempty(t), odir=-1; return; end

  if ischar(t)
     t=strread(t,'%s','whitespace',' ,;|\n\r\t')';
  end

  odir=findstrc(t,'\*$');
  nout=numel(odir);

  if numel(t)==3

     if ~nout, odir=0;
     elseif nout>1, t, error('Wb:ERR',['\n   ERR' ...
        'got unexpected itags (got %g out-indizes !?)'],nout);
     end

  elseif numel(t)==4

     if ~nout || nout>2 || odir(end)~=4, t, error('Wb:ERR', ...
        '\n   ERR got unexpected itags (got %g out-indizes !?)',nout);
     end
     if nout==1, odir=0;
     else odir=odir(1); end

  else t, error('Wb:ERR',...
    '\n   ERR got unexpected itags of rank %g !?',numel(t));
  end

end

