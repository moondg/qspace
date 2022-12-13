function kdir=check_dir(kdir)
% function kdir=check_dir(kdir)
%
%    NB! this function can only be called from /
%    is only within the scope of member functions.
%
% Wb,Apr08,14

% outsourced from updateHam()
% see also QSpace/itags2odir.m // tags: direction

  if ischar(kdir)
     switch kdir
        case {'>>'}, kdir=+1;
        case {'<<'}, kdir=-1;
        otherwise wberr('invalid kdir=''%s''',kdir); 
     end
  else
     if isnumeric(kdir) && numel(kdir)==1 && kdir
        if kdir>0, kdir=1; else kdir=-1; end
     else
        kdir, wberr('invalid orthonormalization direction'); 
     end
  end

end

