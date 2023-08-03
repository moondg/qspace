function kdir=check_dir(kdir)
% function kdir=check_dir(kdir)
%
%    Check sweep direction and return / convert to numeric value.
%
% Wb,Apr08,14

% see also QSpace/itags2odir.m // tags: direction

  if ischar(kdir)
     switch kdir
        case {'>>'}, kdir=+1;
        case {'<<'}, kdir=-1;
        otherwise wbdie('invalid kdir=''%s''',kdir); 
     end
  else
     if isnumeric(kdir) && numel(kdir)==1 && kdir
        if kdir>0, kdir=1; else kdir=-1; end
     else
        kdir, wbdie('invalid orthonormalization direction'); 
     end
  end

end

