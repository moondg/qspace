function ncols=getcols(nc_)
% function ncols=getcols([nc_])
%
%    Get terminal width (`number of columns').
%    where nc_ specifies default if routine cannot determine width
%    e.g., in non-terminal mode such as batch mode.
%
% AW 2011-2025.

% added input default nc_  // Wb,Aug01,25

  persistent nc

  if nargin && isequal(nc_,'-t')
     if nargout, ncols=nc;
     else wblog(1,'TST','got ncols=%s',num2str2(nc)); end
     return
  end

  if ~isempty(nc) && (~nargin || isequal(nc,nc_))
     ncols=nc;
  else
     if ~nargin || nc_<60 || nc_>180, nc_=90; end
     if isdesktop>1
        [e,s]=system('command -v getcols.pl');
        if ~e
           q=system('getcols.pl -q 2>/dev/null'); if isempty(q), q=-1; end
           if q<64
                wblog('WRN','got ncols=%g -> %g',q,nc_);
           else nc_=q; end
        end
     end
     if ~nargin || isempty(nc), nc=nc_; end
     ncols=nc_;
  end

end

