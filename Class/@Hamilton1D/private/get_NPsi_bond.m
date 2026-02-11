function [NPsi,NPsi2]=get_NPsi_bond(Psi)
% function [NPsi,NPsi2]=get_NPsi_bond(Psi)
% 
%    return the number of global multiplets in Psi
%     * assuming bond-picture using lr[g] index order convention.
%     * NPsi=0 if Psi is rank-2 (while Psi still represents
%       a global scalar multiplet!)
%     * NPsi>=1 if Psi is rank-3.
%    NPsi2 returns the corresponding number of states.
% 
% See also get_NPsi_site.m
% Wb,Aug12,16

  if ~isobject(Psi), Psi=QSpace(Psi); end
  qdir=getqdir(Psi,'-s');

  if isequal(qdir,'++'), NPsi=0; NPsi2=1;
  elseif isequal(qdir,'++-')
     d=getDimQS(Psi); NPsi=d(1,end); NPsi2=d(end);
     t=Psi.info.itags{end};
     if isempty(regexpi(t,'Psi'))
        wblog('WRN','got empty itag ''%s'' for Psi',t);
     end
  else wbdie('invalid usage (got qdir=''%s'')',qdir);
  end

end

