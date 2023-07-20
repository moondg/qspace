function [NPsi,NPsi2]=get_NPsi_site(Psi,cflag)
% function [NPsi,NPsi2]=get_NPsi_site(Psi [,'-c'])
% 
%    return the number of global multiplets in Psi
%     * assuming single-site picture using LRs[g] index order convention.
%     * NPsi=0 if Psi is rank-3 (while Psi still represents
%       a global scalar multiplet!)
%     * NPsi>=1 if Psi is rank-4.
%    NPsi2 returns the corresponding number of states.
% 
% Option '-c' indicates that Psi represents the current site.
% 
% See also get_NPsi_bond.m
% Wb,Aug18,16

  if ~isobject(Psi), Psi=QSpace(Psi); end
  qdir=getqdir(Psi,'-s');

  if nargin<2, dbase='**+';
     if numel(qdir)>=2
        if qdir(1)~='+' && qdir(2)~='+', wbdie('invalid Psi'); end
        dbase(1:2)=qdir(1:2);
     end
  elseif isequal(cflag,'-c'), dbase='+++';
  elseif isequal(cflag,'-A'), dbase='**+';
     if numel(qdir)>=3
        if qdir(1)==qdir(2), wbdie('invalid Psi'); end
        dbase(1:2)=qdir(1:2);
     end
  else wbdie('invalid Psi'); end

  if isequal(qdir,dbase), NPsi=0; NPsi2=1;
  elseif isequal(qdir,[dbase '-'])
     d=getDimQS(Psi); NPsi=d(1,end); NPsi2=d(end);
     t=Psi.info.itags{end};
     if isempty(regexpi(t,'Psi'))
        wblog('WRN','got itag ''%s'' for Psi !?',t);
     end
  else wbdie('invalid usage (got qdir=''%s'' !?)',qdir); end

end

