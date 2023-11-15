function [isd,istr]=isdiag(A,varargin)
% function [isd,istr]=isdiag(A [,opts])
%
%    check whether QSpace is in diagonal rank-2 format.
%    The return value isd is as follows:
%
%    1   if block-diagonal in symmetry space (here blocks themselves
%        do not need to be diagonal see option -f to also enforce this)
%    2   if on top of (1), all data{} sets are row-vectors,
%        e.g., compressed format as returned bei eigQS()
%    3   same as 2, but if all entries in data{} are column-vectors.
%
%    0   if none of the above.
%
% Options
%
%  '-d'  If all data{} blocks are scalars, i.e., of dimension 1x1,
%        by default, this returns isd=1; with '-d' this returns
%        isd=2, instead.
%
%  '-f'  check for fully diagonal, i.e., also enforce diagonal
%        blocks data{} within numerical noise (1E-12) for return
%        value 1 (if non-diagonal, this returns 0, instead).
%
% Wb,Apr24,09

% [11/14/2012] added option '-f': this turns off full diagonal
% check by default, since this adds considerable overhead
% (in response to Jeongmin Shim's email stating that QSpace/trace
% (which uses isdiag!) is much slower than QSpace/norm which
% directly wraps to normQS.
% [11/15/2012] also removed error in case of not diagonal!
% rather check return value in caller! also wrapped earlier
% error string istr into info string istr.

  getopt('init',varargin);
     dflag=getopt('-d');
     fflag=getopt('-f');
     deps =getopt('deps',1E-12);
  getopt('check_error');
  vflag=(nargout>1); istr={};

% empty A.Q (r=0) can represent emtpy QSpace but also scalar QSpace!
% NB! both are *not inconsistent* with diagonal, e.g., empty QSpace
%     could be diagonal empty operator, such as empty HK or HD.
% --> return isd=1, after including further checks below.
% this is consistent with matlab's isdiag([]) which also returns 1.
  r=numel(A.Q); isd=0;
  if r
     if r~=2
        if vflag, istr=sprintf('got rank-%d object',r); end
        return
     elseif ~isequal(A.Q{:})
        if vflag, istr=sprintf('not block-diagonal in symmetries'); end
        return
     end
  end

  isd=-99; n=numel(A.data);

  for i=1:n, ai=A.data{i}; s=size(ai);
     if numel(s)~=2 || diff(s) && all(s>1)
        if vflag, istr={'of size %s',size_str(ai)}; end
        isd=0; break
     end
     if all(s==1) || isempty(ai), continue; end

     if isd==1
        if diff(s), if vflag
           istr={'mixed full and compressed blocks (%s)',size_str(ai)}; end
           isd=-isd; break
        elseif fflag && s(1)>1, e=norm(ai-diag(diag(ai)));
           if e>deps
              if vflag, istr={'non-diagonal @ %.3g',e}; end
              isd=0; break
           end
        end
     elseif isd>1
        if s(isd-1)>1, if vflag
           istr={'mixed up compressed format (%s)',size_str(ai)}; end
           isd=-isd; break
        end
     elseif any(s>1) && isd<0
        if ~diff(s), isd=1;
           if fflag
              e=norm(ai-diag(diag(ai)));
              if e>deps, isd=0; end
           end
        elseif s(1)==1, isd=2;
        elseif s(2)==1, isd=3;
        else
           wbdie('unexpected size [%s ] (isd=%d)',sprintf(' %d',s),isd);
        end
     else
        wbdie('unexpected size [%s ] (isd=%d)',sprintf(' %d',s),isd);
     end
  end

  if vflag && ~isempty(istr)
     istr=sprintf(['data(%d/%d): ' istr{1} ' / isd=%d'],i,n,istr{2:end},isd);
  end
  if isd==-99
     if dflag && r, isd=2; else isd=1; end
  end

end

% -------------------------------------------------------------------- %
function s=size_str(a)
   s=sprintf('x%d',size(a));
   if ~isempty(s), s=s(2:end); end
end

% -------------------------------------------------------------------- %

