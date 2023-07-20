function [isd,estr]=isdiag(A,dflag)
% function [isd,estr]=isdiag(A)
%
%    check whether QSpace is in diagonal rank-2 format.
%
% Return value isd:
%
%    1   if block-diagonal in symmetry space
%    2   if, in addition, all data{} sets are row-vectors, e.g.
%        compressed format as returned bei QSpace::EigenSymmetric()
%    3   same as 2, but if all entries in data{} are column-vectors
%    0   otherwise
%
% Options
%
%  '-d'  If all data blocks are scalars, i.e., of dimension 1x1,
%        by default, this returns isd=1; with '-d' this returns
%        isd=2, instead.
%
% Wb,Apr24,09

  n=length(A.data); isd=-1; estr='';

  r=numel(A.Q); if r
     if r~=2, estr=sprintf(...
        '%s requires rank-2 object (%d).',mfilename,r);
     elseif ~isequal(A.Q{:}), estr=sprintf(... 
        '%s requires block-diagonal operator.',mfilename); end
     if nargout<2 && ~isempty(estr), wbdie(estr); end
  end

  for i=1:n, ai=A.data{i}; s=size(ai);
     if length(s)~=2 || s(1)~=s(2) && all(s~=1)
        estr=sprintf('expecting symmetric operator!');
        break;
     end
     if isempty(ai), continue; end

     switch isd
       case 2
         if all(s~=1), estr=sprintf(...
           'confusing block dimensions (compressed !??)'); break; end
       case 1
         if s(1)~=s(2)
            if any(s==1), estr=sprintf(...
              'confusing block dimensions (compressed !??)'); break
            else, estr=sprintf(...
              'got rectangular data block (%s)',vec2str(s)); break
            end
         else e=norm(ai-diag(diag(ai)));
            if e>1E-12
               isd=0; break
            end
         end
       case 0
       otherwise
         if any(s>1)
            if s(1)==s(2), isd=1;
            elseif any(s==1), if s(1)<s(2), isd=2; else isd=3; end
            else, isd=0; estr=sprintf(...
              'got rectangular data block (%s)',vec2str(s)); break
            end
         end
     end

     if isd==1
        e=norm(ai-diag(diag(ai)));
        if e>1E-12, isd=0; end
     end
  end

  if ~isempty(estr)
     if nargout<2, error('Wb:ERR',estr); end
     isd=0;
  elseif isd<0
     if nargin>1 && r
        if ~isequal(dflag,'-d')
           error('Wb:ERR','\n   ERR invalid dflag'); end
        isd=2;
     else isd=1; end
  end

end

