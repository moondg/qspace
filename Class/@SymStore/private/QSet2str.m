function qs=QSet2str(q,varargin)
% function qs=QSet2str(q [opts])
%
%    get formatted q-label string.
%    this routine deals with a specific non-abelian symmetry
%    and hence is specific to RC_STORE (hence a subroutine
%    to the SymStore class!)
%
%   '-f*'  print QStr that can be used with LoadCStore (*-format)
%   '-f;'  print QStr that can be used with LoadCStore (;-format)
%
% Wb,Jan07,16

  getopt('init',varargin);
     if getopt('-f*'), cflag=1;
     elseif getopt('-f;'), cflag=2; else cflag=0; end
  getopt('check_error');

  r=1;
  if ~isempty(regexp(q.type,'^SU'))
      t=regexprep(q.type,'[()]','');
      r=str2num(t(3:end))-1;
  elseif ~isempty(regexp(q.type,'^Sp'))
      t=regexprep(q.type,'[()]','');
      r=str2num(t(3:end))/2;
  end

  n=length(q.qdir);
  l=length(q.qset);

     if n*r~=l, q, error('Wb:ERR',...
     '\n   ERR invalid qset data (%g*%g = %g) !?',n,r,l); end

  qs=cellstr(qmat2char(reshape(q.qset,r,[]))')';

  k=find(diff(q.qdir)); nk=numel(k);
  if nk>1, error('Wb:ERR','\n   ERR invalid qdir !?'); end

  qs(2,:)={''};

  if cflag==1, if nk
     qs(end+1,k+1:end)={'*'}; end
     qs(end+1,1:end-1)={','};
  elseif cflag==2
     qs(2,1:end-1)={','}; if nk
     qs{2,k}=';'; end
  else
     qs(2,1:end-1)={','};
     if nk
        if cflag, s=';'; else s=' | '; end
        qs{2,k}=s;
     end
  end

  qs=[qs{:}];

end

