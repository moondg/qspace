function P=subsasgn(P,S,val)
% SUBSASGN
%
% Reference http://www.cs.ubc.ca/~murphyk/Software/matlabObjects.html
% Classes and objects in matlab: the quick and dirty way
% Kevin Murphy, 2005

  if ~isa(P,'PSet'), P=PSet(P); end

  f=S(1).subs;
  if isequal(S(1).type,'.') && ischar(f) && ~builtin('isfield',P,f)
     for i=1:length(P.vars)
        if isequal(P.vars{i},f)
           if length(S)>1
                P.data{i}=subsasgn(P.data{i},S(2:end),val);
           else P.data{i}=val; end
           P=update_size(P);
           return
        end
     end
     if length(S)==1 && ischar(f)
        wblog(1,'WRN','adding new parameter field %s',f);
        P.vars{end+1}=f; P.data{end+1}=val;
        P=update_size(P);
        return
     else error('Wb:ERR','invalid usage'); end
  end

  P=builtin('subsasgn',P,S,val);

end

