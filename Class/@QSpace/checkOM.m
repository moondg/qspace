function m=checkOM(A,varargin)
% function m=checkOM(A [,opts])
%
%    Check for outer multiplicity in given QSpace
%
% Wb,Feb15,14

  getopt('init',varargin);
     vflag=getopt('-v');
     sflag=getopt('-s');
  getopt('check_error');

  n=numel(A); m=zeros(size(A));

  if sflag
     for k=1:numel(A), Ak=A(k);
        ra=numel(Ak.Q); if ra<=2, continue; end
        rs=getsym(Ak,'-r'); if any(rs) && (ra>3 || any(rs>1)), m(k)=1; end
     end
     return
  end

  for k=1:numel(A), Ak=A(k);
     [qq,I,D]=uniquerows(cat(2,Ak.Q{:}));
     [dd,id,d2]=uniquerows(double(D)'); im=find(dd>1);

     if ~isempty(im), m(k)=d2(im)*dd(im)-numel(im);
        if vflag
           for j=1:numel(im), l=im(j);
              getsub(Ak,[I{id{l}}])
              fprintf(1,['   got multiplicity of %g (%g times) ' ... 
              '[%g/%g records]\n'], dd(l), d2(l), dd(l)*d2(l), size(Ak.Q{1},1));
           end
           fprintf(1,...
           '   got total outer multiplicity of %g/%g records.\n\n', ...
           m(k), size(Ak.Q{1},1));
        end
     end
  end

  if ~nargout && vflag
     if all(m(:)==0), s='got no'; else s='NB! got'; end
     fprintf(1,'\n   %s outer multiplicity.\n\n',s);
     clear m
  end

end

