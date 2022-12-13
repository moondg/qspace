function ndiff=finddiffstruct(a,b,varargin)
% function ndiff=finddiffstruct(A,B)
% 
%    Find fields that differ in between
%    structures or objects A and B.
% 
% Options
% 
%    '-v'  verbose mode
%    '-r'  recursive mode in nested structures
%          i.e. enter fields that contain structs if different
% 
% Wb,Aug08,16

  getopt('init',varargin);
     vflag=getopt('-v');
     rflag=getopt('-r');
     aflag=getopt('-a');
     pre=getopt('--prefix','');
  getopt('check_error');

  if nargin<2 || nargout>1
     helpthis, if nargin || nargout
     wberr('invalid usage'), end, return
  end

  if isobject(a), a=struct(a); end
  if isobject(b), b=struct(b); end

  if ~isstruct(a) || ~isstruct(b)
     helpthis, if nargin || nargout
     wberr('invalid usage'), end, return
  end

  ipre=find(pre=='.');
  lpre=numel(ipre);
  if lpre
       pre_=pre(ipre(end):end);
  else pre_=pre; end

  fa=fieldnames(a);
  fb=fieldnames(b); nl=repmat(' |  ',1,lpre);

  ff=intersect(fa,fb); ndiff=0;

  if numel(ff)~=numel(fa), ndiff=ndiff+1;
     xa=setdiff(fa,fb); s=sprintf(', %s',xa{:});
     fprintf(1,'%s-> got %g extra field(s) in A: %s.%s\n',nl,pre,s(3:end));
  end
  if numel(ff)~=numel(fb)
     ndiff=ndiff+1; if ndiff==1, fprintf(1,'X \n'); end
     xb=setdiff(fb,fa); s=sprintf(', %s',xa{:});
     fprintf(1,'%s-> got %g extra field(s) in B: %s.%s\n',nl,pre,s(3:end));
  end

  na=numel(a);

  if ~isequal(size(a),size(b))
     sa=sprintf('x%g',size(a));
     sb=sprintf('x%g',size(b)); fprintf(1,['%s-> got %g got ' ... 
       'size mismatch with B: %s <> %s\n'],nl,pre,sa(2:end),sb(2:end));
     na=0; ndiff=ndiff+1;
  elseif rflag, o={'-r'};
     if vflag, o{end+1}='-v'; end
     if aflag, o{end+1}='-a'; end
  end

  e2=ndiff;
  for k=1:na, if ~aflag && e2~=ndiff
      fprintf(1,'%s   ...\n',nl); break; end
  for i=1:numel(ff)
     xa=getfield(a,{k},ff{i});
     xb=getfield(b,{k},ff{i});
     if ~isequal(xa,xb)
         ndiff=ndiff+1;
         if rflag && (isobject(xa) || isstruct(xa)) ...
                  && (isobject(xb) || isstruct(xb))
            if vflag, if ndiff==1 && isempty(pre), fprintf(1,'\n'); end
                 fprintf(1,'%s \\ .%s differs (struct)\n',nl,ff{i});
            else fprintf(1,'%s \\ .%s\n',nl,ff{i});
            end
            if na>1
                 pre2=sprintf('(%g/%g).%s',pre_,k,na,ff{i});
            else pre2=[pre_ '.' ff{i}]; end

            ndiff=ndiff+finddiffstruct(xa,xb,'--prefix',pre2,o{:});
         else
            if na>1
                 pre2=sprintf('%s(%g/%g).%s',pre_,k,na,ff{i});
            else pre2=[ '.' ff{i}]; end
            if vflag
                 fprintf(1,'%s-> %s differs\n',nl,pre2);
            else fprintf(1,'%s-> %s\n',nl,pre2);
            end
         end
     end
  end
  end

  if vflag
     if    ~ndiff, fprintf(1,'\n   got matching fields.\n');
     elseif isempty(pre), fprintf(1,' \n'); end
  end

  if ~nargout, clear ndiff; end

end

