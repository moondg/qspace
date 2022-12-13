function i=checkSameSpace(A,ia,B,ib,varargin)
% function i=checkSameSpace(A,ia,B,ib)
%
%    check whether QSpaces A and B share the same vector spaces
%    (simple check in terms of symmetry blocks and respective
%    multiplet [and Clebsch-Gordan dimensions, if applicable]).
%
%    This routine uses getQDimQS(A,ia) [and similar for B]
%    hence ia accepts the same index specification as getQDimQS
%    that is, a specific dimension or 'op'(erator).
%
% Wb,Jan12,12

  if nargin<4
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     vflag=getopt('-v');
  getopt('check_error');

  [qa,da,ca]=getQDimQS(A,ia);
  [qb,db,cb]=getQDimQS(B,ib);

% if ischar(ia), k=1; else k=ia; end
%    [xa,Ja]=matchIndex(A.Q{k},qa,'-s');
% if ischar(ib), k=1; else k=ib; end
%    [xb,Jb]=matchIndex(B.Q{k},qb,'-s');
% keyboard

  if isempty(qa) || isempty(qb)
     if vflag, fprintf(1,'\n   WRN got empty QSpace(s)\n'); end
     if nargout, i=i; else clear i; end
     return
  end

  [Ia,Ib,I]=matchIndex(qa,qb);
  if isempty(Ia)
     if vflag, fprintf(1,'\n   QSpaces do not share any symmetries\n'); end
     if nargout, i=0; else clear i; end
     return
  end

  if ~isequal(da(Ia),db(Ib))
     if vflag, n=0;
        fprintf(1,'\n   QSpaces differ in (multiplet) dimension\n');
        for i=1:numel(Ia)
            if da(Ia(i))~=db(Ib(i)), n=n+1;
               if n>5, fprintf(1,'   ...\n'); break; end
               fprintf(1,'   A{%2d}: %4d    B{%2d}: %4d\n',...
               Ia(i),da(Ia(i)), Ib(i),db(Ib(i)));
            end
        end
     end
     if nargout, i=0; else clear i; end
     return
  end

  if ~isequal(ca(Ia,:),cb(Ib,:))
     if vflag, n=0;
        fprintf(1,'\n   QSpaces differ in Clebsch-Gordan dimension\n');
        fmt=' %3d';
        fmt=[ '  %-3s' ...
            repmat(fmt,1,size(qa,2)) '  %10s | ' ...
            repmat(fmt,1,size(qb,2)) '  %10s %5s\n'
        ];
        for i=1:numel(Ia)
            if isequal(ca(Ia(i),:), cb(Ib(i),:))
               s=cell(1,2); else s={' =>', ' <='}; end
            fprintf(1,fmt, s{1}, ...
            qa(Ia(i),:), ['(' vec2str(ca(Ia(i),:),'sep','x') ')'], ...
            qb(Ib(i),:), ['(' vec2str(cb(Ib(i),:),'sep','x') ')'], s{2});
        end
     end
     if nargout, i=0; else fprintf(1,'\n'); clear i; end
     return
  end

  if nargout, i=1; else, clear i
     fprintf(1,'   QSpaces have matching vector spaces.\n');
  end

end

