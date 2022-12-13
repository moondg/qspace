function ie=isequal(A,B,varargin)
% function i=isequal(A,B [,opts])
% 
%    overloading isequal and ==
%    where input A or B can also be struct( QSpace object ).
%
%    Note that the usage A==B is equivalent to isequal(A,B),
%    i.e. this assumes the default setting without any options.
%
% Options
%
%   '-l'  lenient (ignore difference in itags)
%   '-f'  enforce all fields to be the same (no lenience at all)
%
% Wb,Nov09,09 ; Wb,Jun28,21

  ie=false;

% use leniency on irrelevant info fields
% since the operator == does not permit any options to be specified
  lflag=1;

  q=[ isobject(A), isobject(B) ];
  if any(~q)
     if (~q(1) && ~isstruct(A)) || (~q(2) && ~isstruct(B)), return; end
     if q(1), A=struct(A); end
     if q(2), B=struct(B); end % by builtin('isequal',A,B) below
  elseif ~isequal(size(A), size(B)), return; end

  if nargin>2
     getopt('init',varargin);
        if getopt('-l'), lflag=2;
        elseif getopt('-f'), lflag=0;  end
     getopt('check_error');
  end

  nA=numel(A);
  if lflag
     if ~isfield(A,'info') || ~isfield(B,'info'), lflag=0; end
     for k=1:nA
         if ~isfield(A(k).info,'cgr'), lflag=0; break; end
         if ~isfield(B(k).info,'cgr'), lflag=0; break; end
     end

     for k=1:nA
        if isfield(B(k).info,'ctime')
           A(k).info.ctime=B(k).info.ctime;
        elseif isfield(A(k).info,'ctime')
           B(k).info.ctime=A(k).info.ctime;
        end
     end
  end

  if lflag
     for k=1:nA
        a=A(k).info.cgr;
        b=B(k).info.cgr;
        if ~isequal(size(a),size(b)), return; end

        if isfield(a,'cid') && isfield(b,'cid'), m=0;
           for i=1:numel(a)
              q={ a(i).cid(4:end), b(i).cid(4:end) };
              if ~isequal(q{:}), a(i).cid(4:end)=q{2}; m=m+1; end

              q=[ a(i).nnz, b(i).nnz ];
              if diff(q) && any(q<0), a(i).nnz=q(2); m=m+1; end
           end
           if m, A(k).info.cgr=a; end
        end

        if lflag>1
           a=A(k).info.itags;
           b=B(k).info.itags;
           for i=1:min(numel(a),numel(b))
              if ~isequal(a{i},b{i})
                 b{i}=regexprep(b{i},'\*$','');
                 a{i}=[b{i}, regexprep(a{i},'^[^\*]*','')];
              end
           end
           A(k).info.itags=a;
        end
     end
  end

  ie=builtin('isequal',A,B);
end

