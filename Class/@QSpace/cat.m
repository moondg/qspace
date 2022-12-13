function X=cat(varargin)
% function X=cat(A1,A2,...,k)
%
%    Concatenate data sectors of given QSpaces along dimension k
%    (QSpace vectors are also included in the concatenation procedure).
%    NB! k may be a cell, i.e. {k,qdir} to specify qdir of newly added
%    dimensions.
%
% Wb,Sep17,18

  X=QSpace; k=varargin{end}; varargin(end)=[];
  for i=1:numel(varargin)
  for j=1:numel(varargin{i})
      X=cat_1(X,varargin{i}(j),k);
  end
  end

end

% -------------------------------------------------------------------- %

function X=cat_1(A,B,k)

   if isempty(A.data), X=B; return; end
   if isempty(B.data), X=A; return; end

 % by using X=A+B below, this should also take care of outer multiplicity!
 % if checkOM(A), error('Wb:ERR',...
 %    '\n   ERR outer multiplicity not yet implemented in cat()'); end
 % QA=[A.Q{:}]; QB=[B.Q{:}];
 % if size(uniquerows(QA),1)~=size(QA,1) || size(uniquerows(QB),1)~=size(QB,1)
 %    error('Wb:ERR','\n   ERR outer multiplicity not yet implemented in cat()');
 % end

   ra=numel(A.Q); rb=numel(B.Q);

   if iscell(k)
      if numel(k)~=2 || ~isnumber(k{1}) || ~isnumeric(k{2})
         error('Wb:ERR','\n   ERR invalid usage'); end
      qdir=k{2}; k=k{1};
   elseif ra<rb, qdir=getqdir(B); else qdir=getqdir(A);
   end
   qdir(end+1:k)=qdir(end);

   if k>ra, A=appendSingletons(A,qdir); end
   if k>rb, B=appendSingletons(B,qdir); end

   [qa,da]=getQDimQS(A,k); [~,Ia]=matchIndex(A.Q{k},qa,'-s');
   [qb,db]=getQDimQS(B,k); [~,Ib]=matchIndex(B.Q{k},qb,'-s');

   [ia,ib,Im]=matchIndex(qa,qb);
   n=size(qa,1); Jb=zeros(1,n); Jb(ia)=ib;
   n=size(qb,1); Ja=zeros(1,n); Ja(ib)=ia;

   n=numel(A.data);
   for i=1:n, j=Jb(Ia(i));
      if j
         I=size(A.data{i}); I(end+1:k)=1; I(k)=I(k)+db(j);
         S=struct('type','()','subs',{mat2cell(I,1,ones(1,numel(I)))});
         A.data{i}=subsasgn(A.data{i},S,0);
      end
   end

   n=numel(B.data);
   for i=1:n, j=Ja(Ib(i));
      if j
         s=size(B.data{i}); s(end+1:k)=1; s(k)=da(j);
         B.data{i}=cat(k,zeros(s),B.data{i});
      end
   end

   X=A+B;

end

% -------------------------------------------------------------------- %

