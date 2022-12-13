function C=oplus(A,B,varargin)
% function C=oplus(A,B [r12,opts]
%
%    Direct sum (latex oplus) of two tensors A and B along dimension r.
%    If r12 is not specified, an MPS or MPO like setting is assumed
%    with L/R as first and second index, i.e., r12=[1 2].
%
% Options
%
%   'first'             direct sum only in index r12(2)
%   'last'              direct sum only in index r12(1)
%                       otherwise: direct sum in both indices r12
%
%   'bfac',..           apply bfac to B prior to extending A
%
% E.g. see https://mathworld.wolfram.com/MatrixDirectSum.html
% Wb,Nov23,20

  getopt('init',varargin);
     if getopt('first'); w=-1;
     elseif getopt('last'); w=1; else w=0; end
     bfac=getopt('bfac',1);
  r12=getopt('get_last',[]);

  if isempty(r12), r12=[1 2];
  elseif ~isnumeric(r12) || numel(r12)~=2 || any(r12<1), r12
     wberr('invalid usage (r12)');
  end

  if numel(A)~=1 || numel(B)~=1
     wberr('invalid usage (single QSpace tensors A and B required)'); end
  q={ getqdir(A), getqdir(B) }; if ~isequal(q{:})
     wberr('invalid usage (QSpace qdir mismatch)'); end

  ta=getitags(A); r=numel(ta);
  tb=getitags(A); r(2)=numel(tb);
  if diff(r), wberr('invalid usage (QSpace rang mismatch %g/%g)',r); end

  ra=r(1); if ~ra, return; end

  for i=1:numel(ta)
      if isempty(ta{i}), if ~isempty(tb{i}), ta{i}=tb{i}; end
      elseif isempty(tb{i}), tb{i}=ta{i}; end
  end
  if ~isequal(ta,tb), wberr('invalid usage (QSpace itag mismatch)'); end

  if     w<0, r12=r12(2);
  elseif w>0, r12=r12(1);
  end

  S=struct('type','()','subs',{{}});

  for r=r12
     [qa,Ia,Da]=uniquerows(A.Q{r});
     [qb,Ib,Db]=uniquerows(B.Q{r}); [ia,ib,Ix]=matchIndex(qa,qb);

     for i=1:numel(ia)
        sa=size(A.data{Ia{ia(i)}(1)},r);
        sb=size(B.data{Ib{ib(i)}(1)},r);

        for j=Ia{ia(i)}
           a=A.data{j}; s=size(a); s(end+1:r)=1;
           S.subs=matcell(s); S.subs{r}=s(r)+sb; 
           A.data{j}=subsasgn(a,S,0);
        end
        [p,ip]=initperm(ra,'--2front',r);
        for j=Ib{ib(i)}
           b=permute(B.data{j},p); s=size(b); b=reshape(b,s(1),[]);
           b=[ zeros(sa,size(b,2)); b ]; s(1)=s(1)+sa;
           B.data{j}=permute(reshape(b,s),ip);
        end
     end
  end

  C=plusQS(A,B,bfac);
  C=setitags(QSpace(C),ta);

end

