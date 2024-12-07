function C=oplus(A,B,varargin)
% function C=oplus(A,B [,d12 [,opts]])
%
%    Direct sum (latex oplus) of two tensors A and B along
%    dimensions d12 in a block-diagonal sense (default: d12=[1 2]).
%    If d12 is single index, oplus is performed along that
%    dimension only. For example, the following are equivalent
%
%                oplus(A,B,2)
%                oplus(A,B,'first')
%                oplus(A,B,[1 2,'first')
% Options
%
%   'first'      direct sum only in index d12(2)
%   'last'       direct sum only in index d12(1)
%                where first/last is interpreted from an MPS point of
%                view, where the first (last) matrix shall be extended
%                along the 2nd (f1st) dimension, respectively.
%                By default: direct sum in (both) indices in d12.
%
%   'bfac',..    apply bfac to B as in: A \oplus (bfac*B)
%
% E.g. see https://mathworld.wolfram.com/MatrixDirectSum.html
% Wb,Nov23,20

% [09/06/2023] permitting single index for d12 without first/last

  getopt('init',varargin);
     if getopt('first'); w=-1;
     elseif getopt('last'); w=1; else w=0; end
     bfac=getopt('bfac',1);
  d12=getopt('get_last',[]); n=numel(d12);

  if numel(A)~=1 || numel(B)~=1
     wbdie('invalid usage (single QSpace tensors A and B required)'); end
  if ~isequal(getqdir(A), getqdir(B))
     wbdie('invalid usage (QSpace qdir mismatch)'); end

  ta=getitags(A); r=numel(ta); if ~r, return; end
  tb=getitags(A);

  for i=1:numel(ta)
      if isempty(ta{i})
         if ~isempty(tb{i}), ta{i}=tb{i}; end
      elseif isempty(tb{i}), tb{i}=ta{i}; end
  end
  if ~isequal(ta,tb), wbdie('invalid usage (QSpace itag mismatch)'); end

  S=struct('type','()','subs',{{}});

  if isempty(d12)
     if ~w, d12=[1 2]; elseif w<0, d12=2; else d12=1; end
  elseif ~isnumeric(d12) || any(d12<1) || n>2 || (w && n~=2)
     wbdie('invalid usage (d12)');
  elseif w<0, d12=d12(2); % 'first' / left boundary
  elseif w>0, d12=d12(1); % 'last'  / right boundary
  end

  for l=d12
     [qa,Ia,Da]=uniquerows(A.Q{l});
     [qb,Ib,Db]=uniquerows(B.Q{l}); [ia,ib,Ix]=matchIndex(qa,qb);

     for i=1:numel(ia)
        sa=size(A.data{Ia{ia(i)}(1)},l);
        sb=size(B.data{Ib{ib(i)}(1)},l);

        for j=Ia{ia(i)}
           a=A.data{j}; s=size(a); s(end+1:l)=1;
           S.subs=matcell(s); S.subs{l}=s(l)+sb; 
           A.data{j}=subsasgn(a,S,0);
        end
        [p,ip]=initperm(r,'--2front',l); px=[];
        for j=Ib{ib(i)}, r_=ndims(B.data{j});
         % bug-fix trailing OM index [email Changkai Zhang 11/26/2024]
           if r_<=r, p_=p; ip_=ip;
           elseif r_==r+1
              if isempty(px), [px,ipx]=initperm(r_,'--2front',l); end
              p_=px; ip_=ipx;
           else wbdie('unexpected rank r_=%d/%d',r_,r); end

           b=permute(B.data{j},p_); s=size(b); s2=prod(s(2:end));
           b=[ zeros(sa,s2); reshape(b,s(1),s2) ]; s(1)=s(1)+sa;
           B.data{j}=permute(reshape(b,s),ip_);
        end
     end
  end

  C=plusQS(A,B,bfac);
  C=setitags(QSpace(C),ta);

end

