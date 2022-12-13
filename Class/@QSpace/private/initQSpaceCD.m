function A=initQSpaceCD(C,perm)
% function A=initQSpaceCD(C [,perm])
%
%    By default, outer multiplicity is placed on last index
%    after permutation, as this is the correct behavior for operators.
%
% Wb,Apr17,17

% C is expected to contain CStore data,
% e.g. as in the output SymStore // LoadCData

  if nargin<2, perm=[]; end
  A=QSpace;

  for k=numel(C):-1:1, Ck=C(k);
     r=numel(Ck.qdir); r_=numel(Ck.cgd.S);
     if r_>r, m=Ck.cgd.S(r+1); else m=1; end

     Ak=get_struct(...
        mat2cell(repmat(Ck.qset,m,1),m,repmat(numel(Ck.qset)/r,1,r)),...
        cell(m,1), struct);
     Ak.info.qtype=regexprep(Ck.type,'[()]','');

     for l=m:-1:1
        d=zeros(1,1,m); d(l)=1; if ~isempty(perm), d=ipermute(d,perm); end
        Ak.data{l}=d;

        c.type=Ak.info.qtype;
        c.qset=Ck.qset;
        c.qdir=Ck.qdir;
        c.cid=Ck.cid;
        c.size=Ck.cgd.S;
        c.cgw=repmat(int8(['0',0]'),1,m); c.cgw(1,l)='1';

        Ak.info.cgr(l,1)=c;
     end

     Ak.info.otype='';

     t=char(Ck.qdir-'+'); t(find(t))='*';
     Ak.info.itags=cellstr(t')';

     A(k)=class(Ak,'QSpace');

     if ~isempty(perm)
        A(k)=QSpace(permuteQS(A(k),perm));
     else
        A(k)=QSpace(permuteQS(A(k),1:numel(Ak.Q)));
     end
  end
end

% -------------------------------------------------------------------- %

