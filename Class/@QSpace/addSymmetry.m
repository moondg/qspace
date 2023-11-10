function A=addSymmetry(A,qtype,varargin)
% function A=addSymmetry(A,qtype [,opts])
%
%    Append scalar symmetry (e.g. new symmetry label,
%    where all existing states get a newly added qlabel q=0
%    with the zeros repeated according to the rank of the
%    added symmetry.
%
% Options
%
%   'pos',..  actually insert symmetry at particular position
%             in symmetry order; e.g., 'pos',1 will prepend the symmetry.
%   'q',..    specify particular symmetry labels other than (default q=0).
%             if non-zero, this assumes A to be of operator type in that
%             specified q is applied to the first two indices only *)
%
% *) rather than specifying non-zero q, consider starting with zeros
%    first, followed by mapping of operators into tensor product
%    of state spaces of operators
%
% See also: makeIrop.m
% Wb,Jul12,11 ; Wb,Dec08,14

% tags: extendQ addQ expandQ, prependScalarSymmetry
% removed q from original usage [A=addSymmetry(A,q,qtype)]
% and added pos as optional argument // Wb,Feb06,14

% renamed from appendScalarSymmetry -> addSymmetry
% since name was rather misleading by now // Wb,Nov04,23

  getopt('init',varargin);
     qs =getopt('q',[]);
     pos=getopt('pos',0);
  getopt('check_error');

  if isempty(qs)
     if isequal(qtype,'A') || ~isempty(regexp(qtype,'^Z\d+$')), r=1;
     elseif ~isempty(regexp(qtype,'SU\d+$'))
        r=str2num(qtype(end))-1;
     elseif ~isempty(regexp(qtype,'Sp\d+$'))
        r=str2num(qtype(end))/2;
     elseif ~isempty(regexp(qtype,'Z\d+$')), r=1;
     else
        wbdie('invalid/unknown symmetry (%s)',qtype);
     end
     if r<1 || r>9 || r~=round(r)
        wbdie('invalid usage (got rank=%g !??)',r);
     end
     qs=zeros(1,r);
  end

  for k=1:numel(A), Ak=A(k); Ik=Ak.info; rA=numel(Ak.Q);
     qd_=getqdir(Ak,'-s');

     if ~isempty(Ik.cgr)
        if ~isstruct(Ik.cgr), wbdie('invalid info.cgr field'); end
        Ik.qtype=[Ik.qtype ',' qtype];

        qd=Ik.cgr(1).qdir;
        if isempty(qd), qd=qd_; end

      % by keeping remaining fields (cid, size, cgt, ...) empty,
      % this either copies this information from the existing
      % corresponding CData in gCG.BUF if present, otherwise the
      % corresponding CData is loaded from the RCStore
      % see also QSpace/getvac // INIT_SCALAR // Wb,Jan29,15
        c=emptystruct(Ik.cgr);
          c.type=qtype;
          c.qset=get_qset(qs,rA);
          c.qdir=qd;
          c.cgw=1;

        Ik.cgr(:,end+1)=c;
        Ak.info=Ik;

     elseif isequal(qtype,'A')
        if ~isempty(Ik.qtype)
           Ak.info.qtype=[Ak.info.qtype ',' qtype];
        end
     else
        if isempty(Ik.qtype)
             Ik.qtype=[repmat('A,',1,size(Ak.Q{1},2)), qtype];
        else Ik.qtype=[Ik.qtype ',' qtype];
        end
        qt=strread(Ik.qtype,'%s','delimiter',',;');
        nq=numel(qt);

        c=repmat(getRC('--empty'),1,nq);
          c(end).type=qtype;
          c(end).qset=get_qset(qs,rA);
          c(end).qdir=qd_;
          c(end).cgw=1;
        Ik.cgr=repmat(c,size(Ak.Q{1},1),1);
        Ak.info=Ik;

     end

     Q=Ak.Q; rA=numel(Q); qi=qs;
     for i=1:rA
        if i==3 && any(qi), qi(:)=0; end
        Q{i}=[Q{i},repmat(qi,size(Q{i},1),1)];
     end
     Ak.Q=Q;

     if norm(qs) && ~isempty(Ik.cgr)
        Ak=QSpace(permuteQS(Ak,'12'));
        s=Ak.info.cgr(1,end).size;
        if s(1)>1, wfac=sqrt(double(s(1)));
           if rA<=3
              for i=1:size(Ik.cgr,1)
                  Ak.info.cgr(i,end).cgw=wfac;
              end
           else
              for i=1:numel(Ak.data),
                  Ak.data{i}=wfac*Ak.data{i};
              end
           end
        end
     end

     if pos>0, n=numsym(Ak);
        if pos>n, wbdie('index pos out of range (%g/%g) !?',pos,n);
        elseif pos<n
           perm=[1:pos-1,n,pos:(n-1)];
           Ak=symperm(Ak,perm);
        end
     end

     A(k)=Ak;
  end
end

% -------------------------------------------------------------------- %
function q=get_qset(qs,rA)
   if rA>2 && any(qs)
        q=[ repmat(qs,1,2), repmat(zeros(size(qs)),1,rA-2) ];
   else q=repmat(qs,1,rA);
   end
end

% -------------------------------------------------------------------- %

