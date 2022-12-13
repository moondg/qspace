function [A,nsp]=split_ops(A,varargin)
% function [A,nsp]=split_ops(A [,opts])
%
%    split composite operators back to (irreducible) operator sets
%    that have well defined irop symmetry labels
%    (returns A -> A(:) expanded along dim2, such that
%    A(k) still references a (partial) original QSpace).
%
% Options
%
%    --set-op    set as operator type
%    --s2f       move scalar operators to front
%    --trans     return transpose array
%    --vec       return as serialzed vector
%    --cell      return as cell array
%
% Wb,Aug31,19

  vectr=0; trans=0; cflag=0;

  getopt('init',varargin);
    oflag=getopt('--set-op');
    sflag=getopt('--s2f');

    if     getopt('--trans'), trans=1;
    elseif getopt('--vec'  ), vectr=1;
    elseif getopt('--cell' ), cflag=1; end

  getopt('check_error');

  sA=size(A); A=reshape(A,[],1); nA=numel(A); dd=repmat(-1,1,nA);
  Ac=cell(1,nA); nsp=ones(1,nA);

  for k=1:nA
     if oflag && isempty(A(k,1).info.otype)
        A(k,1).info.otype='operator';
     end

     Ak=A(k,1); rk=numel(Ak.Q);
     if rk==3
        [q3,i3,d3]=uniquerows(Ak.Q{3}); n3=numel(d3); dd(k)=n3;
        if n3>1, nsp(k)=n3;
           for i=1:n3
               A(k,i)=getsub(Ak,i3{i});
               A(k,i).info.otype=Ak.info.otype;
           end
           if sflag
              iz=find(sum(q3.^2,2)<1E-8);
              if ~isempty(iz), i=1:n3; i(iz)=[];
                 A(k,1:n3)=A(k,[iz,i]);
              end
           end
           if cflag || vectr, Ac{k}=A(k,1:n3); end
        end

     elseif rk~=2, wblog('WRN','got rank-%g QSpace(%g)',rk,k);
     end
  end

  e=norm(diff(dd));
  if e
     wblog('WRN','got composite operators with different symmetry labels');
  end

  if     trans, A=builtin('transpose',A);
  elseif vectr, A=[Ac{:}];
  elseif cflag, A=Ac;
  end

end

