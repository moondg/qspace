function E=get_local_id(HAM,k,varargin)
% function E=get_local_id(HAM,k [,itag])
% function E=get_local_id(HAM,k [A,j])
%
%    Get identity operator for local state space at site k,
%    assigning itag if specified.
%
% Wb,May19,17

  i=get_stype(HAM,k);
  E=HAM.oez(1,i).op;

  if nargin<3, return
  elseif nargin==3, t=varargin{1};
  elseif nargin>4, wbdie('invalid usage');
  else
     A=varargin{1}; t=A.info.itags;
     j=varargin{2}; if j<=numel(t), t=t{j};
     else wbdie('index j=%g/g out of bounds',j,numel(t)); end

     if ~isempty(A.Q)
        [i1,i2,Im]=matchIndex(A.Q{j},E.Q{1});
        if ~isempty(Im.ix1)
           wbdie('got local state space mismatch !?'); 
        end
     end
  end

  t=regexprep(t,'\**$',''); % skip trailing conj flags '*'
  E.info.itags=regexprep(E.info.itags,'^\w*',t,'emptymatch');

end

