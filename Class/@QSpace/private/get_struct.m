function A=get_struct(Q,data,info)
% function A=get_struct(Q,data,info)
%
%    Ensure presence of required fields in correct order.
%    for QSpace object. With no arguments, this returns
%    empty QSpace structure.
%
% Wb,Apr17,17

% tags: get_empty_QSpace_struct, get_empty_struct_QSpace

  if nargin==3
     if ~isempty(info) && ~isstruct(info)
        n=[ numel(Q), numel(info) ]; q3=diff(n);
        if ~iscell(info) || ~iscell(Q)  || q3<0 || q3>1
           wbdie('unexpected info data'); end
        info=struct('qtype',{'A'},'otype',{''},'itags',{info},'fdir',{''},'cgr',[]);
        if q3
           Q(end+1:n(2))={zeros(size(Q{1}))};
        end
     end
  else, info=[];
    if nargin<2, data=cell(0,1);
       if nargin<1, Q=cell(1,0); end
    end
  end

  A=struct('Q',{Q}, 'data',{data},'info',{info});

end

