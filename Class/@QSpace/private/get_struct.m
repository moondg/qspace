function A=get_struct(Q,data,info)
% function A=get_struct(Q,data,info)
% Wb,Apr17,17

   if nargin<3, info=[];
      if nargin<2, data=cell(0,1);
         if nargin<1, Q=cell(1,0); end
      end
   end

   if ~isempty(info) && ~isstruct(info)
      n=[ numel(Q), numel(info) ];
      if iscell(info) && iscell(Q) && n(1)>=2 && (n(2)==n(1) || n(2)==n(1)+1)
         info=struct('qtype',{'A'},'otype',{''},'itags',{info},'cgr',[]);
         if diff(n)
            Q(end+1:n(2))={zeros(size(Q{1}))};
         end
      else
         wblog('ERR','unexpected info data');
      end
   end

   A=struct('Q',{Q}, 'data',{data},'info',{info});

end

