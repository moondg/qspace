function s=get_sym_info(HAM)
% function s=get_sym_info(HAM)
% Wb,Oct27,21

  s={};
  for i=1:numel(HAM.info.IS)
     s{i}=strhcat(HAM.info.IS(i).SOP.info,'-s',', ');
  end
  s=strjoin(s,' & ');

end

