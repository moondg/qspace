function s=get_sym_info(HAM)
% function s=get_sym_info(HAM)
% Wb,Oct27,21

  IS=HAM.info.IS;

  if isfield(IS,'SOP'), s={};
     for i=1:numel(IS)
        s{i}=strhcat(HAM.info.IS(i).SOP.info,'-s',', ');
     end
     s=strjoin(s,' & ');
  elseif numel(IS)==1 && isfield(IS,'E')
     s=regexprep(IS.E.info.qtype,',',' & ');
  else IS
     wbdie('invalid usage (unexpected content of IS)'); 
  end

end

