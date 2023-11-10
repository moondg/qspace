function i=got_gauge(HAM)
% function i=got_gauge(HAM)
%
%    return zero if HAM.info.param.gauge is not set;
%    otherwise
%      2 if gtype contains string '..op'
%      1 else
%
% Wb,Apr30,23

  i=0;

  if isfield(HAM.info.param,'gauge'), Ig=HAM.info.param.gauge;
     if numel(Ig)==1 && isfield(Ig,'gtype') && ...
        isfield(Ig,'iq') && isfield(Ig,'g2')
        if isempty(regexp(Ig.gtype,'op$')), i=1; else i=2; end

        iq=Ig.iq; g2=Ig.g2; r=size(HAM.oez(1).op.Q{1},2);
        if iq<1 || iq>r, wbdie('invalid param.gauge.iq = %d / %d',iq,r); end
        if numel(g2)~=1 || ~isreal(g2) || g2<0, disp(g2)
           wbdie('invalid param.gauge.g2');
        end
     else wblog('WRN','got unexpected gauge setting - ignore'); 
     end
  end

end

