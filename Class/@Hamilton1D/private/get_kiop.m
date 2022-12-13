function [k,iop]=get_kiop(x)
% function [k,iop]=get_kiop(x)
%    outsourced from dispMPO()
% Wb,Jan19,16

    k=floor(x);
    iop=(x-floor(x))*2^10;

end

