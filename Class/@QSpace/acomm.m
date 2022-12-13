function C=acomm(varargin)
% function C=acomm(A [,B,cflag])
%
%    cflag if set, must equal 'conj'
%    and is interpreted as conjA, i.e. C = {A',B}
%
% Wb,Feb06,16

  C=comm(varargin{:},'+');

end

