function str = basename(fname)
% str = basename(fname)
%
%       returns fname up to the first "."
%
% Wb,Nov29,00

    [p,f,e]=fileparts(fname); str=[f,e];

end

