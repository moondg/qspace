function [f]=QSet2file(C)
% function [f]=QSet2file(f)
% Wb,Sep20,16

  f=[ C.qdir '/(' QSet2str(C,'-f;') ').cgd'];

end
