function p=permparity(P)
% function p=permparity(P)
%
%    Obtain parity of given permutation P.
%
% Wb,Nov25,20

% https://stackoverflow.com/questions/16351039/..
% is-there-a-built-in-function-in-matlab-to-determine-if-a-permutation-is-even-or

  D=length(P);
  E=speye(D); p=det(E(:,P));

end

