function A=permute(A,varargin)
% Usage: A=permute(A)
% permute QSpace with given permutation
% Wb,Sep08,06 ; Wb,Nov24,14

 % NB! need to be careful about CGRef data!
 % e.g. if the same q-label appears multiple times
 % cgp+conj is not necessarily the same as obtained through
 % QSpace::CGRef::Sort(cgp,conj)! // Wb,Dec04,14
 % => use permuteQS() instead of modifying info.cgr data

   for k=1:numel(A)
      A(k)=QSpace(permuteQS(A(k),varargin{:}));
   end

end

