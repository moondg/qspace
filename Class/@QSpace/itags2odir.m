function odir=itags2odir(A,varargin)
% function odir=itags2odir(A [,opts])
%
%    find index that points outwards of A (outgoing index;
%    expecting single outgoing index, by default).
%
% Options: see usage for private/itags_to_odir.m
% Wb,May15,17

   odir=cell(size(A)); nA=numel(A);
   for k=1:nA
      odir{k}=itags_to_odir(A.info.itags,varargin{:});
   end

   if nA==1, odir=odir{1};
   elseif nargin==1, odir=reshape([odir{:}],size(A));
   end

end

