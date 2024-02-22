function varargout=subsref(A,S)
% function Aout=subsref(A,S)
%
%   This manages access to the fields
%   in the underlying structure of QSpace A.
%
% AW (2012)

% NB! This routine is required to access fields within a QSpace tensor
% outside this class environment. If not defined, matlab issues error:
% "Access to an object's fields is only permitted within its methods."

% WRN! subsref() behaves like a regular matlab function. This behaves
% different e.g., from A.data{:} which returns a `comma separated list'.
% Therefore if multiple output arguments are requested, as in A.Q{:}
% or A.data{I}, then *all* output arguments need to be assigned:
%
%    q=cell(size(I)); [q{:}]=A.data{I};   % ok
%
%    q={A.data{:}}  % likely not ok: only returns {A.data{1}} = A.data(1);
%    Q=[A.Q{:}]     % similar problematic syntax, hence not supported
% 
% Wb,Nov09,23

  if numel(A)==1 && numel(S)==2 && isequal([S.type],'.{}')
     A=builtin('subsref',A,S(1)); S_2=S(2); S_2.type='()';
     S_2.subs={[S_2.subs{:}]};
     varargout=subsref(A,S_2);
  else
     varargout={builtin('subsref',A,S)};
  end

end

