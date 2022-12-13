function [S,e,istr]=addfields(S,Sref)
% function [S,e,istr]=addfields(S,Sref)
%
%    adds extra fields in Sref to S.
%    if the final structure has exactly the same fields,
%    the fields in S are (re)ordered exactly the same way as in Sref.
%
% Return values
% 
%    e=0   got structure of exactly the same type.
%    e=nan got the same number of fields, yet these needed to be reordered.
%    e>0   got (+e) extra fields in S.
%    e<0   got (-e) extra fields in Sref, with no extra fields in S itself
%          => output S will have the same order of fields as Sref.
%
%    istr  string that describes difference in between input structures
%          (empty if none)
%
% Wb,Apr19,14

  f0=fieldnames(Sref);  n0=numel(f0);
  f1=fieldnames(S);     n1=numel(f1);

  if isequal(f0,f1), e=0; istr=''; return; end

  fx=setdiff(f0,f1)'; nx=numel(fx);
  for f=fx
     S=setfield(S,f{1},getfield(Sref,f{1}));
  end

  n2=n1+nx; e=n2-n0;
  if e==0
     S=orderfields(S,Sref); e=-nx;
     if nargout>2
        istr=sprintf('added fields {%s}',strhcat(fx{:},'-s',', '));
     end
  elseif nargout>2
     f2=setdiff(f1,f0);
     istr=sprintf('significant difference in fields {extra %s; added %s}',...
     strhcat(f2{:},'-s',', '), strhcat(fx{:},'-s',', '));
  end

end

