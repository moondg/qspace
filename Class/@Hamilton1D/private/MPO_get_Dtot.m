function [Dtot,dloc]=MPO_get_Dtot(HAM)
% function Dtot=MPO_get_Dtot(HAM)
%
%    return total Hilbert space dimension
%    to normalize operator norm2().
%
%    Note however, do not normalize MPOs by default based on
%    Dtot=dloc^(-L), as this quickly leads to confusion
%    when taking powers // DONT_NORM_OPS
%
%    For example, taking the identity operator Id, then Id^p = Id
%    only holds without additional normalization factors!
%    Therefore rather make sure that truncation as with orthoQS()
%    uses rtol *= Dfac [stol *= sqrt(Dfac)] relative to Dfac >> 1.
%
% Wb,Jul26,23

  if ~using_full_MPO(HAM), wbdie('invalid usage (full MPO required)'); end
  if ~isfield(HAM,'info') || ~isfield(HAM.info,'mpo')
     wbdie('invalid usage (full MPO required)');
  end

  dd=MPO_getDim(HAM.mpo);
  dloc=dd(:,3:4,end); d2=diff(dloc,[],2); i=find(d2);
  if ~isempty(i)
     d2=sort(d2(i,:),2);
     wblog('WRN','inconsistent dloc [%d .. %d]',min(d2(:,1)),max(d2(:,2)));
  end

  dloc=max(dloc,[],2);
  Dtot=prod(dloc);

  if Dtot>1E150
   % keep safety margin from double precision overflow at ~1E308
   % e.g. for dloc=4 => Lmax = log(150)/log(dloc) ~ 250
   % or Lmax ~ 500 for dloc=2 (e.g., spin chain)
   % e.g. consider |Ham|^2 ~ |Id|^2 = Dtot, etc.
   % NB! there is no simple overall normalization of MPOs
   % since if say |Ham|^2=tr(H'*H)=1, then |Ham^2|^2 ~ 1/Dtot <<<< 1.
   % generally though, |MPO|^2 < ~ Dtot
   % ==> see also safeguards (ERR/WRN) in this regard
   %     in MPO_trace() and MPO_norm2().
   % Wb,Aug02,23
     wblog('WRN',['MPO operations will likely fail\n' ... 
     'having D_tot = dloc^L = %.3g'],Dtot); 
  end

  if any(diff(dloc))
     wblog('WRN','got various dloc=[%d .. %d]',min(dloc),max(dloc)); 
  else
     dloc=dloc(1);
     if dloc<2, wbdie('got dloc=%d',dloc); end
  end

end

