function [c2,qd,qa,qr]=CasimirC2(sym,varargin)
% function [c2,qd,qa,qr]=CasimirC2(sym,varargin)
%
%    Get quadratic Casimir for non-abelian symmetries
%
%                   dim  CasimirC2
%        SU2 1           3/4   = ( 3)/4  \
%        SU3 10          4/3   = ( 8)/6   \
%        SU4 100         15/8  = (15)/8    > consistent with C2=(N^2-1)/2N
%        SU5 1000                         /
%        SU6 10000       35/12 = (35)/12 /
%
%        SU2 2           2               \
%        SU3 11          3                \
%        SU4 101         4                 > consistent with C2=N
%        SU5 1001                         /
%        SU6 10001       6               /
%
%        SU3 20         10/3   = 3.333    consistent with Fuhringer08
%
%                                       |Li|^2 = C2*dim_q/dim_adj
%        Sp4 10     4    5/4   = 1.25   (5/4)*4/10 = 1/2
%        Sp4 01     5    2     = 2.     (  2)*5/10 = 1.
%        Sp4 20    10    3     = 3.     (  3)*10/10= 3.
%
% The quadratic Casimir for SU(3) with q=(p,q) [C1; e.g. see in
% https://en.wikipedia.org/wiki/Clebsch-Gordan_coefficients_for_SU(3)]
%
%    S2 = [ p^2 + q^2 + 3(p+q) + pq) ] / 3
%    e.g. defining or anti-defining: p or q=1 => S2 = 4/3 - ok.
%    e.g. for p=q: S2 = (3p^2 + 6p) / 3 = p(p+2)
%         --> p=q=1 (adjoint): S2=3=N - ok.
%
% See also PHYS//Casimir
% Wb,Oct05,15

% SymStore('Sp4','Casimir','10')

  qr=LoadRData(sym,varargin{:});
  [qd,qa]=get_defining(sym);

  dr=size(qr.Z,1); r=cgs2double(qr.Sp(1));
  dd=size(qd.Z,1); d=cgs2double(qd.Sp(1)); 
  da=size(qa.Z,1);

  c2 = norm(r(:))^2 / norm(d(:))^2;

  c2 = 0.5 * (da/dr) * c2;

end

