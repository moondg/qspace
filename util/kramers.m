%  Usage: [yr,I]=kramers(xi,yi [,OPTS])
%
%    calculate Kramers-Kronig transform from discrete data
%    using linear interpolation scheme.
%
%      xi         x-data
%      yi         y-data
%      yr         kramers kronig transformed data (default: imag->real)
%
%  Options
%
%     '-disc'     yi data is already discrete
%     'eps',..    for continuous data: regularization of principal
%                 value data for xi==xj (xi += eps*(xi-x[i-1]; 1E-10)
%   If isdisc
%
%     'faceta',.. width of Lorentzian delta function relative to
%                 data point spacing dx (1.)
%  else
%
%     'g1flag'    include tails to infinity by assuming that
%                 input yi decays like 1/x for large x
%     'g2flag'    same as g1flag but assuming 1/x^2 decay
%     'r2i'       do real->imag transform instead (=overall minus sign)
%
%  Wb,Dec12,06
