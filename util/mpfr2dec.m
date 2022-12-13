%  Usage: [data vector]=mpfr2dec([char matrix] [,base]);
%
%     Convert mpfr data strings to decimals.
%     Since matlab is col-major, the data string for each
%     value is expected to represent an entire column.
%
%  Example
%
%     mpfr2dec(['0':'9', 'A':'Z', 'a':'z'],62)
%
%  Wb,Jul01,14
