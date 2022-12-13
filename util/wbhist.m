%  Usage: ybin=wbhist([strid,] x,y [,xbin]);
%
%     mex routine to collect binned data. If y is not specified,
%     this generates a simple histogram counting occurances.
%
%     NB! The binning occurs with respect to the mid-points of
%     the xbin data, where the first bin contains everything for
%     x < mean(xbin(1:2)), and the last bin everything for
%     x >= mean(xbin(end-1:end)).
%
%     NB! data is stored globally within mex file, i.e. is
%     cleared when using 'clear functions'. However, this allows
%     subsequent calls without last argument xbin to add to
%     existing binned data. The last argument must be specified
%     for the first call given some strid, and whenever the
%     histogram data shall be initialized.
%
%     default hist set has strid 'default'.
%
%  Further tags:
%
%    'stat'         current status of wbhist (global)
%    'clear'        clear entire wbhist (global)
%    'id,'remove'   remove wbhist data set associated with 'id'
%
%  Examples
%
%     wbhist([],[],xbin);
%     wbhist(x, []);   plain count (assumes y=ones(size(x))).
%     [a,b,c]=wbhist(['default']);
%
%  Wb,Feb12,11
