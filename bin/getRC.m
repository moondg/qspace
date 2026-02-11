% Usage: C=getRC(cgr_in [,opt])
%
%     Get RCStore data for specified input cgr data.
%     If cgr_in is an array, e.g., taken from  QSpace X.info.cgr,
%     this also returns a respective array as output.
%     This routine accesses the RCStore located as specified
%     by the environmental variable RC_STORE.
%
%   Options
%
%     -f  returns CGC data in double format (by default, just CG reference)
%     -F  returns CGC data in QSpace-internal multiprecision format (MPFR)
%
% Usage #2: C=getRC(<sym>,'--info');
%
%     Return listing of currently loaded RC data
%     specifically also including <sym>.
%
% Usage #3: getRC(<sym>,'--ping')
%
%      This just `pings' the RCStore for specified symmetry
%      to ensure that an RCStore for specified symmetry is present.
%      If not, it will be generated.
%
% Wb,Oct 2015 ; Nov 2016 ; Wb,Jun12,25
