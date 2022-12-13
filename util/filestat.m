% Usage: i=filestat(f,tag [,'-s'])
%
%    get system file information
%
%    If '-s' is specified, a time string is returned;
%    otherwise (by default) the returned date is in Matlab datenum
%    format, e.g. can be used with datestr(t,*).
%
% Available tags
%
%   'm'     modification time stamp
%   'c'     creation time stamp
%   'a'     last access time stamp
%
% All these time stamps return split-second resolution.
%
% Wb,Jan28,16
