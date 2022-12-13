function [ff,I]=getnrgfiles(nrg,varargin)
% function [ff,I]=getnrgfiles([nrg])
%
%    get filenames of NRG file structure (default nrg: './NRG/NRG')
%    if requested, I returns data of info file.
%
% Options
%
%   -f   full name, i.e.,
%        prepend every file name with path as returned by dir()
%
% Wb,Aug22,12

  getopt('init',varargin);
     fflag=getopt('-f');
  nrg=getopt('get_last','./NRG/NRG');

  if ~ischar(nrg), nrg, error('Wb:ERR',...
     'invalid NRG data specification (string expected)'); end

  ff=dir2([nrg '_\d+.mat']);

  if nargout>1
     i=[ ff(1).folder '/' regexprep(ff(1).name,'_\d+.mat','') '_info.mat'];
     if exist(i,'file'), I=load(i);
     else wberr('failed to find info file (%s)',nrg); end
  end

  if fflag
     for i=1:numel(ff), p=ff(i).folder;
        if ~isempty(p) && ~isequal(p,'.')
            ff(i).name=[ p '/' ff(i).name ];
        end
     end
  end

end

