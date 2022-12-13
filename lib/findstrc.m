function [im,nn,ii]=findstrc(ss,q,varargin)
% function [im,nn,ii]=findstrc(ss,q [,opts])
%
%    Find string q in cell string ss based on regexp[i].
%    The eturned data is:
% 
%     im  index to cells in ss that match the input pattern q
%     nn  number of matches for each cell
%     ii  position of first match (unless -f is specified; zero if none)
%
% Options
%
%     -n  return nn as only argument
%     -f  return full cell array in ii for matched positions
%    
%     -i  ignore case in pattern match
%     -x  extend to objects of arbitrary type where only strings are
%         going to be checked for matches
%
%     All remaining options, if any, are forwarded to rexexp[i].
%
% Wb,Aug07,14

% in order to just get the index to the entries that match
% i = find(cellfun(@numel,regexp(ss,pat))); // Wb,Jun09,19

  if nargin<2 || ~iscell(ss)
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  getopt('init',varargin);
     nflag = getopt('-n'); if nflag, fflag=0; else
     fflag = getopt('-f'); end
     iflag = getopt('-i');
     xflag = getopt('-x');
  args=getopt('get_remaining');

  if iflag, 
       my_regexp = @(varargin) regexpi(varargin{:});
  else my_regexp = @(varargin) regexp(varargin{:}); end

  ms=reshape(cellfun(@ischar,ss),1,[]);

  if xflag
     ii=cell(size(ss));
     for i=reshape(find(ms),1,[]), ii{i}=my_regexp(ss{i},q,args{:}); end
  elseif all(ms)
     try ii=my_regexp(ss,q,args{:});
     catch me
        [e,i]=lasterr; e=regexprep(e,char(10),'\\N');
        wberr('invalid input (cellar array of strings required)\N\N%s',e);
     end
  else wberr('invalid usage (char cell array required)'); end

  nn=reshape(cellfun(@numel,ii),size(ss));
  if nflag
     if nargout>1, wberr(...
       'invalid usage (one output argument only with option -c'); end
     im=nn; 
  else
     im=reshape(find(nn),1,[]);
     if nargout>2 && ~fflag
        q=zeros(size(nn)); for j=im, q(j)=ii{j}(1); end
        ii=q;
     end
  end

end

