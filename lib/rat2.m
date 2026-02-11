function [sr,Ir]=rat2(x,varargin)
% function str=rat2(x [,fmt])
%
%    Get fractional representation of x, if any.
%    Based on wbrat mex routine.
%
%    The format specifier fmt only affects the string width
%    if x is a matrix (string fmt is enforced).
%
% Options
%
%    '-u',{uval,ustr}  specificy units
%    '-c'              return cell array of strings (rather than single string)
%    '~r'              skip sqrt() check
%    '~a'              do not replace sqrt(..) by ASCII √
%  
% See also wbrat, rats()
% Wb,Mar31,05 ; Wb,Jul04,22

  if ~nargin
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  opts={};

  getopt('init',varargin);
     udat =getopt('-u',{});
     cflag=getopt('-c');
     if getopt('-q'), opts{end+1}='-q'; end

     if getopt('~r'), rflag=0; opts{end+1}='~r'; else rflag=1; end
     if getopt('~a'), aflag=0; else aflag=1; end

  fmt=getopt('get_last','');

  if ~isempty(udat)
     if ~iscell(udat) || numel(udat)~=2 || numel(udata{1})~=1
     wbdie('invalid usage: ''-u'',{uval,ustr}'); end
     if udat{1}>0, x=x/udat{1};
     else
        wblog('WRN','ignorig unit %g',udat{1}); 
        udat={};
     end
  end

  if isempty(x), str=''; e=[]; return; end
  [sr,Ir]=wbrat(x,opts{:}); if ~iscell(sr), sr={sr}; end

  if aflag
     srd=char(8730);
     for i=1:numel(sr)
        sr{i}=regexprep(sr{i},'sqrt\(\s*(\d+)\s*\)',[srd '$1']);
        sr{i}=regexprep(sr{i},'sqrt\((.*)\)',[srd '($1)']);
     end
  end

if cflag, return; end

  if numel(sr)>1
     if isempty(fmt)
        n=4+max(cellfun(@numel,sr),[],'all');
        for i=1:numel(sr), l=length(sr{i});
           if l<n
              s=repmat(floor((n-l)/2),1,2); if mod(l-n,2), s(1)=s(1)+1; end
              sr{i}=[ repmat(' ',1,s(1)), sr{i}, repmat(' ',1,s(2)) ];
           end
        end
     elseif ~isempty(fmt)
       if isempty(regexp(fmt,'%[^a-z]*s'))
          wbdie('invalid usage (string fmt expected)'); 
       end
       for i=1:numel(sr), sr{i}=sprintf(fmt,sr{i}); end
     end
  end

  sr=cell2mat(sr);

end

