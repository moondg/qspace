function [ll, name] = lineno(varargin)
% function [ll, name] = lineno([id=1, opts])
%
%    get line-number and name of calling function
%    returns line number and function name of calling function,
%    (0 and '' if called from the Matlab command-line prompt.)
%
%    id sets the offset:
%    id=1 = current mfile, id=2 = calling mfile, etc.
%
% Options
%
%   'nx',..  number of highest stack entries to skip
%    -s      short (do not print line numbers)
%    -n      number only
%   'all'    display full stack inline
%   'ALL'    display full stack inline (including subroutine names)
%   'len',.. max. length of output string (for nargout==1 only)
%   '<FMT>'  explicit format specification
%            '%.20F -> %S:%L'  (F=file, S=subroutine, L=line)
%
% Output - if nargout
%
%   =0: dbstack is displayed only
%   =1: returned object is single string "name:line"
%       an optional length argument shortens the string using shortfstr()
%
% Ref. http://woodshole.er.usgs.gov/...
% staffpages/cdenham/public_html/snackbar/lineno.m
%
% Wb,Nov13,05  Wb,Jun16,07

  if nargin && isnumeric(varargin{1})
  id =varargin{1}; varargin(1)=[]; else id=1; end

  getopt('INIT',varargin);
     sflag=getopt('-s');
     nflag=getopt('-n');
     len  =getopt('len',-1);
     nx   =getopt('nx',0)+2;

     if getopt('ALL'), aflag=2; 
     elseif getopt('all'); aflag=1; else aflag=0; end

  varargin=getopt('get_remaining');

  if length(varargin)+aflag>1, wbdie('invalid usage'); end

  if length(varargin), fmt=varargin{1}; else fmt=[]; end

  [stack, index] = dbstack;

  if aflag
     if aflag==2
        ll={stack.file; stack.name; stack.line};
        ll=fliplr(ll(:,nx:end));
        if sflag, ll=ll(1:2,:);
             ll=sprintf('%s::%s -> ',ll{:});
        else ll=sprintf('%s::%s:%d -> ',ll{:});
        end; ll=ll(1:end-4);
     elseif aflag==1
        ll={stack.file; stack.line}; ll=fliplr(ll(:,nx:end));
        if sflag, ll=ll(1,:);
             ll=sprintf('%s -> ',ll{:});
        else ll=sprintf('%s:%d -> ',ll{:});
        end; ll=ll(1:end-4);
     end
     return
  end

  if numel(stack)>1
     id = max([min([1+id, length(stack)]), 1]);
     stack=stack(id);
  else
     stack=struct('line',0, 'file','', 'name','');
  end

  if ~isempty(fmt)
     [i1,i2,m]=regexp(fmt,'%[^FSL]*[FSL]','start','end','match'); args={};
     for l=1:length(i1)
        switch fmt(i2(l))
           case 'F', fmt(i2(l))='s'; args{end+1}=stack.file;
           case 'S', fmt(i2(l))='s'; args{end+1}=stack.name;
           case 'L', fmt(i2(l))='d'; args{end+1}=stack.line;
           otherwise wbdie('invalid fmt=`%s''',m{l});
        end
     end
     ll=sprintf(fmt,args{:});
     return
  end

  if ~nargout && ~sflag && ~nflag, disp(stack)
  elseif nargout<=1 && ~nflag
     if ~isequal(stack.file, [stack.name '.m']);
          ll = sprintf('%s>%s:%d', stack.file, stack.name, stack.line);
     else ll = sprintf('%s:%d', stack.file, stack.line); end
     if len>1, ll=shortfstr(ll,len); end
  else
     ll=stack.line;
     name=stack.name;
  end

end

