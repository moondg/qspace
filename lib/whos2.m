function [S,sz]=whos2(varargin)
% Function [S,sz]=whos2(varargin)
%
%    To be used exactly the same way as whos but accepts further
%    sort instructions. Also, input may be a simple structure,
%    which then generates a `whos' for all its fields.
%
%    -s    sort with respect to byte-size
%    -s#   sort with respect to field number # in struct whos
%    -r    same as -s, but reverse sort
%
% Wb,Apr06,08

  sopt=''; sflag=0; vflag=0; mark=zeros(1,nargin);
  for i=1:nargin
     if ~ischar(varargin{i}), continue
     elseif ~sflag && ~isempty(regexp(varargin{i},'^-[rs]'))
        sflag=varargin{i}(2); sopt=varargin{i}(3:end); mark(i)=1;
     elseif ~vflag && isequal(varargin{i},'-v')
        vflag=1; mark(i)=2;
     end
  end
  varargin(find(mark))=[];
  narg=numel(varargin); w='whos';

  if narg && isstruct(varargin{1})
     if narg>1, w=[ w '(' strhcat('-a',varargin(2:end)) ')']; end
     setuser(0,'whos2__',w);

     if sflag || nargout
        S=whos_struct(varargin{1});
     else whos_struct(varargin{1}); return
     end
  else
     if narg, w=[ w '(' strhcat('-a',varargin) ')']; end
     if sflag || nargout
        s=['u=get(0,''UserData''); ' ...
           'u=setfield(u,''whos2__'', ' w '); set(0,''UserData'',u)'];
        evalin('caller',s);
        S=getuser(0,'-rm','whos2__');
     else
        evalin('caller',w);
        return
     end
  end

  if ~sflag && ~vflag
     if nargout>1, sz=sum([S.bytes]); end
     return
  end

  if isempty(sopt), sopt=3;
  elseif ~isempty(regexp(sopt,'^\d+$')), sopt=str2num(sopt); 
  else wberr('invalid sort specification'); end

  nS=numel(S);
  if nS>1
     ff=fieldnames(S);
     if numel(sopt)~=1 || sopt>numel(ff)
        wberr('invalid sort specification'); end
     fs=ff{sopt}; ss=cell(nS,1);
     for i=1:nS, ss{i}=getfield(S(i),fs); end
     if isnumeric(ss{1})
        for i=1:nS, ss{i}=prod(ss{i}); end
        ss=[ss{:}];
     end
     [sx,is]=sort(ss);

     if isequal(sflag,'r'), is=flip(is); end
     S=S(is);
  end

  if ~nargout || vflag
     fprintf(1,'\n   Name             size     bytes  class\n\n');
     for i=1:nS
        s=S(i).size; if isempty(s), s='[]'; 
        elseif prod(s)==1, s='1';
        else s=sprintf('x%g',s); s=s(2:end); end

        x=''; if S(i).global,     x=[x 'g']; end
              if S(i).sparse,     x=[x 's']; end
              if S(i).complex,    x=[x 'z']; end
              if S(i).persistent, x=[x 'p']; end
        c=S(i).class;
        if isequal(c,'double'), c=''; end

        fprintf(1,'%-3s%-12s %8s %9d  %s\n',x,S(i).name, s, S(i).bytes,c);
     end
     fprintf(1,'\n'); if ~vflag, clear S, end
  end

end

% -------------------------------------------------------------------- %

function S=whos_struct(S__)
   structexp(S__); clear S__
   if nargout
      eval(['S=' getuser(0,'-rm','whos2__') ';']);
   else
      eval(getuser(0,'-rm','whos2__'));
   end
end

% -------------------------------------------------------------------- %

