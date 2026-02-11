function [S,sz]=whos2(varargin)
% Function [S,sz]=whos2(varargin)
%
%    Wrapper routine to Matlab's `whos' routine
%    which accepts further functionality like sort instructions. 
%
% Options
%
%    -s    sort with respect to byte-size
%    -s#   sort with respect to field number # in struct returned by whos
%    -r    same as -s, but reverse sort
%
% Usage #2: whos2(S,...)
%
%    If input S is a plain structure,
%    this generates a `whos' output for all its fields where as with whos,
%    subsequent options may contrain matching field names
%
% Usage #3: whos2(x,field)
%
%    If input x is not a string, then this returns the specified
%    field out of `whos' for given object x. For example, then
%    whos2(x,'class') returns the class name of object x.
%
% Wb,Apr06,08 ; Wb,Jun16,25

% Examples (usage #3): whos2(Q,'class') gets the class name for object Q
% tags: getClassName, classname, get_class_name

  sopt=''; sflag=0; vflag=0; mark=zeros(1,nargin);
  for i=1:nargin
     if ~ischar(varargin{i}), mark(i)=-1; continue
     elseif ~sflag && ~isempty(regexp(varargin{i},'^-[rs]')), mark(i)=1;
        sflag=varargin{i}(2);
        sopt=varargin{i}(3:end);
     elseif ~vflag && isequal(varargin{i},'-v')
        vflag=1; mark(i)=2;
     end
  end

  i=find(mark>0); if ~isempty(i)
     varargin(i)=[]; mark(i)=[];
  end
  narg=numel(varargin); cmd='whos';

  if isequal(mark,[-1 0]) && ~sflag && ~vflag
     q=varargin{1};
     q=whos('q');
     S=getfield(q,varargin{2}); sz=q.bytes;
     return
  elseif narg && isstruct(varargin{1})
     if narg>1, cmd=[ cmd '(' strhcat('-a',varargin(2:end)) ')']; end
     setuser(0,'whos2__',cmd);

     if sflag || nargout
        S=whos_struct(varargin{1});
     else whos_struct(varargin{1}); return
     end
  else
     if narg, cmd=[ cmd '(' strhcat('-a',varargin) ')']; end
     if sflag || nargout
        s=['u=get(0,''UserData''); ' ...
           'u=setfield(u,''whos2__'', ' cmd '); set(0,''UserData'',u)'];
        evalin('caller',s);
        S=getuser(0,'-rm','whos2__');
     else
        evalin('caller',cmd);
        return
     end
  end

  if ~sflag && ~vflag
     if nargout>1, sz=sum([S.bytes]); end
     return
  end

  if isempty(sopt), sopt=3;
  elseif ~isempty(regexp(sopt,'^\d+$')), sopt=str2num(sopt); 
  else wbdie('invalid sort specification'); end

  nS=numel(S);
  if nS>1
     ff=fieldnames(S);
     if numel(sopt)~=1 || sopt>numel(ff)
        wbdie('invalid sort specification'); end
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
        else s=sprintf('x%d',s); s=s(2:end); end

        x=''; if S(i).global,     x=[x 'g']; end
              if S(i).sparse,     x=[x 's']; end
              if S(i).complex,    x=[x 'z']; end
              if S(i).persistent, x=[x 'p']; end
        c=S(i).class;
        if isequal(c,'double'), c=''; end

        fprintf(1,'%-3s%-12s %8s %9s  %s\n',x,S(i).name,...
        s, num2str2(S(i).bytes,'--bytes'), c);
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

