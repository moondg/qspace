function [ss,s]=getsym(A,varargin)
% function [ss,s]=getsym(A [,dim,opts])
% Options
%
%   -c        return symmetries split as cell array
%   -d        get dimensions (rank) of each symmetry [returns 1 for Abelian]
%   -r        get symmetry-rank for each symmetry [returns 0 for Abelian]
%   'ir',dim  index range in q-labels for given dimension
%   '-I',dim  get extended info structure for given symmetry
%
% Wb,Oct24,11

  iflag=0;
  getopt('INIT',varargin);
     vflag=getopt('-v');
     cflag=getopt('-c');

     if     getopt('-d'), dflag=1;
     elseif getopt('-r'), dflag=2; else dflag=0; end

     ir=getopt('ir',[]); if isempty(ir)
     ir=getopt('-I',[]); if ~isempty(ir), iflag=1; end; end
  dim=getopt('get_last',[]);

  if ~isempty(dim)
     if ~isnumber(dim) || cflag
        wblog('ERR','%s() invalid usage',mfilename);
        return
     end
     cflag=1;
  end
  if dflag || ~isempty(ir), cflag=1; if ~dflag, dflag=1; end; end

  if dflag, ss=[]; elseif cflag, ss={}; else ss=''; end
  if isempty(A), return; end

  for k=1:numel(A), Ak=struct(A(k));
     if isfield(Ak.info,'qtype') && ~isempty(Ak.info.qtype), s=Ak.info.qtype;
        q=strread(s,'%s','delimiter',',');
        if numel(q)>size(Ak.Q{1},2), wbdie(...
          'inconsistent symmetry type (''%s''; got rank-%d)',s,size(Ak.Q{1},2));
        end
        if cflag, s=q; end
     else
        if ~isempty(Ak.Q), ns=size(Ak.Q{1},2);
           if cflag, s=repmat({'A'},1,ns);
           else s=repmat(',A',1,ns); s=s(2:end); end
        else s={}; end
     end
     if k>1, if ~isequal(s,ss)
        wbdie('inconsistent symmetry type (%s <> %s)',s,ss); end
     else ss=s; end
  end

  if dflag
     for i=1:numel(ss), s=ss{i};
        if regexp(s,'^([AP]|Z\d+)$') % 'A','P','Z2', etc.
           if dflag>1, ds(i)=0; else ds(i)=1; end
        elseif regexp(s,'^SU\d+$'), ds(i)=str2num(s(3:end))-1;
        elseif regexp(s,'^Sp\d+$'), ds(i)=str2num(s(3:end))/2;
        elseif regexp(s,'^SO\d+$')
           ds(i)=floor(str2num(s(3:end))/2);
        else 
           wbdie('unknown symmetry `%s''',s);
        end

        if ds(i)~=round(ds(i))
           wbdie('invalid symmetry rank %g (%s)',ds(i),s);
        end
     end

     if ~isempty(ir)
        if numel(ir)>1, wbdie('SINGLE index expected for ''ir''');
        elseif ~isempty(dim), wbdie(...
          'dimension shall be specified with options `ir'' itself');
        end

        s=ss{ir};
        ic=cumsum([1 ds]); ss=ic(ir):(ic(ir+1)-1);

        if iflag
           ss=struct('sym',s,'is',ir,'j',ss,'r',numel(ss));
           clear s
        end
     else
        if isempty(dim), s=ss; ss=ds; else s=ss{dim}; ss=ds(dim); end
     end
  else
     if ~isempty(dim), ss=ss{dim}; end
     if vflag
        if iscell(ss)
             for i=1:numel(ss), ss{i}=verbose_sym(ss{i}); end
        else ss=verbose_sym(ss);
        end
     end
     clear s
  end

end

% -------------------------------------------------------------------- %

function s=verbose_sym(s)
   if regexp(s,'^S[pU](\d+)$'), s=[s(1:2) '(' s(3:end) ')']; end
end

% -------------------------------------------------------------------- %

