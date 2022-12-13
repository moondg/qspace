function dispopts(varargin)
% Function: dispopts(opts_set [,OPTS])
%
%    display option set in varargin format
%
% Options
%
%   'vid',..   value id, i.e. help which entries represent values (helper set)
%   'istr',..  info string
%
% Wb,Oct23,06

  getopt('init',varargin);
     vid =getopt('vid', []);
     vpat=getopt('vpat',[]);
     fid =getopt('fid', []);
     fpat=getopt('fpat',[]);
     istr=getopt('istr',[]);
  varargin=getopt('get_remaining'); n=length(varargin);

  if ~isempty(istr), fprintf(1,'\n%s\n\n', istr); else inl 1; end
  if n==0, return; elseif n>1
     eval(['help ' mfilename]);
     wberr('invalid usage');
  end

  op=varargin{1}; n=length(op);

  mark=zeros(1,n);

  if ~isempty(vid), mark(vid)=1; mark(vid+1)=2; mark(vid-1)=2; end
  if ~isempty(fid), mark(vid)=3; mark(vid+1)=2; end

  if ~isempty(vpat)
     for i=1:n, if ischar(op{i}) && ~isempty(regexp(op{i},vpat))
     mark(i)=2; mark(i+1)=1; end, end
  end

  if ~isempty(fpat)
     for i=1:n, if ischar(op{i}) && ~isempty(regexp(op{i},fpat))
     mark(i)=3; mark(i+1)=2; end, end
  end

  for i=1:n
     if ~ischar(op{i}) || ~isvarname(op{i}) || numel(op{i})>16
        if i<=1 || none(mark(i-1)==[0 2]), op
        wblog('ERR','expecting name field !??'); end;

        mark(i)=1; if i<n && ~mark(i+1), mark(i+1)=2; end
        mark(i-1)=2;
     elseif ~isempty(regexp(op{i},'-+[a-zA-Z0-9]+')), mark(i)=3;
     end
  end

  i=1; while i<n
     if mark(i), i=i+1; continue; end

     q=1; for k=i+1:n+1
        if k>n || mark(k), break; end
        if mod(k-i,2) && ~isempty(op{k}) && op{k}(1)=='-', q=0; end
     end

     m=k-i; k=k-1;
     if q && mod(m,2)
     mark(i:2:k)=1; mark(i+1:2:k)=2; end

     i=k+1;
  end

  ilast=1; q=1; FF={}; S.FLAGS={};

  i=0; while i<n, i=i+1;
     if mark(i)>1 && ilast<i
        eval(sprintf('S.opts%d=op(%d:%d);',q,ilast,i-1))
        q=q+1;
     end

     if mark(i)==2
        if i<n && mark(i+1)==1
           eval(sprintf('S.%s=op{%d};',op{i},i+1))
           ilast=i+2;
        elseif i==n
           FF{end+1}=op{i}; ilast=i+1;
        end
     end

     if mark(i)==3
        FF{end+1}=op{i}; ilast=i+1;
        continue;
     end
  end

  if ~isempty(FF), S.FLAGS=FF;
  else S=rmfield(S,'FLAGS'); end

  if ilast<i, 
  eval(sprintf('S.opts%d=op(%d:%d);',q,ilast,i)); end

  disp(S)

end

