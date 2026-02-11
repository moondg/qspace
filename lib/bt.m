function bt(varargin)
% Function bt([n,opts])
%    
%    Show backtrace of current caller stack, e.g., when debugging
%    skipping lowest n entries (default: 0).
%
% Options
%
%    -v|--show [m]  show source of last m stack entries (default: m=3)
%
% Wb,Aug27,08 ; Wb,Jul16,23

  persistent gStack iStack

  n=0; iS=[]; nargs=0; df=[]; vflag=1;
  if nargin
     if isnumber(varargin{1}), n=varargin{1}; varargin(1)=[]; end
     nargs=numel(varargin);
     if nargs, q=varargin{1};
        if ~isempty(regexpi(q,'^(--show|-v$)'))
           vflag=2; if q(end)=='V', vflag=3; end
           iS=varargin(2:end);
        elseif ~isempty(regexp(varargin{1},'^--(up|down)(?@q=$1;)')) && nargs<=2
           df=[varargin{2:end}];
           if isempty(df), df=1; elseif ischar(df), df=str2num(df); end
           if isequal(q,'down'), df=-df; end
        else disp(varargin), wbdie('invalid usage'); end
     end
  end

if isempty(gStack) || isempty(df)
  [S,i]=dbstack('-completenames');

  if 2+n>1
     S=S(2+n:end); if isempty(S), return; end
  end

% -------------------------------------------------------------------- %
% get current frame
% => useless, since for the call to dbstack here will always
%    return *this as current frame!
% -------------------------------------------------------------------- %
%    [q]=evalc('dbstack'); q=[10 q(1:end-1)];
%     i=find(q==10); i=i(max(1,2+n):end);
%   % [q(i+1); q(i+2); q(i+3);; q(i+4); q(i+5)]'
%     j=find(q(i+1)=='>');
%     if ~isempty(j), S(j(1)).curr=1; end
% -------------------------------------------------------------------- %

  if vflag~=2, dispstack(S); else printf('\n'); end
  gStack=S; iStack=1;
end

  if ~isempty(df)
     N=numel(gStack); i=iStack+df;
     if i<1
        fprintf(1,'\n  already at lowest frame\n\n');
        return
     elseif i>N
        fprintf(1,'\n  already at highest frame\n\n');
        return
     end
   % WRN! not allowed by matlab:
   % ---> ERR Debug commands only allowed when stopped in debug mode.
   % if     df>0, evalin('caller',repmat('dbup; ',  1, df));
   % elseif df<0, evalin('caller',repmat('dbdown; ',1,-df));
   % end
     iStack=i; iS=i;
  else
     N=numel(S);
     if isempty(iS) && nargs
        if vflag>2 || N<=3
          iS=1:N;
        else 
          iS=1:3;
        end
     end
     if isempty(iS), return; end
  end

  if ~iscell(iS), iS={iS}; end
  for i=1:numel(iS), k=iS{i};
     if ischar(iS{i}), iS{i}=str2num(iS{i}); end
  end
  iS=flip(unique([iS{:}])); nS=numel(iS);
  if nS>1, dl=1; else dl=3; end

  for k=iS
     if ischar(k), k=str2num(k); end
     if i<1 || i>N
        wblog('ERR','index out of bounds (%d/%d)',i,N); 
        continue
     end

     e0=[char(27) '[38;5;8m' ];
     e2=[char(27) '[38;5;12m'];
     em=[char(27) '[0m'];

     s=sprintf('%2d: %s ',k, regexprep(S(k).file,'.*\/',''));
     if isempty(regexp(S(k).file,[S(k).name '\.m']))
        s=[s ':: ' S(k).name ' '];
     end

     if nS==1, s=[s ' ' char(10)];
     else s(end+1:80)='-'; end

     fprintf(1,[e0 '%s' em '\n'],s);

     l=S(k).line; l=sprintf('%d:%d',max(1,l-dl),l+dl);
     cmd=['dbtype(''' S(k).file ''',''' l ''')'];

     if ~wblog('--hl-check'), eval(cmd);
     else
        s=textscan(evalc(cmd),'%s','whitespace','\n'); s=s{1};
        s=regexprep(s,'^(\d)  ','  $1');
        s=regexprep(s,'^(\d\d) ',' $1');

        lpat=['^\s*' num2str(S(k).line) '\s'];
        fmt=['     %s\n'];
        for j=1:numel(s)
           if ~isempty(regexp(s{j},lpat))
              s{j}=[ e2 s{j} em];
           end
           fprintf(1,fmt,s{j});
        end
     end
  end
  fprintf(1,'\n');

end

