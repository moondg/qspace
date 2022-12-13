function q=query_yesno(istr,varargin)
% function q=query_yesno(istr [,opts])
% Wb,Sep17,20

  q=0;

  getopt('init',varargin);
     defYes=getopt('--yes'); if defYes, q=1; end
     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end
  getopt('check_error');

  if isbatch, fprintf(1,...
    '\n  WRN ignoring query due to batch mode\n  WRN %s\n\n',istr);
     return
  end

  fs='';
  if vflag, S=dbstack(1);
     if ~isempty(S)
        fs=sprintf('%s:%d',regexprep(S(1).file,'.*\/',''),S(1).line);
        if vflag>1
             fs=['\n   ' fs '\n   '];
        else fs=[fs ' :: ']; end
     elseif vflag>1, fs='(workspace) ';
     end
  end

  if defYes
       s=['\n   ' fs istr '  { [1y] | 0n }  '];
  else s=['\n   ' fs istr '  { 1y | [0n] }  '];
  end

  i=input(s,'s');
  if ~isempty(i)
     if regexpi(i,'^([0n]|no)$$'), q=0;
     elseif regexpi(i,'^([1y]|yes)$'), q=1;
     elseif regexpi(i,'^\d'), q=str2num(i);
     else q=i; end
  end

end

