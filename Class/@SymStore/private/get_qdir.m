function qdir=get_qdir(qstr);
% function qdir=get_qdir(qstr);
% Wb,Aug27,16

  qdir=regexprep(qstr,'-','');
  qdir=regexprep(qdir,'([\w\d]+)\*,*','-');
  qdir=regexprep(qdir,'([\w\d]+),*','+');

  i=find(qdir==';'); n=numel(i);

  if n==1
     if ~isempty(find(qdir=='-'))
        wbdie('unexpected QSet string\n    ERR %s',qstr);
     end
     qdir=[ qdir(1:i-1), regexprep(qdir(i+1:end),'\+','-') ];
  elseif n>1
     qstr, wbdie('invalid input string'); 
  end

end

