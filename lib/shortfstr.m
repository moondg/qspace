function str = shortfstr(str, len)
% expects an input string of the sort fname>subroutine:linenumber

  if length(str)>len
     if len<6, str=str(1:len); return; end

     i3 = findstr(str,':'); if ~isempty(i3), i3=i3(end); else i3=length(str)+1; end
     s3 = str(i3:end);

     i2 = findstr(str,'>');
     if ~isempty(i2), i2=i2(end);
     else
        i2 = findstr(str,'_');
        if ~isempty(i2), i2=i2(end)+1; else i2=i3; end
     end
     s2 = str(i2:i3-1);

     s1=str(1:i2-1);

     n=len-length(s3);
     n1=min([length(s1), max([ceil(n/3),n-length(s2)]) ]); n2=n-n1;

     if length(s2)==0
          m=ceil(n1/3);
          if length(s1)>n1, s1=[s1(1:m) '..' s1(end-n1+3+m:end)]; end
     else if length(s1)>n1, s1=[s1(1:n1-1) '''']; end
     end; if length(s2)>n2, s2=[s2(1:n2-1) '''']; end

     str = [s1 s2 s3];
  end

return

