function s=sysctl(q)
% function s=sysctl(parameter)
% Wb,Dec05,24

  if ~nargin || ~ischar(q), wbdie('invalid usage'); end

  is_unix=isunix();
  if is_unix
       cmd=['getconf -a | egrep ' q];
  else cmd=['sysctl ' q];
  end

  [e,s]=system(cmd); if e, wbdie(['''' cmd ''' returned e=%d'],e); end

  ss=textscan(s,'%s','whitespace','\n'); ss=ss{1};
  n=numel(ss); ff=cell(2,n);
  for i=1:n
     s=ss{i}; if is_unix, s=regexprep(s,'\s{2,}',': '); end
     j=regexp(s,': ');
     if ~isempty(j)
        s={ s(1:j-1); s(j+2:end) }; ss(i,2:3)=s;
        if ~isempty(regexp(s{2},'^\d[\.\d]*$')), s{2}=str2num(s{2}); end
        ff(:,i)=s;
     else wblog('WRN','skipping unexpected entry %s',s); 
     end
  end

  if n>1
     if nargout, 
        for i=1:n, ff{1,i}=regexprep(ff{1,i},'\.','_'); end
        s=struct(ff{:});
     else
        fprintf(1,'\n');
        for i=1:n
             fprintf(1,'   %-32s %12s\n',ss{i,2:3});
        end; fprintf(1,'\n');
        clear s
     end
  else s=ff{end};
  end

end

