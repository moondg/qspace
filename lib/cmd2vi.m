function cmd2vi(s)
% function cmd2vi(s)
% see also help2vi
% Wb,Feb03,10

  if nargin~=1, error('Wb:ERR','invalid usage'); end

  ftmp='/tmp/cmd2vi.tmp';
  fid=fopen(ftmp,'w');

  for i=1:nargin
     if i>1, fprintf(fid,'\n%s\n',repmat('=',1,72)); end
     cmd=[ 'evalin(''caller'',''' strrep(s,'''',''''''), ''')'];
     fprintf(fid,'%s',evalc(cmd));
  end

  fclose(fid);

  eval(sprintf(...
   '! sed ''s/<a[^<>]*>//g'' %s | sed ''s/<\\/a>//g'' > %s2 ; /bin/mv %s2 %s',...
    ftmp,ftmp,ftmp,ftmp...
  ));

  system(['unset LD_LIBRARY_PATH ; ' 10 ...
     sprintf('vim -R -c ''set ignorecase'' %s',ftmp)]);
  system(['/bin/rm ' ftmp]);

end

