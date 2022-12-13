function help2vi(s)
% function help2vi(s)
% Wb,Jan18,10

  if nargin~=1, error('Wb:ERR','invalid usage'); end

  ftmp='/tmp/help2vi.tmp';
  fid=fopen(ftmp,'w');
     for i=1:nargin
        if i>1, fprintf(fid,'\n%s\n',repmat('=',1,72)); end
        fprintf(fid,'%s',evalc([ 'help ' s ]));
     end
  fclose(fid);

  eval(sprintf(...
   '! sed ''s/<a[^<>]*>//g'' %s | sed ''s/<\\/a>//g'' > %s2 ; /bin/mv %s2 %s',...
    ftmp,ftmp,ftmp,ftmp...
  ));

  system(sprintf('unset LD_LIBRARY_PATH ; vim -R -c ''set ignorecase'' %s',ftmp));
  system(['/bin/rm ' ftmp]);

end

