function check_finished(nfin,task,cmd)
% function check_finished(nfin,task,cmd)
% Wb,May06,25

  persistent task_ cmd_ ncall err

  if ~isequal(cmd,cmd_) || isempty(err), ncall=0; err=0; end; cmd_=cmd;
  ncall=ncall+1;

% wblog('TST','nfin=%d, ncall=%d, err=%d, task=''%s'', cmd=''%s''',...
% nfin,ncall,err,task,cmd); 

  if ~nfin, task_=task; err=err+1;
     if err==1 % safeguard on infinite recursive calls
        evalin('caller',['clear task; ' cmd]) % rerun with default
     else task
        f='tstEnv'; u=getuser(0,f);
        if isempty(u), wbdie('invalid default task');
        else
           u.nerr=u.nerr+1;
           setuser(0,f,u);
        end
     end

  elseif err
     f='tstEnv'; u=getuser(0,f);

     if ~isempty(u)
        if ncall<=2
             u.nwrn=u.nwrn+1;
        else u.nerr=u.nerr+1; end
        setuser(0,f,u);
     else
        s={'',''};
        if ischar(task_), s{1}=['=''' task_ '''']; else disp(task_); end
        if ischar(task ), s{2}=[' ''' task  '''']; end

        [~,q]=wblog('--hl-check','WRN');
        if ncall<=2
           s={['reset task to default' s{2}],
              ['Hint: use ''clear task'' when switching between ' ... 
               'different example scripts']};
           fprintf(1,[q{1} '\n  WRN %s' q{2} '\n  %s\n\n'],s{:});
        else
           fprintf(1,[q{1} '\n  WRN invalid task%s, ' ... 
              'switched to default%s, instead.\n\n' q{2}],s{:});
        end
     end
     task_=''; err=0;
   end

end

