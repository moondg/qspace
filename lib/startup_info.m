function startup_info(varargin)
% function startup_info(Iml)
% Wb,Jul18,19

  vflag=1;
  getopt('init',varargin);
     if getopt('-q'), vflag=0;
     elseif getopt('-v'), vflag=vflag+1; end
  I=getopt('get_last',[]);

  if isempty(I), I=mlinfo; end

  s=sprintf('pid=%d (%x)',I.pid,I.pid);

  n4={ setNumThreads(), 
       str2num(getenv('OMP_NUM_THREADS'))
       str2num(getenv('MKL_NUM_THREADS'))
       str2num(getenv('QSP_NUM_THREADS'))
  };
  n4(cellfun(@isempty,n4))={-1}; n4=[n4{:}];

  if any(n4>0) || isdeployed()
     s={s, sprintf(', nthreads=%g',n4(1)), '',''};

     q=max(n4(2:3));
     if q<0, s{3}=' ''''';
     elseif q~=n4(1), s{3}=sprintf(' %g',q); end

     if n4(4)>0
        s{4}=sprintf(' x%g_qsp',n4(4));
     end
     if ~isempty([s{3:4}]), s{2}=[s{2} ':']; end
     s=[s{:}];
  end

  q=str2num(getenv('ML_USING_PARPOOL'));
  if ~isempty(q) && q>1
     s=[s sprintf(' x%g_parpool',q)];
  end

  q=I.cgs.verbose;
  if ~isempty(q) && any(q)
     if numel(q)>1, v=q(1);
       q=dec2bin(q(end)); m=4; l=length(q);
       q=[ repmat('0',1,mod(-l,m)), q]; l=mod(size(q,2)+1,2);
       q=reshape(q,m,[]); q(end+1,:)='_'; q(end,l:2:end)=' ';
       v=sprintf('\n   verbose = %g ''%s''',v,reshape(q(1:end-1),1,[]));
     else s=[s, sprintf(', verbose=%g',q)]; v=''; end
  else v=''; end

  if isdeployed(), i={};
     if I.status.ismcc,      i{end+1}='mcc';       end
     if I.status.isdeployed, i{end+1}='deployed';  end
     if I.status.isbatch,    i{end+1}='batch';     end
     if usejava('desktop'),  i{end+1}='desktop';   end
     if ~I.display,          i{end+1}='nodisplay'; end

     if ~isempty(i)
        q=sprintf(', %s',i{:});
        l=[length(s), length(q)]; l(3)=sum(l);
        if l(3)<80, s=[s q]; else s=[s 10 '>> ' q(3:end)]; end
     end

  else
     if ~I.display, s=[s ', nodisplay']; end
  end

  q=get(0,'DefaultFigureRenderer');
  if ~isempty(q), s=[s ', ' q]; end

  q=getenv('MYMATLAB');
  if ~isempty(q)
     fprintf(1,'\r>> matlab startup with base %s\n>> having %s%s\n',...
     repHome(q),s,v);
  else
     fprintf(1,'\n>> matlab startup in pwd %s\n>> having %s%s\n\n',...
     repHome(pwd),s,v);
  end

end

