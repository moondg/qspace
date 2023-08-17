function rval=wbtic(task)
% function t=wbtic(task)
%
%    wbtic --start:#   start clock # = 1..9
%    wbtic --stop:#    stop clock #
%    wbtic --resume:#  resume clock #
%    wbtic --time:#    show current time 
%
%    wbtic --start:###  permits to start group of clocks
%    where only one clock can run at any given time;
%    (that is, a toc only is takend relative to the last
%    toc of any other clock in the group).
%
% Wb,Aug15,23

% see also wbtic.cc

  persistent IT

  q=(nargin && isequal(task,'--clear-all'));
  if q && (nargout || nargin>1), wbdie('invalid usage (%s)',task); end
  if isempty(IT) || q
     IT=struct('id',[],'fline','','group',[],'started',0,'tic',[],'toc',0,'tt',0);
     %   id       index, while redundant here, useful when returning clock
     %            via GET as this then also shows the id then.
     %   fline    file:line of last call to (re)start or toc
     %   group    if clocks are started simultaneously they act as group
     %            with effectively only one clock running at any time
     %   started  set to 0 if not started, or now() from the time started
     %   tic      reference used with toc (set when the clock is started)
     %   tt       always stores time differences
     IT=repmat(IT,1,9);
     for k=1:numel(IT), IT(k).id=k; end
     if q, return; end
  end

  nargs=nargin;
  if nargout && ~nargs, task='--get-all'; nargs=1; end

  if ~nargs && ~nargout
     kk=find([IT.started]);
     if isempty(kk)
        fprintf(1,'\n   (no clocks stared yet)\n\n');
     else
        for k=kk, q=IT(k);
           n={ numel(q.tt)-1, 's' }; if n{1}==1, n{2}=''; end
           fprintf(1,'   clock %d: %s (%d click%s)\n',k, ...
           datestr(q.started),n{:});
        end
     end
     return
  end

  if nargs~=1 || ~ischar(task) && ~isnumeric(task)
     if nargs, o={task}; else o={}; end
     if ~helpthis(nargout,o{:}), wbdie('invalid usage'); end
     return
  end

  if isnumeric(task) q={'get',task};
  else
     q=regexp(task,'^--(\w[\w-]+):(\d+)$','tokens');
     if ~isempty(q)
        q=q{1}; q{2}=double(q{2}-'0');
        if ~isempty(q{2})
           if any(q{2}<1 | q{2}>9)
              wbdie('invalid clock ids (%s)',task);
           end
        end
     else
        q=regexp(task,'^--(\w[\w-]+)$','tokens');
        if ~isempty(q), q=q{1}; q{2}=[];
        else wbdie('invalid task %s',task); end
     end
  end

  task_=task; task=q{1}; kk=q{2};

  if ~isempty(kk)
     switch task
        case {'start','restart'}
           k=unique([IT(kk).group]);
           if ~isempty(k), for k=k, IT(k).group=[]; end; end

           S=dbstack; if numel(S)>1, S=S(2);
                fline=sprintf('%s:%d',S.file,S.line);
           else fline=''; end

           for k=kk
              if ~IT(k).started || task(1)=='r'
                   IT(k).group=kk; IT(k).fline=fline;
                   IT(k).started=now; IT(k).tic=tic; IT(k).tt=0;
              else wblog('WRN','clock %d already started',k); 
              end
           end

        case {'toc'}, k=kk;
           if numel(k)~=1, wblog('ERR',...
              'one clock expected with toc (got%s)',sprintf(' %d',k)); end
           if isempty(IT(k).tic)
              wblog('ERR','clock %d not yet started',k);
           else
              S=dbstack; if numel(S)>1, S=S(2);
                   fline=sprintf('%s:%d',S.file,S.line);
              else fline=''; end

              kk=IT(k).group;
              tlast=max([IT(kk).toc]);
              IT(k).toc=toc(IT(k).tic);
              IT(k).tt(end+1)=IT(k).toc-tlast;
              IT(k).fline=fline;
           end

        case 'get'

           kk=unique(IT(kk).group);
           rval={IT(kk).tt}; if numel(rval)==1, rval=rval{1}; end

        case 'GET'

           kk=unique(IT(kk).group);
           rval=IT(kk);

        otherwise wbdie('invalid task %s',task_);
     end
  else
     switch task
        case 'get-all', k=find([IT.started]);
           rval={IT(k).tt};

        case 'GET-ALL', k=find([IT.started]); rval=IT(k);

        otherwise wbdie('invalid task %s',task_);
     end
  end

end

