function varargout=SymStore(sym,task,varargin)
% SymStore
%
%   Matlab class to analyze the RCStore for a specfic symmetry
%
% Usage by example:
%
% # SymStore SU4 dim   # show all generated irreducible multiplets
%   q=SymStore('SU4','dim');   # same as previous but also returns data
%
% # analysis of outer multiplicity for standard rank-3 CGCs
%   q=SymStore('SU4','mp3');
%
% # check mp3 integrity (e.g. in terms of IDs etc)
%   q=SymStore('SU4','--check-mp3');
%
% # check XData integrity (e.g. in terms of IDs etc)
%   q=SymStore('SU4','--check-X');
%
% # load specific RCStore data
%   q=SymStore('SU3','-d','...')      % show dimension for every multiplet mentioned
%   q=SymStore('SU3','-C','-f','...') % accesses CStore
%   q=SymStore('SU3','-C','++-','..') % load all CData that match qstr
%   q=SymStore('SU3','-M','...')      % load map3 data (in CStore)
%   q=SymStore('SU3','-R','...')      % accesses RStore
%   q=SymStore('SU3','-X','...')      % accesses XStore
%   q=SymStore('SU3','XOM'        )   % print XStore OM statistics
%   q=SymStore('SU3','rmX',qstr,ID)   % remove XStore entries that match qstr+ID
%
% # Permutation group
%   SymStore('PERM','213');           % dimension of given Young tableau
%
% # Get Casimir for given representation // Wb,Oct05,15
%   q=SymStore('SU3','Casimir', '...') 
%
% # generate list to rebuild RCStore e.g. on Cluster // Wb,Aug18,16
%   SymStore('SU4','--rebuildR');
%   SymStore('SU4','--rebuildC','++--',qloc);
%
%   SymStore('--xtr',f);      % translates xstore file into referenced CStore file
%   SymStore('--prc',f1,...); % print qset
%
% Examples
%
%    SymStore('SU3','rmX','(56,65)','-t')         % '-t' to test (dry-run)
%    SymStore('SU3','rmX','++','#ffffa50a','-t'); % '-t' to test (dry-run)
%
% Wb,Jan30,15

  if ~nargin, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  varargout=cell(1,max(1,nargout));

  if ischar(sym) && ~isempty(regexp(sym,'^-'))
     if nargin>1
          varargin={sym,task,varargin{:}};
     else varargin={sym}; end

     if isequal(varargin{1},'--check-gCX')
        check_gCX(varargin{2:end});
        varargout={}; return
     elseif isequal(varargin{1},'--xtr')
        [varargout{:}]=xfile2qout(varargin{2:end});
        return
     elseif isequal(varargin{1},'--prc')
        for i=2:numel(varargin)
           n=inputname(i); if isempty(i), n=sprintf('#%g',i); end
           print_qset(n,varargin{i},'-v');
        end
        varargout={}; return
     else wbdie('invalid usage');
     end
  end

  if isstruct(sym) && ...
     isfield(sym,'type') && isfield(sym,'qset') && isfield(sym,'qdir')
     if nargin>1, varargin={task,varargin{:}}; end
     [varargout{:}]=LoadCData(sym,varargin{:});
     return
  elseif isequal(upper(sym),'PERM') || isequal(upper(sym),'SYM')
     if nargin<2, wbdie('invalid usage'), end
     [varargout{:}]=sdim(task,varargin{:});
     return
  end

  if ischar(sym) && length(regexp(sym,'[()\s@,;]'))>2

     if nargin>1, varargin={task, varargin{:}}; end

     s=sym; ext=''; sym='';
     i=max(find(s=='.')); if ~isempty(i), ext=s(i+1:end); s=s(1:i-1); end
     s_=s; xpat=0;

     s=regexprep(s,'[\w\d]+:[\w\d]+\s*',''); % skip 'd:c3' etc..

     if regexp(s,'[_@]'), xpat=1;
        s=regexprep(s,'\s*@\s*','_');
     elseif regexp(s,'\)\s*\d+\**\s*(\(|$)'), xpat=1;
        s=regexprep(s,'\)\s*(\d+\**)\s*\(',')_$1 (');
        s=regexprep(s,'\)\s*(\d+\**)\s*$',')_$1');
     end

     if regexp(s,'S[Up](?(\d+)?')
        s=regexprep(s,'(S[Up])\(?(\d+)\)?',' $1$2 ');
        q=strread(s,'%s');

        I=regexp(q,'S[Up](?(\d+)?');
        mark=ones(size(q));

        for i=1:numel(q), if ~isempty(I{i}), mark(i)=0;
           s=regexprep(q{i},'[()]','');
           if isempty(sym), sym=s;
           elseif ~isequal(sym,s)
             wbdie('inconsistent symmetry !?\n   ERR %s\n',s_);
           end
        end, end

        q=q(find(mark));

        if regexp(sym,'^Sp'), r=str2num(sym(3:end))/2;
        elseif regexp(sym,'^SU'), r=str2num(sym(3:end))-1;
        else wblog('ERR','invalid sym=%s !?',sym); end

        q=q'; q(2,1:end-1)={' '}; q=[q{:}];
     else
        q=regexprep(s,'\s+',' ');
     end

     s=regexprep(q,'_[\d\w]+\**\s*',' ');
     s=strread(regexprep(s,'[()\*,;_]',' '),'%s');
     rq=size(cat(1,s{:}),2);
     if isempty(sym), r=rq;
        sym=sprintf('SU%g',rq+1);
     elseif ~isequal(r,rq)
        wbdie('invalid input (%g/%g) !?\n   ERR %s',s_,rq,r);
     end

     if xpat
        q=regexprep(q,'\*\s+(','*(');
        [varargout{:}]=LoadXData(sym,q,varargin{:});
        return
     end

     if regexp(ext,'^m')
        [varargout{:}]=LoadMData(sym,q,varargin{:});
        return
     end

     if regexp(q,'\*')
        s=strread(regexprep(q,'[,()]',' '),'%s');
        I=regexp(s,'\*$');
        for i=1:numel(I), if isempty(I{i}), I{i}=0; end, end
        I=[I{:}]; i=find(diff(I));
        if numel(i)==1
           s=regexprep(s,'\*$','')';
           s(2,1:end-1)={','}; s(2,i)={';'};
           q_=q; q=[s{:}];
        elseif numel(i)>1, q, wbdie('got unsorted QSet !?'); 
        end
     end

     [varargout{:}]=LoadCData(sym,q,varargin{:});
     return
  end

  if nargin<2, task='dim'; 
     varargout=cell(1,nargout);
  end

  if nargin<1 || ~ischar(task)
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  if ~isempty(regexp(task,'\.mp3'))
     varargin=[{task}, varargin]; task='-M';
  elseif ~isempty(regexp(task,'\.cgd'))
     varargin=[{task}, varargin]; task='-C';
  end

  sym=regexprep(sym,'(S[pU])\(*(\d+)\)','$1$2');
  [D,s,err]=get_dir(sym);
  if err || isempty(D)
     r=get_rank(s); s='';
     if r<=3
        s=sprintf('\nhint: consider running ''setupRCStore %s''',sym);
     elseif r>5
        s=sprintf('\nwarning: asking for large rank-%d symmetry',r);
     end
     wbdie('-q',['symmetry %s not yet setup in RC_STORE' s],sym);
  end

  switch task
     case 'mp3', [varargout{:}]=getOMmp3 (sym,varargin{:});
     case 'XOM', [varargout{:}]=getSymXOM(sym,varargin{:});
     case {'-C','loadC'}, [varargout{:}]=LoadCData(sym,varargin{:});
     case {'-M','loadM'}, [varargout{:}]=LoadMData(sym,varargin{:});
     case {'-R','loadR'}, [varargout{:}]=LoadRData(sym,varargin{:});
     case {'-X','loadX'}, [varargout{:}]=LoadXData(sym,varargin{:});
     case 'Casimir', [varargout{:}]=CasimirC2(sym,varargin{:});

     case {'-d'}, [varargout{:}]=getSymDim(sym,varargin{:});
     case 'dim', [varargout{:}]=getSymDim(sym,varargin{:});
                 if ~nargout, varargout={}; end

     case '--plotR', plotRStat(sym,varargin{:}); varargout={};

     case 'rmX',        [varargout{:}]=remXData(sym,varargin{:});
     case 'findID',     [varargout{:}]=findID(sym,varargin{:});

     case '--rebuildR', rebuildRStore(sym,varargin{:}); varargout={};
     case '--rebuildC', rebuildCStore(sym,varargin{:}); varargout={};

     case '--check-mp3',[varargout{:}]=check_mp3(sym,varargin{:});
     case '--check-X',  [varargout{:}]=check_XStore(sym,varargin{:});

     case '--tensor-prod',  [varargout{:}]=tensorprod_QSet(sym,varargin{:});
     case '--contract',     [varargout{:}]=contractSym(sym,varargin{:});

     otherwise
     wbdie('invalid task ''%s''',task);
  end

end

