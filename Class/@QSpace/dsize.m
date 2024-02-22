function [D,DD]=dsize(A,varargin)
% function [D,[DD|I]]=dsize(A [,opts])
% get maximum D for every dimension of given MPS state
%
%     where D is a vector of length rank(A)
%
% Options
%
%   '-d'  returns detailed block sizes of single input A in DD with
%         CG space considered if exists (overall size returned in D)
%
% Wb,Aug08,06

  getopt('init',varargin);
     dflag=getopt('-d');
  varargin=getopt('get_remaining'); narg=length(varargin);

  k=[];
  if narg==1 && isnumber(varargin{1}), k=varargin{1};
  elseif narg
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if nargin<1
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  if ~isa(A,'QSpace')
     if isempty(A), D=[]; DD=[]; return; end
     if ~isstruct(A) || ~isfield(A,'Q') || ~isfield(A,'data')
        eval(['help ' mfilename]); return
     end
  end

  if dflag
     if numel(A)~=1
        wbdie('invalid usage (scalar A expected with -d)'); end

     r=numel(A.Q); if isempty(k), k=1:r; end
     for i=1:length(k)
        [QQ{i},dd,dc]=getQDimQS(A,k(i)); DD{i}=[dd'; prod(dc,2)];
        D{i}=[ sum(DD{i}(:,1),1), sum(prod(DD{i},2),1) ];
     end
     D=cat(1,D{:});

     nd=numel(A.data);
     if ~isempty(A.info), nq=size(A.info.cgr,2); else nq=0; end
     sd=ones(nd,r); sc=ones(nd,r,nq);
     for i=1:nd
        s=size(A.data{i}); sd(i,1:length(s))=s; for j=1:nq
        s=cgr_size(A,i,j); sc(i,1:length(s),j)=s; end
     end

     if any(sd(:)==0) || any(sc(:)==0)
     wblog('WRN','zeros encountered in size data'); end

     DD=add2struct('-',DD,QQ,sd,sc); clear QQ

     s=whos('A');
     if s.bytes>2^30, sfac=2^30; t='G';
     elseif s.bytes>2^20, sfac=2^20; t='M';
     elseif s.bytes>2^10, sfac=2^10; t='k';
     end

     f=sprintf('size_%s',t); DD=setfield(DD,f,s.bytes/sfac);
     f=sprintf('mem_%s',t); DD=setfield(DD,f, ...
       ( sum(prod(sd,2)) + sum(sum(prod(sc,2))) ) * (8/sfac) );
     f=sprintf('exp_mem_%s',t); DD=setfield(DD,f, ...
       ( sum(prod(sd.*prod(sc,3),2)) ) * (8/sfac) );
     return
  end

  N=length(A); r=length(A(1).Q);
  DD=zeros(N,r);

  for i=1:length(A)
     di=getDimQS(A(i)); s=size(di);
     DD(i,1:s(2))=di(1,:);
     if s(1)==2, DD(i,1:s(2),2)=di(2,:); end
  end

  D=squeeze(max(DD));

  if nargout==0
     n=inputname(1); if ~isempty(n)
     fprintf(1,'\nSize distribution of QSpace set %s = \n\n',n);
     else inl 1; end

     sA=size(A); m=prod(sA);
     for i=1:m
        si=vec2str(ind2sub(sA,i),'sep',',','fmt','%2g');

        sd=permute(DD(i,:,:),[3 2 1]);
        sd=mat2str2(sd,'fmt','%3g','nofac','rowsep','; ');

        fprintf(1,'  (%s): %s\n',si,sd);
     end
     inl 1

     clear D DD
  end

end

