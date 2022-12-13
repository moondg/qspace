function [Il,RR,isw]=fix_Il_data(HAM,varargin)
% function [Il,RR,isw]=fix_Il_data(HAM,Il,RR)
%
%    This is quick fix for restarted sweeps with incomplete Il data.
%    => set empty Il data based on HAM(k).info
%
%    If Il is not specified or empty, all of Il will be reloaded.
%
% Wb,Jun16,17 : added RR data
% Wb,Sep08,15

   L=numel(HAM.mpo); Il=[]; RR=[]; k=[]; ISW=cell(L-1,1); e=0;

   for i=1:numel(varargin)
      if isfield(varargin{i},'e0')
         if isempty(Il), Il=varargin{i}; else e=bitxor(e,1); end
      elseif isfield(varargin{i},'data')
         if isempty(RR), RR=varargin{i}; rflag=1; else e=bitxor(e,2); end
      else k(end+1)=i; end
   end

   if e
      helpthis, if nargin || nargout
      wberr('invalid usage (%g)',e), end, return
   end

   varargin=varargin(k); aflag=0; is_=0;
   if ~isempty(k)
      getopt('init',varargin);
         aflag=getopt('-a');
      is_=getopt('get_last',is_);
   end

   rflag=bitxor(~isempty(RR),2*(nargout>1));
   if rflag==2, RR=QSpace(L-1,2);
   elseif rflag==1
      rflag=0;
   end

   aflag=bitxor(isempty(Il),2*aflag);
   if aflag
      if aflag<2, Il=repmat(struct(...
         'E0',[],'e0',[],'se',[],'rr',[],'r2',[],...
         'NK',[],'Dk',[],'Eg',[],'ex',[]),L-1,2);
      end
   elseif size(Il,1)+1<L, Il(L-1,2).e0=[]; end

   nfix=0;
   NPsi=load_dmrg_info(HAM,'NPsi');

   for k=1:L-1
      if aflag || numel([Il(k,:).e0])<2
         I=load_dmrg_data(HAM,k,'info');

         q=I.ek;
         if iscell(q) && ~isempty(q)
            q=q{end};
            if numel(q)>1, q=sum(q(:,end)); end
         end
         if isempty(q), q=nan; end

         R=I.Rho(end);
         if NPsi, R=R/NPsi; end

         I.E0=min(I.dsw(:,3));
         I.e0=q;
         I.se=SEntropy(R);
         I.Dk=getDimQS(R); I.Dk=I.Dk(:,2);
         I.NK=I.Dk(1);

         nfix=nfix+1;
         if nfix==1 && ~aflag
            wblog('NB!','\Nfixing missing Il data ...');
         end

         if isfield(I,'isw'), isw=I.isw; ISW{k}=isw; else isw=0; end

         fprintf(1,'\r  %4g Il(%g/%g, 2; %g) ... \r',nfix,k,L,isw);

         for j=1:2
            if aflag>1
               Il(k,j)=I;
            elseif isempty(Il(k,j).e0)
               Il(k,j)=setfields(Il(k,j),I,'-k');
            end
         end

         if k<=1, continue; end

         if isfield(I,'Eg')
            l=numel(I.Eg);
            if isw && fix((l+1)/2)~=isw, wblog('WRN',...
               'got isw inconsistency (k=%g): Eg @ %g/%g ',k,l,I.isw);
            end
            if l>1
               if is_, j=2*is_; else j=l; end
               Il(k,1).Eg=I.Eg(j-1);
               Il(k,2).Eg=I.Eg(j);
            end
         end

         if rflag && l>1
            l=numel(I.Rho);
            if isw && fix((l+1)/2)~=isw, wblog('WRN',...
               'got isw inconsistency (k=%g): Rho @ %g/%g ',k,l,I.isw);
            end
            if is_, j=2*is_; else j=l; end
            if mod(l,2), j=[j,j-1];
            else         j=[j-1,j]; end

            RR(k-1,1)=I.Rho(j(1));
            RR(k,  2)=I.Rho(j(2));
         end
      end
   end

   if nfix, fprintf(1,'\n'); end

   isw=unique([ISW{:}]);
   if numel(isw)>1
      wblog('WRN','got mixed sweep data (nsw=[%s])',vec2str(isw));
   end

   if rflag
      for k=1:L-1
         q=[ 0, isempty(RR(k,2)) ];
            if k>1 && isempty(RR(k-1,1)), q(1)=1; end
         if all(~q), continue; end

         nfix=nfix+1;
         fprintf(1,'   ..  RR(%2g,%g)\n',k,L);

         I=load_dmrg_data(HAM,k,'info');

         if q(1)
            l=length(I.Rho); if mod(l,2)==0, l=l-1; end
            if l>0, RR(k-1,1)=I.Rho(l); end
         end
         if q(2)
            l=length(I.Rho); l=l-mod(l,2);
            if l>0, RR(k,2)=I.Rho(l); end
         end
      end
   end

   if nfix, fprintf(1,'--> loaded/fixed %g iterations.\n\n',nfix); end

end

