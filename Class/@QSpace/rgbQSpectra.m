function [cc,Iout]=rgbQSpectra(HH,varargin)
% function [cc,Iout]=rgbQSpectra(H,q3 [,opts])
%
%   Encode given energy spectrum to rgb specturm
%   by using thermal weights in specified (or otherwise dominant)
%   symmetry sectors
%
% Options
%
%   'beta',.. effective inverse `temperature' to use to compute
%            `thermal' weight in symmetry sectors (2)
%   '-R'      assume input is density matrix which already gives
%             weights (with beta ignored here)
%
% Wb,Sep09,21

% NB! cc may be rotated into any other color basis in caller
% e.g. by cc -> cc * c3 with c3 some 3x3 other color specs.

   vflag=1;
   getopt('init',varargin);
      Rflag=getopt('-R');     if Rflag, beta=[]; else
      beta =getopt('beta',2); end

      if     getopt('-q'), vflag=0;
      elseif getopt('-v'), vflag=2; end

      q3=getopt('q3',[]);

   if isempty(q3)
        q3=getopt('get_last',[]);
   else getopt('check_error');
   end

   nH=numel(HH); w=0; e=0; mQ=ones(1,nH);
   for k=1:nH, Hk=HH(k); if isempty(Hk), mQ(k)=0; continue; end
      if isdiag(Hk,'-d')~=2, if vflag
         fprintf(1,'\r %4g/%g %s calling eig() ... \r',k,nH,mfilename); end
         [~,q]=eigQS(Hk); HH(k)=QSpace(q.EK); w=w+1;
      end
      Q=Hk.Q;
      if numel(Q)~=2 || ~isequal(Q{:}), e=e+1; end
   end
   if w && vflag, fprintf(1,'\r%60s\r',''); end

   if e
      if nH>1, s=sprintf(' (%g/%g)',e,nH); else s=''; end
      wbdie(['expecting block-diagonal scalar input QSpace operators' s]);
   elseif w && vflag
      s={'diagonalized',' input QSpace operator'};
      if nH>1, s{1}=[s{1} sprintf(' %g/%g',w,nH)]; s{2}(end+1)='s'; end
      wblog('WRN',[s{:}]);
   end

   if isempty(q3), n=4;
      Q={}; for i=find(mQ,n,'last'), Q{end+1}=HH(i).Q; end
      Q=[Q{:}]; Q=uniquerows(cat(1,Q{:}));
      [q,i]=sort(sum(Q.^2,2));
      q3=Q(i([2 1 3]),:);
   end

   cc=zeros(nH,3);

   for k=1:nH, Hk=HH(k); if isempty(Hk), continue; end
      [qk,~,dc]=getQDimQS(Hk,1); dc=prod(dc,2);
      if any(dc>1)
         [i,j]=matchIndex(Hk.Q{1},qk,'-s');
         dc=dc(j);
      end

      q=Hk.data;
      if Rflag
         for i=1:numel(q), q{i}=dc(i)*sum(q{i}); end
      else
         qmin=min([q{:}]);
         for i=1:numel(q), q{i}=dc(i)*sum(exp(-beta*(q{i}-qmin))); end
      end

      w=[q{:}]'; q=sum(w);
      if Rflag
         if abs(q-1)>1E-9, w=w/q; wblog('WRN',...
           'got non-normalized density matrix (1%+.3g)',q-1); 
         end
      else w=w/q;
      end

      qk=Hk.Q{1}; 
      if isempty(q3)
         if size(qk,1)<3, wbdie(...
           'fewer than 3 symmetry sectors encountered (need to specify q3)');
         end
         [w,is]=sort(w,'desc');
         q3=qk(is(1:3),:);
         cc(k,:)=w(1:3);
      else
         [i,j]=matchIndex(q3,qk);
         cc(k,i)=w(j);
      end
   end

   if nargout>1
      Iout=add2struct('-',q3,Rflag,beta);
   end
end

% -------------------------------------------------------------------- %

