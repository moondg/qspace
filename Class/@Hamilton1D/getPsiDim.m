function [dd,dstr]=getPsiDim(HAM,varargin)
% function [dd,dstr]=getPsiDim(HAM [,opts])
%
%    get LR dimensions throughout entire MPS crrently linked to HAM.
%
% Options
%
%    -x    extended dimension D (i.e., states rather than multiplets)
%    -a    all, i.e., both [D*, D] (default: just D*)
%
% Wb,Apr11,14

  getopt('init',varargin);
     xflag=getopt('-x'); if ~xflag
     aflag=getopt('-a'); end
     qflag=getopt('-q');
  getopt('check_error');

  L=numel(HAM.mpo);
  dd=repmat(-1,L+1,2);

  for k=1:L
     q=load_dmrg_data(HAM,k,'HK','info');
     if ~isempty(q.HK.data)
        d=getDimQS(q.HK); o=itags2odir(QSpace(q.info.itags));
        if o==1, dd(k,:)=d(:,1);
        elseif o==2, dd(k+1,:)=d(:,1);
        end
     end
  end

  for k=find(dd(:,1)<0)'
     if dd(k)>0, continue; end

     q=load_dmrg_data(HAM,min(k,L),'AK');
     if isempty(q), continue; end

     d=getDimQS(q); if size(d,1)==1, d=[d;d]; end

     if k<=L
        dd(k,:)=d(:,1);
        if dd(k+1)>=0
           if ~isequal(dd(k+1,:),d(:,2)') && ~qflag, wblog('WRN',...
              'dimension mismatch of H and A-tensors (k=%g: %g/%g)',...
              k,dd(k+1),d(1)); 
           end
        else
           dd(k+1,:)=d(:,2);
        end
        if k==L, break; end
     else
        if any(dd(k-1,:)>d(:,1)')
           wbdie('dimension mismatch of A-tensors !?'); 
        end
        dd(k,:)=d(:,2);
     end
  end

  if nargout>1, dstr=get_dim_str(max(dd,[],1)'); end

  if xflag, dd=dd(:,2);
  elseif ~aflag, dd=dd(:,1);
  end

end

