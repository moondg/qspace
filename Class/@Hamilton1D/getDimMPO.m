function [dd,dstr]=getDimMPO(HAM,varargin)
% function [dd,dstr]=getDimMPO(HAM)
%
%    get largest LRs dimension for (pseudo) MPO HAM.mpo.
%
% Options
%
%    -x  return detailed dimensions for every tensor in the MPO.
%
% Wb,Aug28,15

  getopt('init',varargin);
     xflag=getopt('-x');
  getopt('check_error');

  L=numel(HAM.mpo);

  if using_full_MPO(HAM)
     dd=cell(1,L); for k=1:L, dd{k}=getDimQS(HAM.mpo(k)); end
     for k=[1 L]
        if size(dd{k},2)==3
           dd{k}=dd{k}(:,[1 1 2:end]); if k==1, i=1; else i=2; end
           dd{1}(:,i)=1;
        end
     end
     dd=cat(3,dd{:});
        if xflag, return; end
        if isempty(dd), dd=0; dstr={'undef','()'}; return; end

     dd={ max(dd(:,1:2,:),[],3), max(dd(:,3:end,:),[],3) };
     dd=[ dd{1}, max(dd{2},[],2) ];

  else
     dd=zeros(L,3);
     for k=1:L, dd(k,:)=size(HAM.mpo(k).M); end
     if xflag, return; end

     dd=max(dd,[],1);
  end

  if nargout>1
     dstr=get_dim_str([ max(dd(:,1:2),[],2), dd(:,3) ]);
  end

end

