function [fdir,isym]=get_fdir(A)
% function [fdir,isym]=get_fdir(A)
% Wb,May31,25

  if numel(A)~=1
     wbdie('invalid usage (singe QSpace required)'); end

  fdir='';
  if isfield(A.info,'fdir'), fdir=A.info.fdir; end

  if ~ischar(fdir)
     r=numel(A.Q); nsym=numsym(A);
     if numel(fdir)<=r || ~isint(fdir) || any(~fdir), fdir
        wbdie('invalid fdir');
     end

     is=fdir(r+1:end); fdir=fdir(1:r);
     if any(is<1) || any(diff(is)<=0) || nargin>2 && any(is>nsym), fdir
        wbdie('invalid fdir'); end
     is=['@' sprintf('%d',is) ];

     if all(abs(fdir)<=1) % keep compact by default
        fdir(fdir>0)='+';
        fdir(fdir<0)='-'; fdir=[fdir is];
     else
        fdir=sprintf(',%d',fdir);
        fidr=['[' fdir(2:end) ']' is];
     end
  end

  if nargout>1, isym=[];
     i=find(fdir=='@'); n=numel(i);
     if   n==1, isym=double(fdir(i+1:end)-'0'); fdir=fdir(1:i-1);
     else n || ~isempty(fdir), disp(fdir); wbdie('invalid fdir');
     end
  end

end

