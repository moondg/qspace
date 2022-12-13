function [HH,mm]=convert_HS(HS,mmax)
% function [HH,m]=convert_HS(HS [,mmax])
%
%    Convert Hamiltonian encoding by structure HS
%    to older numeric format HH with 5 columns that only can
%    deal with up to 2-site terms in the Hamiltonian.
%
%    m returns number of terms in individual Hamiltonian term.
%    Show warnings for m>mmax (default: mmax=2).
%
% See also convert_HH.m
% Wb,Nov25,20

% outsourced from @Hamilton1D/plot.m // Wb,Mar08,22

  if ~isstruct(HS) || ~isfield(HS,'iop'), wbdie('invalid usage'); end

  nH=numel(HS); mm=zeros(1,nH);
  for i=1:nH, q=HS(i); m=size(q.iop,1); mm(i)=m;
     if     m==1, HS(i).HHrec=[q.iop, q.iop, q.hcpl];
     elseif m==2, HS(i).HHrec=[reshape(q.iop',1,[]), q.hcpl];
     end
  end
  HH=cat(1,HS.HHrec);

  if nargin<2, mmax=2; end
  if any(mm>mmax), [m,I,D]=uniquerows(mm');
     for i=find(m>mmax)'
        if D(i)>1, wblog('WRN',['skipping %d ' ... 
            '%d-site terms in HS(%d..%d) for HH'],D(i),m(i),I{i}([1,end]));
        else wblog('WRN','skipping %d-site term HS(%d) for HH',m(i),I{i});
        end
     end
  end
end

