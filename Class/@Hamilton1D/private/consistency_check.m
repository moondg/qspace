function [e,s]=consistency_check(HAM,varargin)
% function [e,s]=consistency_check(HAM [,opts])
%
%    this includes minor formatting corrections,
%    if output argument is specified
%
% Wb,Apr09,14

% outsourced from Hamilton1D
% NB! many consistency checks are already performed in setup_mpo.m
% Wb,Aug26,15

  getopt('init',varargin);
     tflag=getopt('-t');
     qflag=getopt('-q'); if qflag, tflag=tflag+1; end
  getopt('check_error');

  if ~isvector(HAM.mpo)
     wbdie('inconsistent H (mpo not a vector)'); end
  if nargout>1, qflag=1; end

  L=numel(HAM.mpo); nops=size(HAM.oez,1)+size(HAM.ops,1);

  if isstruct(HAM.mpo)
     for k=1:L
       s=size(HAM.mpo(k).M); s(end+1:3)=1;
       if s(3)~=nops, wbdie(['inconsistent H ' ... 
          '(invalid mpo-dimension [%g: dim3=%g/%g])'],k,s(3),nops);
       end
       if k>1, s1=size(HAM.mpo(k-1).M); else s1=s([2 1]); end
       if k<L, s2=size(HAM.mpo(k+1).M); else s2=s([2 1]); end

       if k>1
          if s(3)~=d || length(s)>3, wbdie( ...
            'invalid mpo (got changing d_op=%g/%g @ k=%g)',s(3),d,k);
          end
       else
          d=s(3); if d<2
            wbdie('invalid mpo (d_op>1 required; got %g)',d); 
          end
       end

       if s1(2)~=s(1) || s(2)~=s2(1), wbdie(['inconsistent H ' ... 
          '(invalid mpo-dimension [%g: (..x%g) * (%gx%g) * (%gx..)])'],...
          k,s1(2),s(1),s(2),s2(1));
       end
     end
  end

  ie=[isempty(HAM.store); isempty(HAM.mat)]; e=0; s='';
  if all(ie)
     ie=[ ie
        isempty(HAM.mpo)
        isempty(HAM.oez)
        isempty(HAM.ops)
     ]; if all(ie), e=-1; return, end
  end

  if ~diff(ie(1:2))
     wbdie('invalid usage (either HAM.store or HAM.mat must be set)');
  elseif ~ie(1)
     if ~isvarname(HAM.store)
       wbdie('got invalid name for global storage ''%s''',HAM.store); 
     end
  else
     i=find(HAM.mat=='/'); s='';
     if isempty(i), e=1; s=sprintf(...
        'invalid mat storage specification ''%s'' (at least one / expected)',...
        HAM.mat);
     elseif ~isdir(HAM.mat(1:i(end))), e=2; s=sprintf(...
        'directory for mat storage does not yet exist (%s [%s])', ...
         HAM.mat(1:i(end)), HAM.mat(i(end)+1:end)); 
     end

     if ~isempty(s) && ~qflag
        if ~tflag, wbdie('%s',s); end
     end
  end

end

% -------------------------------------------------------------------- %

