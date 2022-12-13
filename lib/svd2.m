function varargout=svd2(A,varargin)
% function svd2(A,varargin)
%
%    svd() got a problem with matrix of dimensions 130x65
%    e.g. a=randn(130,65); a=a/norm(a(:));
%    [u,s,v]=svd(a); e=norm(u*s*v'-a,'fro');
%    => e as large as 1E5 for matlab/2013b!
%
%    matlab/2013a does not have this problem
%    also all later versions (<= matlab/2015b)
%    have a similar problem! [see ML/Bugs/svd_bug.m and svd_bug.pdf].
%
% Wb,Feb16,16
% -------------------------------------------------------------------- %

% -------------------------------------------------------------------- %
% function svd2(A,varargin) // OLD VERSION
%
%    For ill-conditioned matrices, matlab/2013a can issue
%    error "SVD did not converge" @_*(&_!#@
%    Apparently, a problem that was not there with Matlab 2012!
%
% => see Archive/svd2_160216.m
%
% Wb,Mar27,15
% -------------------------------------------------------------------- %

  if nargout<=1, 
       varargout=cell(1,1);
  else varargout=cell(1,max(3,nargout));
  end

  e=0;

  try
     [varargout{:}]=svd(A,varargin{:});

     if nargout<=1
          ns=norm(varargout{1});
     else ns=norm(diag(varargout{2}));
     end

     nA=norm(A,'fro');
     e=abs(1-ns/nA); if e<1E-10, return; end

  catch l
     s=strrep(l.message,char(10),', ');
     wblog('ERR',['\N''%s''\nsaving data to file, ' ...
       'trying transpose ...\n%s/tmp-1.mat'],s,pwd);
     dispstack(l); whos A, save ./tmp-1.mat
  end

  if e
     s=size(A); s=sprintf(['SVD failed having %gx%g (@ e=%g !?) - ' ...
      'retrying with transpose ...'],s,e);
     warning('Wb:WRN','WRN\n\n> %s',s);
  else
     nA=norm(A,'fro');
  end

  try
     [varargout{:}]=svd(A',varargin{:});
  catch l
     wblog('ERR','saving data to file\n%s/tmp-2.mat',pwd);
     dispstack(l); save ./tmp-2.mat
     rethrow(l);
  end

  if nargout<=1
     ns=norm(varargout{1});
     e=abs(1-ns/nA);
  else
     varargout={varargout{3}, varargout{2}', varargout{1}};

     A_=varargout{1}*varargout{2}*varargout{3}';
     e=norm(A_-A,'fro')/nA;
  end

  if e>1E-8, wberr('SVD also failed for transpose (e=%g !?)\n',e);
  else fprintf(1,'> SVD now ok (@ e=%.3g)\n\n',e); end

end

