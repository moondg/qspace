function [dfac,sw]=get_cgr_fac(A,i,sflag)
% function [dfac,sw]=get_cgr_fac(A,i [,opts])
%
%    CGData is no longer stored with the QSpaces themselves,
%    only their weights cgw within outer multiplicity space.
%    Therefore CGData is globally normalized to 1. up to
%    outer multiplicity.
%
% Options
%
%    '-s'  also return range as string info, or 
%    '-S'  same as -s, but return rats/wbrat expresssions
%          for info.cgr.cgw normalization
%
% Wb,Oct17,14

% adapted from fromer get_cgr_range.m

  if nargin<2 || numel(A)~=1 || numel(i)~=1
     wbdie('invalid usage');
  end

  if nargin<3, sflag=0;
  elseif isequal(sflag,'-s'), sflag=1;
  elseif isequal(sflag,'-S'), sflag=2;
  else wbdie('invalid usage'); end

  if ~gotCGS(A), q=1;
  else
     if ~isfield(A.info,'cgr') || ~isfield(A.info.cgr,'size')
        wbdie('invalid info.cgr data'); end

     cgr=A.info.cgr(i,:); q=ones(1,numel(cgr));
     for j=1:numel(cgr)
        w=cgr(j).cgw; if ~isa(w,'double'), w=mpfr2dec(w); end
        if isempty(w)
           if ~isempty(cgr(j).type) || ~isempty(cgr(j).qset) || ...
              ~isempty(cgr(j).qdir), wbdie(...
              'unexpected CGRef data (assuming CGR_ABELIAN)');
           end
           q(j)=1;
        elseif numel(w)==1, q(j)=w;
        elseif norm(w-w(1)*eye(size(w)))<1E-12, q(j)=w(1);
        elseif isvector(w), q(j)=norm(w);
        else
           s=size(w); if numel(s)>2 || diff(s)>0
              s=sprintf('x%d',s);
              wberr('unexpected cgw(%d,%d) data size %s',i,j,s); end

           w2=w'*w;
           e=norm(w2-w2(1)*eye(size(w2)));
           if e>1E-12, disp(w2), wblog('WRN',...
              'unexpected normalization cgw''*cgw(%d,%d) @ %.3g',i,j,e); 
           end
           q(j)=sqrt(w2(1));
        end
     end
  end

% dfac=1/sqrt(prod(q(1,:)));
% dfac is used by display.m => don't apply cgc factor on data sector
% as this may lead to confusion when actually looking up A.data{i};
% NB! also adapted NormCGC such that rank-2 CGC's are always stored
% as identity matrix (still), yet rank>2 CGC's have |cgc|^2 = 1
% (or delta_i,j in outer multiplicity i and j). % Wb,Oct18,14
  dfac=1;
  srd=char(8730); % = surd = char(hex2dec('221A')) // see also rat2.m

  if nargout>1 || nargin>=3 && sflag
     if norm(q-1)>1E-12, q=prod(q);
        if numel(q)==1
           if sflag>1 && q && abs(q)~=1
              [sw,~]=wbrat(q); r=0;
              sw=regexprep(sw{1},'sqrt\((.*)\)(?@r=1;)',[srd '$1']);
              if ~r, sw_=sw;
                 sw=[ srd num2str(q*q)]; if q<0, sw=['-' sw]; end
                 if length(sw)>4, sw=sw_; end
              end
           else sw=sprintf('%.5g',q);
           end
        elseif norm(diff(q))>1E-12
             sw=sprintf(' * %.4g',q); sw=sw(4:end);
        else sw=sprintf('%.4g (x%d)',q(1),numel(q));
        end
     else sw='';
     end
     if nargin<2, q=sw; end
  end

end

