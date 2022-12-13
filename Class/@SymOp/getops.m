function op=getops(G,varargin)
% function op=getops(G)
% Options
%
%   -f     ensure full matrices in returned object
%   -hc    return *.hc matrices instead of default *.op 
%   -cat3  concatenate matrices along 3rd dimension
%   -catV  concatenate vectorized matrices along 2nd dimension
%
% Wb,Nov17,11

  getopt('init',varargin);
     fflag=getopt('-f');
     hcflag=getopt('-hc');
     cflag=getopt('-cat3');   if cflag, fflag=1; else
     cflag=getopt('-catV')*2; if cflag, fflag=1; end, end
  getopt('check_error');

  op=cell(size(G)); nG=numel(G);

  if hcflag
       for i=1:nG, op{i}=G(i).hc; end
  else for i=1:nG, op{i}=G(i).op; end, end

  if fflag
     for i=1:nG, op{i}=full(op{i}); end
  end

  if cflag
     op=cat(3,op{:}); if cflag>1, s=size(op);
     op=reshape(op,s(1)*s(2),[]); end
  end

end

