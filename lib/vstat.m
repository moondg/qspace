function sout = vstat(v, varargin)
% function sout = vstat(v [,opts])
% 
%     statistics on data vector
% 
% Wb,Oct12,02 ; Wb,Aug18,18

  n=inputname(1); if ~isempty(n), n=[n ' = ']; else n=''; end
  s=size(v);
  v=reshape(v,[],1);

  av=mean(v); nrm=max(1,mean(abs(v)));
  dv=std(v);

  if dv==0
     fprintf(1,'\n   %srepmat(%g%s)\n\n',n,v(1),size_str(s));
     return
  elseif dv/nrm<1E-12
     fprintf(1,'\n   %srepmat(%g%s) @ %.3g\n\n',n,v(1),size_str(s),dv);
     return
  end

  q=std(diff(v));
  if q/nrm<1E-8
     fprintf(1,'\n   %slinspace(%g,%g,%g) @ %g\n\n',n,...
        v(1),v(end),numel(v),diff(v(1:2)));
     return
  end

  if all(v>0)
     q=std(diff(log(v)));
     if q/nrm<1E-8
        fprintf(1,'\n   %slogspace(%g,%g,%g) = %g * %.3g * %g\n\n',...
          n,log10(v(1)),log10(v(end)),numel(v),v(1),v(2)/v(1),v(end));
        return
     end
  elseif all(v<0)
     q=std(diff(log(-v)));
     if q/nrm<1E-8
        fprintf(1,'\n   %s-logspace(%g,%g,%g) = %g ** %.3g ** %g\n\n',...
          n,log10(-v(1)),log10(-v(end)),numel(v),v(1),v(2)/v(1),v(end));
        return
     end
  end

  vr = [min(v), max(v)]; vr(3)=diff(vr);

  s=num2str2(av,dv);
  if regexp(s,'+/-'), s=''; else s=['= ' s]; end

  fprintf(1,'\n   avg   : %.4g +/- %.4g %s\n',av,dv,s);
  fprintf(1,'   range : [ %.4g, %.4g ] @ %.3g\n\n', vr);

end

% -------------------------------------------------------------------- %

function s=size_str(sz)

   s=['[' sprintf(' %g',sz) ' ]'];

end

% -------------------------------------------------------------------- %

