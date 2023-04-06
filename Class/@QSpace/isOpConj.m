function q=isOpConj(A,varargin)
% function i=isOpConj(A [,output_type])
%
%    Check whether input A is compatible with irop structure
%    in which case the direction of the irop index is returned
%    (with optional 2nd argument, a string is returned, instead).
%
% Wb,Jun06,19

  getopt('init',varargin);
     qflag=getopt('-q');
  w=getopt('get_last',[]);

  r=numel(A.Q); q=0;

  if r>=2 && r<=3, d=getqdir(A);
     if d(1)==d(2), q=0;
        if ~qflag
           t=getitags(A); s=sprintf(';%s',t{:});
           wblog('WRN','unexpected operator (%s)',t(2:end));
        end
     elseif r==3, q=d(3); end
  end

  if ~isempty(w)
     switch w
        case '-dag', if q<=0, q=''; else q='^dagger'; end
        case '-s', if q<0, q='*'; else q=''; end
        otherwise, w, wblog('WRN','invalid switch');
     end
  end

end

