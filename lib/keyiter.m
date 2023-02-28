function k=keyiter(k,varargin)
% function k=keyiter(k[,kmax])
% get next k (e.g. useful for interactive plot session)
%
%     '-c'                 clear preexisting keys from getchar()
%
%     up or right arrow  : k++
%     left or down arrow : k--
%     <Home>             : kmax
%     <End>              : 1
%     number <cr>        : k=number
%     q or Q             : k=-1 (quit)
%
%     otherwise          : k++
%
% Further usages:
%
%   if next call in loop to keyiter should proceed automatically:
%   if (...), keyiter('next'); continue;  end
%
% Example usage:
%
%   while 1, ..., k=keyiter(k); if k<1 || k>kmax, break; end
%   k=1:15; while keyiter(k), ...; end  % k is altered here in caller!
%
% Wb,Jan09,09

   persistent kold kmin kmax krange kflag fflag; k__=kold;

   if nargin==1 && ischar(k)
      if ~nargout && isequal(k,'next'), kflag=+1; clear k
      elseif isequal(k,'--status')
         k=add2struct('-',krange,kmin,kmax,kold,kflag,fflag);
      end
      return
   end

   if isempty(kmin), kmin=1; fflag=0; end

   if nargin && isnumeric(k) && numel(k)>1
      if nargin>1 && isnumeric(varargin{1}), wberr('invalid usage'); end

      krange=k; kmin=1; kmax=numel(k); kold=1; k=krange(kold);
         vn=inputname(1);
         if ~isempty(vn), assignin('caller',vn,k); end
      return
   end

   if ~isempty(krange)
      if krange(kold)==k, k=kold;
      else wberr('invalid usage'); end
   end

   if nargin==1 && ischar(k) && nargout==0
      varargin={k}; k=kmin; fprintf(1,...
      '\n   clear keyiter: press <ctrl-c> ...\n\n');
   else
      if isempty(kold), fprintf(1,...
      '\n   keyiter: press keys to navigate (see help) ...\n\n'); end
      kold=k;

      if nargin>1 && isnumeric(varargin{1})
      kmax=varargin{1}; varargin(1)=[]; end
   end

   drawnow

   u=getuser(groot,'NO_PAUSE');
   if ~isempty(u) && u, k=k+1; return; end

   tic; dt=0;
   if ~isempty(varargin)
      if ~isequal(varargin{1},'-c')
      wberr('invalid usage'); end
      dt=1E-2;
   end

   if ~isempty(kflag), c=kflag; kflag=[];
   elseif fflag, c='flush';
   else
      try
         while 1
            c=getkey;
            t=toc; if t>=dt, break; else tic; end
         end
      catch l
         if dt, return, end
         dispstack(l); rethrow(l);
      end
   end

   if isequal(c,'Left') || isequal(c,'Down') || isequal(c,127), k=k-1;
   elseif isequal(c,'Home'), k=kmin;
   elseif isequal(c,'End'), k=kmax;
   elseif isequal(c,'Esc'), k=0;
   elseif numel(c)==1 && (c==double('a') || c==double('f')), fflag=1;
   elseif numel(c)==1 && (c==double('q') || c==double('Q')), k=0;
   elseif numel(c)==1 && (c==double('p') || c==double('t'))
      if ~isempty(k__), k=k__; end
   elseif numel(c)==1 && c>=double('0') && c<=double('9')
      cc=c;
      while 1, c=getkey;
         if c==13, break; end
         if numel(c)==1 && c>=double('0') && c<=double('9')
            cc=[cc,c];
         else
            wblog('ERR','keyiter - invalid number (%g)',c(1));
            return
         end
      end
      k=str2num(char(cc));

   else k=k+1; end

   if ~isempty(kmax) && k>kmax, k=0; end
   if ~isempty(kmin) && k<kmin, k=0; end

   if ~isempty(krange) && k>=kmin && k<=kmax
   kold=k; k=krange(k); end

   if isempty(kmax)
      wblog('WRN','got k=%d empty/initialized kmax',k);
   end

      vn=inputname(1);
      if ~isempty(vn), assignin('caller',vn,k); end

 % while 1
 %    c=getkey;
 %    if isequal(c,'Left') || isequal(c,'Down'), k=k-1; break;
 %    elseif isequal(c,'Right') || isequal(c,'Up'), k=k+1; break;
 %    else wblog('ERR','invalid key (use cursors)');
 %    end
 % end

end

