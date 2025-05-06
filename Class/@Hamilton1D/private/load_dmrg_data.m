function [S,t]=load_dmrg_data(HAM,k,varargin)
% function [S,t]=load_dmrg_data(HAM,k [,vars])
%
%   optional 2nd argument specifies time stamp
%   if DMRG data file if stored on file (-1 otherwise).
%
% see also @Hamilton1D/private/setupStorage.m
% Wb,Apr10,14

  if nargin<2 && ~isnumber(k), wbdie('invalid usage'); end
  if k<0 || k>numel(HAM.mpo)
     wbdie('invalid usage (k out of bounds: %g/%g)',k,numel(HAM.mpo));
  end

  tflag=0;
  for i=1:numel(varargin)
     if isequal(varargin{i},'-t'), tflag=1; varargin(i)=[]; break; end
  end

  if using_mem(HAM), s=HAM.store;
     S=struct; % initialize! // INIT_VARS_EVAL `i ml' 

     if k==0
          eval(sprintf('global %s_info; S=%s_info;',s,s));
     else eval(sprintf('global %s; S=%s(%g);',s,s,k));
     end

     if ~isempty(varargin)
        if numel(varargin)==1
             S=getfield (S,varargin{:});
        else S=getfields(S,varargin{:});
        end
     end

     if nargout>1, t=-1; end

  else mat=HAM.mat;

     if k==0
          f=sprintf('%s_info.mat',mat);
     else f=sprintf('%s_%03g.mat',mat,k); end

     if ~exist(f,'file')
        if tflag, S=[]; return; end
        wbdie('asking for non-existing data file %s',f);
     end

     S=save_load(f,varargin{:});

     if nargout>1
        t=getfield(dir(f),'datenum');
     end

     if numel(varargin)==1
        S=getfield (S,varargin{:});
     end

  end

end

% -------------------------------------------------------------------- %
% NB! if calculation is still running, file that is currently
% written, may not be readable here (message "got corrupt file")
% => pause, and retry. // Wb,Sep26,15

function S=save_load(f,varargin)

  n=4;
  for i=1:n
     try
        S=load(f,varargin{:});
        return
     catch
        l=lasterror; dispstack(l);
        if i<n, nsec=(2^i);
           fprintf(1,['\n ERR got error loading dmrg data ' ... 
             '(%g/%g; waiting for %g sec):\n\n> %s\n'],...
              i,n,nsec,regexprep(l.message,'\n','\n> '));
           pause(nsec);
        else rethrow(l); end
     end
  end

end

% -------------------------------------------------------------------- %
