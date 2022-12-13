function addto_dmrg_info(HAM,varargin)
% function addto_dmrg_info(HAM [,vars])
% see also @Hamilton1D/private/save_dmrg_info.m
% Wb,Apr14,14

  if nargin<2, wberr('invalid usage'); end

  nvars={};
  for i=1:numel(varargin), nvars{i}=inputname(i+1);
     if isempty(nvars{i})
        if iscell(varargin{i}) && numel(varargin{i})==2
           nvars{i}=varargin{i}{1};
           if ~ischar(nvars{i}) || ~isvarname(nvars{i})
               wberr(['invalid usage ' ... 
              '(name for input must be accessible or specified)']);
           end
           varargin{i}=varargin{i}{2};
        end
     end
  end

  if ~isempty(HAM.store), s=HAM.store;
     S=struct; % initialize! // INIT_VARS_EVAL `i ml' 

     eval(sprintf('global %s_info; S=%s_info;',s,s));

     for i=1:numel(varargin)
        S=setfield(S,nvars{i},varargin{i});
     end

     eval(sprintf('%s_info=S;',s));

  elseif ~isempty(HAM.mat), mat=HAM.mat;

     for i=1:numel(varargin),
        eval(sprintf('%s=varargin{%g};',nvars{i},i));
     end
     save(sprintf('%s_info.mat',mat),nvars{:},'-append');

  else wberr('invalid storage specification'); end

end

