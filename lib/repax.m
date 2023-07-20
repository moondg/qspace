function ah = repax(varargin)
% function ah = repax(ah1,ah2)
%
%     replace axis set a1 by a copy of axis set a2
%
% Wb,Dec12,01

  nflag=0; ah1=[]; ah2=[];

  if nargin<2, helpthis
     if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  if nargin==2
     ah1=varargin{1};
     ah2=varargin{2};
  elseif nargin==3 && ~isempty(inputname(1)), nflag=1;
     ah1=varargin{1}(varargin{2});
     ah2=varargin{1}(varargin{3});
  end

  if ~isaxis(ah1) || ~isaxis(ah2)
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  set(gca,'selected','off')

  ph=get(ah1,'parent');
  set(ah1,'Units','normalized');

  save_disp(ah2);

  ah=copyobj(ah2,ph);

  set(ah,'position', get(ah1,'position'), 'Units', 'Normalized');
  delete(ah1)

  if numel(findall(gca,'Type','Line'))>512, drawnow; end

  restore_user_disp(ah2);
  restore_user_disp(ah);

  setax(ah);

  if nflag && ~nargout
     h=varargin{1}; h(varargin{2})=ah;
     assignin('caller',inputname(1),h);
  end

end

% -------------------------------------------------------------------- %

function save_disp(ah)
  for h=findall(ah,'type','line')'
     if isfield(get(h),'DisplayName')
        u.u=get(h,'UserData');
        u.d=get(h,'DisplayName'); set(h,'UserData',u);
     end
  end
end

% -------------------------------------------------------------------- %

function restore_user_disp(ah)
  hh=findall(ah,'type','line'); nh=0; tic;
  for h=hh'
     if isfield(get(h),'DisplayName'), u=get(h,'UserData');
        if ~isempty(u.d), set(h,'DisplayName',u.d); nh=nh+1; end
        set(h,'UserData',u.u);
     end
  end

  t2=toc; if t2>20
     a=get(hh(1),'Parent'); h=get(a,'Legend'); s={'',''};
     if ~isempty(h), s{1}=get(h,'AutoUpdate'); end
     if isequal(a,gca), s{2}='gca'; else s{2}='other axis'; end
     wblog('WRN',['restore_user_disp() took %.1f sec\n' ... 
       'having %d/%d non-empty Disp; AutoUpdate=%s in %s'],...
       t2,nh,numel(hh),s{:});
  end
end

% -------------------------------------------------------------------- %

