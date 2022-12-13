function h=addfinfo(f)
% Usage:  h=addfinfo([f = global savemat])
% Wb,Dec18,12

  f=[]; H=hostname;

  getbase Imain

  if isempty(Imain), q0='---xvar-tmp---';
     setuser(groot,'xvar',q0);
     evalin('caller','if exist(''Imain'',''var''), setuser(groot,''xvar'',Imain); end');
     q=getuser(groot,'xvar','-rm'); if ~isequal(q,q0), Imain=q; end
  end

  if ~nargin, getbase savemat mat
     if ~isempty(savemat), f=savemat;
     elseif ~isempty(mat), f=mat;
     elseif isfield(Imain,'mat'), f=Imain.mat; end
  end

  if ~isempty(f), setuser(gcf,'mat',f);
     i=find(f=='/',1,'last'); if ~isempty(i), f=f(i+1:end); end
  else f=''; end

  if isfield(Imain,'host') && ~isempty(Imain.host)
     H_=Imain.host; if iscell(H_), H_=H_{1}; end
     h=H_; t='';
     try
        dt=datediff(Imain.finished, Imain.started);
        if all(~isnan(dt)) && ~isempty(dt)
           h=[h ' // ' dt];
        else
           t=evalc(['! date -d ''' Imain.started ''' ''+%y%m%d''']);
           if isempty(regexp(f,t))
           h=[h ' // ' dt]; end
        end
        if ~isequal(H_,H) || ~isempty(f), h=[H 10 h]; end
        if ~isempty(f), f=[ h 10 f]; else f=h; end
     catch l
        wblog('WRN'); f, t
        disp(l.message); dispstack(l);
     end

  elseif ~isempty(f), f=[ H 10 f];
  else f=H; end

  if ~isempty(f), f=strrep(f,'_','\_');
    h=header('fright', f);
    set(h,'FontSize',7,'VerticalAlign','bottom','Pos',[0.995 0.001]);
  else
    h=[];
  end

  if ~nargout, clear h; end

end

