function Qs=getQ0str(H,varargin)
% function Qs=getQ0str(H)
% Wb,May30,11

  getopt('init',varargin);
     istr=getopt('istr','Q');
     eps =getopt('eps',0.01);
  getopt('check_error');

  if isfield(H,'HK') && isfield(H,'E0')
     Inrg=H;
     l=numel(Inrg.E0)-2;
     e0=Inrg.E0(l:l+1); if e0(1)>e0(2), l=l+1; end
     H=Inrg.HK(l);
  end

  H=QSpace(H);
  if isdiag(QSpace(H),'-d')<2
     [ee,H]=eig(H);
  end

  qf=getqfmt(H,'sep',', ');

  dd=H.data; e0=min(cat(2,dd{:}))+eps; n=numel(dd); Qs={};
  for i=1:n
     if ~isempty(find(dd{i}<=e0,1))
        Qs{end+1}=['(' sprintf(qf,H.Q{1}(i,:)) ')'];
     end
  end

  if numel(Qs)>1
       Qs=sprintf(',%s',Qs{:}); Qs=sprintf('\\{%s\\}',Qs(2:end));
  else Qs=Qs{1}; end

  if ~isempty(istr), Qs=sprintf('%s=%s',istr,Qs); end

end

