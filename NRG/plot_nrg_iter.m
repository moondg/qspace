function ee=plot_nrg_iter(varargin)
% function ee=plot_nrg_iter(varargin)
% 
%    plot HK and HT data for a given iteration.
%
% Usage
%
%   (1) plot_nrg_iter(q);          struct with entries q.HK and q.HT 
%   (2) plot_nrg_iter(HK,HT);      HK,HT as specified
%   (3) plot_nrg_iter([NRG,]k);    load data from file
% 
% Wb,Sep01,12

  usage=[]; istr='';

  if nargin && isstruct(varargin{1}) && ...
     isfield(varargin{1},'HK') && isfield(varargin{1},'HT'), usage=1;
     varargin={ varargin{1}.HK, varargin{1}.HT, varargin{2:end} };
  elseif nargin && isstruct(varargin{1}) && ...
     isfield(varargin{1},'EK') && isfield(varargin{1},'ET'), usage=1;
     varargin={ varargin{1}.EK, varargin{1}.ET, varargin{2:end} };
  elseif nargin==2 && ...
     isQSpace(varargin{1}) && isQSpace(varargin{2}), usage=2;
  elseif nargin && nargin<=2 && isnumber(varargin{end}), usage=3;
     if nargin==1
        varargin={ [getenv('LMA') '/NRG/NRG'], varargin{:} };
     end
     try 
        if ischar(varargin{1})
           f=sprintf('%s_%02g.mat',varargin{1},varargin{2});
           q=load(f,'HK','HT');
           varargin={ q.HK, q.HT, varargin{3:end} };
        else usage=[]; end
     catch
        wblog('ERR','invalid NRG data (%s)',f);
        return
     end
     istr=untex(repHome(f));
  end

  if isempty(usage) || numel(varargin)>2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
setax(ah(1,1)); box on

  k=0; kk=[]; dy=[]; yn=nan; emin=inf;
  nq=numel(varargin); DD=cell(1,nq); it={};

  for i=1:nq
     data=varargin{i}.data;
     DD{i}=sort(cat(2,data{:})); 
     if ~isempty(DD{i}), emin=min(emin,min(DD{i})); end

     if isfield(varargin{i}.info,'itags')
        it{i}=varargin{i}.info.itags;
     end
  end

  if ~isempty(it)
     it=reshape([it{:}],numel(it{1}),[])'; s=sprintf(', %s',it{:,1});
     header('%M :: \{%s\}',s(3:end));
  end

  for i=1:nq
     dd=DD{i}; n=numel(dd); if emin, dd=dd-emin; end
     oc={'color',getcolor(i)};
     ot={oc{:},'FontSize',8};

     h=plot(k+1:k+n, dd,'o-',oc{:}); mv2back(h); hold on
     if i==1, set(h,'Disp','kept');
     elseif i==nq, set(h,'Disp','discarded'); end

     if n && i<nq
        h=ymark(dd(end),'k:',oc{:}); mv2back(h);
        text(n,dd(end),sprintf('   ^{\\leftarrow}  %.4g',dd(end)),...
          'VerticalAl','top',ot{:}); % 'E_n=...'
        yn=dd(end);
     end
     if n && i>1
        h=ymark(dd(1),'k:',oc{:}); mv2back(h);
        text(k+1,dd(1),sprintf('%.4g  _{\\rightarrow}   ',dd(1)),...
          'HorizontalAl','right',ot{:}); % 'E_1=...'
        dy(i-1)=dd(1)-yn;
     end
     kk(end+1)=k+0.5; k=k+n;
     if nargout, ee{i}=dd; end
  end

  label('state index s','E_s');
  legdisp('Location','NW','erase');

  ymark -fix

  h=findall(gca,'type','text','tag','ymark'); n=numel(h);
  x=xlim; dx=diff(x)/n; x=dx/10:dx:x(2);
  for i=1:n
      p=get(h(i),'Pos'); p(1)=x(n-i+1); p(3)=10; set(h(i),'Pos',p);
  end
  set(h,'FontW','normal','FontSize',10,'Margin',1E-3,'BackgroundC','w');

  sms(4); xtight;
  if numel(kk)>1
     title2('Nkeep=[%s]',vec2str(kk(2:end)-0.5));
     xmark(kk(2:end),'k:');
  end
  if ~isempty(istr), title(istr); end

  if emin, postext('SE','NB! subtracted E_0=%.4g',emin); end

end

