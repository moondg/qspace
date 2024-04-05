function [oo,aa,ah,I] = rsmoothSpec(Om,Aa,varargin)
% Usage: [oo,aa,ah,I] = rsmoothSpec(Om,Aa [,opts])
%
%   subsequent smoothening of spectral function
%
% Options:
%
%   Ifdm    if arg[3] is structure, T will be read out of it (if exists)
%   'T',..  temperature (for plot only)
%   'adisp',{..}  set of labels used as Disp in aa data
%   'afac'  factor when plotting smooth spectral function (pi*Gamma)
%   'nfac'  extra factor (1/nfac) when using SU(N) symmetries
%   'skip'  skip smooth data for |om|<skip (0)
%
%   'nofig' do not show data in figure
%   'xli'   limits for x axis on inset
%   'yli'   y-limits for main panels
%   'raw'   input data Aa is considered raw data
%   'RAW'   same as 'raw', but also plots raw data
%
% NB! remaining options are handed over to getSmoothSpec(Om,Aa,....)
%
% Wb,May12,06

  global N param

  if nargin<2, eval(['help ' mfilename]); return; end

  Ifdm=[];
  if length(varargin)>0
     if isstruct(varargin{1})
     Ifdm=varargin{1}; varargin=varargin(2:end); end
  end

  if ~isempty(param)
     TK=TKondo; xli=5*TK;
     if isfield(param,'B') && param.B>TK, xli=5*param.B; end
     if isequal(xli,0), xli=[]; end
  else
     xli=2E-4;
  end

  getopt ('INIT', varargin);
     T     = getopt('T',     []);
     eps   = getopt('eps',   []);
     reps  = getopt('reps',  []);
     nofig = getopt('nofig'    );
     afac  = getopt('afac',  []);
     adisp = getopt('adisp', {});
     nfac  = getopt('nfac',   1);
     xli   = getopt('xli',  xli);
     yl    = getopt('yli',   []);
     skip  = getopt('skip',   0);
     sigma = getopt('sigma', -1);
     alpha = getopt('alpha', []);
     keepA0= getopt('keepA0', 0);

     if getopt('raw'), rflag=1; 
     elseif getopt('RAW'); rflag=2; else rflag=0; end

  varargin=getopt('get_remaining');

  if length(xli)==1, xli=[-1 1]*xli; end

  if isempty(T) && isfield(Ifdm,'T') && ~isnan(Ifdm.T)
  T=Ifdm.T; end

  if sigma<=0 && ~isempty(param) && isfield(param,'Lambda')
     if param.Lambda==2, sigma=0.6; else
        sigma = log(param.Lambda) * (0.6/log(2));
     end
  end

  if isempty(eps) && ~isempty(T) && T>0
     if ~isempty(reps), eps=T*reps; else eps=T; end
  end
  if ~isempty(eps), varargin(end+1:end+2)={'eps',eps}; end

  if keepA0, o={'-a0'}; else o={}; end
  if sigma>0, setopts(o,sigma); end
  if ~isempty(alpha), setopts(o,alpha); end

  if nfac~=1, Aa=Aa*(1/nfac); end

  [oo,aa,I] = getSmoothSpec(Om,Aa,varargin{:},o{:});

  if skip>0
     i=find(abs(oo)<skip);
     oo(i)=[]; aa(i,:)=[];
     I.skip=skip;
  end

  if nofig, ah=I; return; end

  if isempty(afac)
     if ~isempty(param)
        if isfield(param,'fJ') || isfield(param,'J')
           afac=(pi^2)/2; facstr='\pi^2/2';
        elseif isfield(param,'Gamma')
           afac=pi*max(abs(param.Gamma(:))); facstr='\pi\Gamma';
        else afac=1; facstr=''; end
     else
        wblog('WRN','param not available - scale max_A to 1.');
        afac=1/max(abs(aa(:))); facstr=sprintf('%.4g',afac);
     end
  else
     if iscell(afac)
        if numel(afac)~=2 || ~isnumeric(afac{1}) || ~ischar(afac{2})
        afac, wbdie('invalid afac'); end

        facstr=afac{2}; afac=afac{1};
     elseif afac~=1
        facstr=sprintf('%.4g',afac);
     else facstr=''; end
  end

  I.afac={afac,facstr};
  I.nfac=nfac;

  if rflag
     a=sum(Aa); a=a(find(abs(a)>0.5 & abs(a-round(a))<1E-3));
     if ~isempty(a), a=real(a(1)); end
  end

ah=smaxis(2,1,'tag',mfilename,'DY',0.02);
setax(ah(1));

  ah(1,2)=inset({'NE',[0 -.05]},'scale',[1.2 1]); % see tag 'Spec01i' below

  if ~isempty(param), TK=TKondo; else TK=nan; end

  if isempty(yl)
     yl=[-0.05,  1.2]*afac*max(real(aa(:)));
     if diff(yl)<=0, yl=[]; end
  end

  lfmt0 = {'color', [.7 .7 .7],'tag','raw'};
  lfmtp = {'color', [.7 .7  1],'tag','raw'};
  lfmtn = {'color', [ 1 .7 .7],'tag','raw'};
  if rflag, lg1='cc raw data / \omega'; else lg1='A0'; end
  if rflag, lg2='cc raw data / \omega'; else lg2='A0'; end

  astr='A(\omega)';
  if ~iscell(adisp) && ischar(adisp)
     switch adisp
        case {'spin_ud','spin'}
        adisp={'A_{\uparrow}(\omega)','A_{\downarrow}(\omega)'};
        case 'spin_du'
        adisp={'A_{\downarrow}(\omega)','A_{\uparrow}(\omega)'};
     otherwise; end
  end
  if isempty(adisp)
     n=size(aa,2); adisp=cell(1,n);
     for i=1:n
         adisp{i}=sprintf('A_{%g}(\\omega)',i);
     end
  elseif ~iscell(adisp)
  adisp, wbdie('invalid adisp'); end

  ox=Om(:); ox(find(ox==0))=nan; ip=find(ox>0); in=find(ox<0);

  oa={'tag','Adata'};

  if isreal(aa)
     q=norm(imag(aa(:)))/norm(real(aa(:)));
     if q>1E-8
        wblog('-->','got complex spectral data @ %.3g',q); 
     end
  end
  ar=real(aa);
  cflag = (norm(imag(aa))/norm(aa) > 1E-3);
  copt={':','LineWidth',1.5};

setax(ah(1))

  h2=plot(oo,afac*ar,oa{:});
  if cflag, coloridx --reset
     hi=plot(oo,afac*imag(aa),copt{:},oa{:});
     set(hi(1),'Disp','imag. part');
  else hold on; end

  if ~isempty(yl), ylim(yl), end

     n=min(numel(adisp),numel(h2));
     for i=1:n, set(h2(i),'Disp',adisp{i}); end

  if rflag>1
     Ar=afac*Aa; or=Om;
     Ar=Ar./repmat(abs(ox),1,size(Aa,2));

     nl=log(10)/abs(mean(diff(log(abs(ox(2:5))))));
     m=round(nl/250);
     if m>1
        l=mod(size(Ar,1),m);
        if l, l=m-l;
           Ar(end+(1:l),:)=0;
           or(end+(1:l))=or(end);
        end
        s=size(Ar); Ar=permute(sum(reshape(Ar,[m,s(1)/m,s(2)]),1),[2 3 1]);
        s=length(or)/m; or=permute(sum(reshape(or,[m,s]),1),[2 1]);
     end
     Ar=real(Ar);

     h1=plot(or,Ar,lfmt0{:});

     if ~isempty(I.a0), y=ylim;
     plot([0 0], [0 y(2)],'-', lfmt0{:},'LineW',2); end

     mv2front(h2);

     for i=1:length(h1)
        set(h1(i),'Color',0.7+0.2*getcolor(i));
     end
     set(h1(1),'Disp',lg1);
  else 
  end

  if ~isempty(afac) && isempty(a==1), a=[a 1]; end

  xlim([-1 1]*0.2);
  ymark([0 a],'k--');

  if ~isempty(facstr), astr=[astr '  \ast ' facstr]; end

  label('\omega',astr);
  legdisp({'SE',[0 0.03]},'-flip'); 

setax(ah(3))

  h=plot(oo,afac*ar,oa{:});

    n=min(numel(adisp),numel(h));
    for i=1:n, set(h(i),'Disp',adisp{i}); end

  mv2front(h(1)); hold on

  if diff(xli)>0, xlim(xli); else
     if ~isempty(xli) && ~all(isnan(xli))
     wblog('WRN','invalid xli=[%s]',vec2str(xli)); end
     xtight
  end

  ytight('view'); l=ylim; ytight(1.2,'view'); l2=ylim; l2(1)=l(1); ylim(l2);
  ymark([0 a],'k--');

  if ~isempty(param) && isfield(param,'B') && param.B~=0
  xmark([-1 1]*param.B,'k--'); end

  ok={'Color',[.5 .75 .5],'LineW',2}; % ,'tag','TK'

  if isset('TK'), xmark([-1 1]*TK,ok{:}); end

setax(ah(2))

  if rflag>1
     jp=find(or>0); jn=find(or<0);
     h2=semilogx(abs(or(jn)),Ar(jn,:), lfmtn{:}); coloridx --reset
     h1=semilogx(abs(or(jp)),Ar(jp,:), lfmtp{:}); coloridx --reset
     set(h1(1),'Disp',lg1);
  else
  end

  ip=find(oo>0);
  in=find(oo<0);

  hn=semilogx(-oo(in), afac*ar(in,:),oa{:}); coloridx --reset
  h2=semilogx( oo(ip), afac*ar(ip,:),oa{:});
  set(hn,'LineSt','--'); set(hn(1),'Disp','\omega<0 data');

  if cflag, coloridx --reset
     hn_=semilogx(-oo(in), afac*imag(aa(in,:)),copt{:},oa{:}); coloridx --reset
     h2_=semilogx( oo(ip), afac*imag(aa(ip,:)),copt{:},oa{:}); blurl(hn_);
     set(h2_(1),'Disp','imag. data');
  end

    n=min(numel(adisp),numel(h2));
    for i=1:n, set(h2(i),'Disp',adisp{i}); end

  mv2front(h2); hold on; xtight

  if ~isempty(yl), ylim(yl); else
  ylim(get(ah(1),'YLim')); end

  topt = { 'HorizontalAlignment','Center', ...
  'VerticalAlignment','middle' };

  if ~isempty(param) && isfield(param,'Lambda') && ~isempty(N)
     xl=param.Lambda^(-N/2);
     lh=plot([xl xl], ylim, 'LineWidth',4, 'color', [1 .95 .5]);
     text(xl,.8, sprintf('\\Lambda=%g, N=%g', param.Lambda, N), topt{:});
     mv2back(lh);
  end

  if ~isempty(param)
     xmark(TK,'top','istr','T_K',ok{:});
  end

  if isfield(I,'eps')
     c2=[1 1 1]*0.8;
     xmark(I.eps,'top','istr','\epsilon','Color',c2,'LineW',3);
  end

  if ~isempty(T)
  xmark(T,'top','istr','T'); end
  ymark([0 a],'k:');

  s={};
  if isfield(I,'sigma')
     s{end+1}=sprintf('\\sigma=%.3g',I.sigma);
  end
  if isset('alpha')
     s{end+1}=sprintf('\\alpha=%.3g',alpha);
  end

  if isfield(param,'ALambda') && ~isempty(param.ALambda)
     s{end+1}='using A_\Lambda';
  end
  if ~isempty(s)
     s=sprintf('(%s)',strhcat(s{:},'sep','; '));
     text(.88,.70-0.05*numel(adisp),s,'Units','norm','HorizontalAl','center');
  end

  if ~isempty(I.a0)
     for i=1:2:size(I.a0,2)
         s=sum(sum(abs(I.a0(:,i:i+1))));
         if s>1E-10, j=(i+1)/2;
            o={'Color',get(h2(j),'Color')};
            postext(.04,1-i*.1,'A_%d: %.6g \\delta(\\omega)',j,s,o);
         end
     end
  end

  label('|\omega|',astr);
  legdisp -flip

  set(ah(1),'tag','Spec01');
  set(ah(2),'tag','Spec02');
  set(ah(3),'tag','Spec01i');

  if ~nargout, clear oo aa; end

end

