function plotSS(HAM,varargin)
% function plotSS(HAM [,opts])
%
% Examples
%
%    plotSS(HAM,wsys,{'S5',SX,'dk',3}); % having wsys='HBLX'
%
% Wb,Jan24,18

  if nargin>1 && iscell(varargin{end})
       opts=varargin{end}; varargin(end)=[];
  else opts={};
  end

  [s1,s2,I]=getSS(HAM,varargin{:});

  getopt('init',varargin);
     if getopt('HBLX'),
          plot_HBLX(HAM,s1,s2,I,opts{:});
     else 
        wbdie('invalid usage (need to specify system)');
     end
  getopt('check_error');

end

% -------------------------------------------------------------------- %
function plot_HBLX(HAM,sl,sr,I,varargin)

  getopt('init',varargin);
     sfac = getopt('sfac',1);
     tflag= getopt('-t');
     xfac = getopt('xfac',1);
     dk   = getopt('dk',[]);
     ah   = getopt('ah',[]);
     xl   = getopt('xl',[]);
     S5   = getopt('S5',[]);

     if getopt('--rbox')
          sh={'shape','rbox','--noedge'};
     elseif getopt('--rbox2')
          sh={'shape','rbox2','--noedge'};
     else sh={}; end

  varargin=getopt('get_remaining');
  if isempty(varargin), varargin={'-pn'}; end

  if isempty(HAM) || ~isfield(HAM.info,'HSS')
     wbdie('invalid usage (missing field HAM.info.HSS)'); end
  L=numel(HAM.mpo);

  ss={sl,sr,I.sx};
  S5flag=(~isempty(S5)); if S5flag, ss{end+1}=S5(2:end,:); end
  ns=numel(ss);

  if ~isempty(dk) || dk~=1
     if dk<1 || mod(dk,2)~=1, dk
        wbdie('invalid dk (odd value required)');
     end
     for i=1:ns, ss{i}=ss{i}(1:dk:end,:); end
     l=length(ss{2});
     for i=[1 3:ns], ss{i}(l:end,:)=[]; end
  end

  q=[]; for i=1:3, q=[q;ss{i}(:)]; end
  smean=mean(abs(q)); q=1/smean; for i=1:ns, ss{i}=q*ss{i}; end

if isempty(ah)
ah=smaxis(1,1,'tag',mfilename);
  header('%M :: %s',HAM.info.istr); addt2fig Wb
end

setax(ah(1))

  Sr=ss{2};
  for k=1:numel(Sr)

     x=[k k]; y=[0 1];
     if ~tflag
          plotbond(x,y,Sr(k),sh{:},varargin{:});
     else plot(x,y,'Color',getcolor(1+floor(100*rand)));
     end

     hold on
  end

  if S5flag
     Sl=ss{end}(:,2:3);
  else
     wblog('WRN','only got averaged <S.S>_leg data! (specify {''S5'',..} option)');
     Sl=repmat(ss{1},1,2);
  end

  for k=1:size(Sl,1)

     x=[k,k+1];
     for j=1:2, y=[j,j]-1;
        if ~tflag
             plotbond(x,y,Sl(k,j),sh{:},varargin{:});
        else plot(x,y,'Color',getcolor(1+floor(100*rand)));
        end
     end

     hold on
  end

  o={'EdgeC','none','FaceC',[1 .94 .79]};

  Sx=ss{3};

  for k=1:numel(Sx)

     x=[k k+1]; i=[1 2 2 1 1];
     y=[0 1  ]; j=[1 1 2 2 1];
     if ~tflag
        patch(x(i),y(j),repmat(-1,1,numel(i)),'y',o{1:end-1}, ...
        max(1-((1-o{end})*xfac)*Sx(k),0));
     else
        x=[k+.2,k+.8]; y=[.2,.8];
        plot(x(i),y(j),'Color',getcolor(1+floor(100*rand)));
     end

     hold on
  end

  title(sprintf('<S.S> \\sim %.4g',smean));
  if ~isempty(xl), xlim(xl); else xl=xlim; end

end

% -------------------------------------------------------------------- %

