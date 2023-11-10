function info3(A,dim,varargin)
% function info3(A,dim [,opts])
% graphical analysis of size distribution of CGC spaces

% example usage:
if 0
% U(1) SU(2) SU(3) symmetry
  sl={{'U(1) charge'}, {'SU(2) spin','2S_z'}, {'SU(3) channel','q_1','q_2'}};
  info3(A3,3,'sfac',0.005,'-x','syml',sl);

  sl={{'SU(2) spin','2S_z'},...
      {'SU(2) channel 1','2C_z'},{'SU(2) channel 2','2C_z'},{'SU(2) channel 3','2C_z'}};
  info3(A3,3,'sfac',0.004,'-x','syml',sl);
end

  if ~gotCGS(A)
      wblog('ERR','QSpace with CGC spaces required');
      return
  end
  if nargin<2 || numel(A)~=1
     helpthis, if nargin || nargout
     wbdie('invalid usage'), end, return
  end

  getopt('init',varargin);
     sfac=getopt('sfac',1);
     xflag=getopt('-x');
     syml=getopt('syml',{});
     if getopt('-i'), info2(A,dim); end
  getopt('check_error');

  [qq,D]=getQDimQS(A,dim);
  [s,ss]=getqdim(A); iq=[0, cumsum(s)]; ss0=ss;
  if ~isempty(syml), ss(1:numel(syml))=syml; end

  cgr=A.info.cgr;
  ns=size(cgr,2); m=ceil(sqrt(ns));

  y1fac=1;

if isequal(s,[1 1 2])
  ah=smaxis(3,1,'tag',mfilename,'fpos',[1095 350 500 770]);
  ah=[ splitax(ah(1),[2,1]); mergeax(ah(2:3))];
  mvaxis(ah(2),[0 -0.02]);
  mvaxis(ah(3),[0 -0.02]);
  y1fac=0.15;
elseif isequal(s,[1 1 1 1])
  ah=smaxis(4,1,'tag',mfilename,'fpos',[990 350 606 770]);
  y1fac=0.22;
else
  ah=smaxis(m,m,'tag',mfilename)';
end

  set(ah,'LineW',1.5); addt2fig Wb

  ph=linspace(0,2*pi,32);
  xc=cos(ph);
  yc=sin(ph);
  zc=zeros(size(ph));

  for k=1:ns, setax(ah(k)); hold on; view(2)
     [q,ik,dk]=uniquerows(qq(:,iq(k)+1:iq(k+1)));
     [dk,is]=sort(dk); q=q(is,:); ik=ik(is);

     rq=size(q,2);
     if rq==1, q(:,2)=0;
     elseif rq>2
        wbdie('don''t know how to plot rank-%g object',rq); 
     end

     m=size(q,1); d=zeros(1,m); dd=zeros(m,2);
     for i=1:m
         di=D(ik{i}); d(i)=sum(di);
         x=q(i,1);
         y=q(i,2); rfac=sfac;

         if xflag
            for j=1:numel(di)
               r=rfac*sqrt(di(j));
               plot3(x+r*xc,y+r*yc,-di(j)+zc,'-','LineW',1);
            end
            view(2)
         else
            r=rfac*sqrt(d(i));
            plot3(x+r*xc,y+r*yc,zc,'-','LineW',1);
         end
         dd(i,:)=[max(di), d(i)];
     end
     DD(k,:)=[ max(dd(:,1)), sum(dd(:,2))];

     if iscell(ss{k})
        ts=ss{k}{1};
     else ts=ss{k}; end

     if xflag
        title(sprintf('%s',ts));
     else
        title(sprintf('%s: d_q\\leq%g [%g]',ts,max(dd(:,2)), max(dd(:,1)) ));
     end

     if iscell(ss{k}) && numel(ss{k})>1, xlabel(ss{k}{2}); end
     if iscell(ss{k}) && numel(ss{k})>2, ylabel(ss{k}{3}); end

     xtight(1.1);
     if rq>1; axis equal,
     else
        set(gca,'YTick',[]); axis equal
        set(gca,'YLim',y1fac*[-.5 +.5]*diff(get(gca,'XLim')));
     end
     set(gca,'XTick',unique(round(get(gca,'XTick'))));
     set(gca,'YTick',unique(round(get(gca,'YTick'))));
  end

  set(ah(k+1:end),'Visible','off');

  if norm(diff(DD(:,2))), DD(:,2)'
     wbdie('inconsistency in total multiplet dimension');
  end
  header('%M [d_q\leq%g; D^\ast=%g]',max(DD(:,1)),unique(DD(:,2)));

end

