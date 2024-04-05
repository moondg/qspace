function [xz,yz,Iout]=findpeak(varargin)
% function [xz,yz,Iout]=findpeak(xd,yd,...)
% function [xz,yz,Iout]=findpeak('-gca',...)
% function [xz,yz,Iout]=findpeak(lh,...)
%
% Options
%
%    '-p'     plot fits to extrema
%    'plot',. plot options
%    '--min'  look for minimas only
%    '--max'  look for maximas only
%    'xr',[x1 x2] limit x-range to search
%    'ymin',. min. absolute height of peak
%
%    'dk',[k1 dk] consider subset (k1:dk:end) of data only
%    'dk',dk  take data subset (*:dk:end) with largest values if --max
%             take data subset (*:dk:end) with smallest values if --min
%             else average data over dk data points
%    'n',..   order of polynomial to use for fitting (3)
%    '-gca'   use all line handles in current data set
%
% Wb,May24,17

  getopt('init',varargin);

     if     getopt('--min'), wp=-1;
     elseif getopt('--max'), wp=+1;
     else                    wp= 0; end

     if getopt('-p')
          popts={'default'};
     else popts=getopt('plot',{});
     end

     ymin =getopt('ymin',0);
     dk   =getopt('dk',[]);
     xr   =getopt('xr',[]);
     np   =getopt('n',3);

     if getopt('-gca'), ah=gca;
     else ah=getopt('ah',[]); end

  varargin=getopt('get_remaining');
  opts={xr,np,dk,ymin,wp,popts};

  if ~isempty(ah), t=mfilename; xz={}; yz={};
     hh=findall(ah,'type','line'); nh=numel(hh); mark=ones(1,nh);
     for i=1:numel(hh), h=hh(i);
      % do not include earlier fitted data by *this!
        if isequal(get(h,'Tag'),t), mark(i)=0;
        else
         % also exclude marker sets or contour lines
           xd=get(h,'XData');
           if numel(xd)<=2 || any(diff(xd)<=0), mark(i)=0; end
        end
     end
     hh=hh(find(mark)); nh=numel(hh); xz=cell(1,nh); yz=xz; pp=xz;

     for i=1:nh, h=hh(i);
        [xz{i},yz{i},pp{i}]=findpeak_1(h,[],opts{:},varargin{:});
     end
     if nh==1, xz=xz{1}; yz=yz{1}; end
  elseif ishandle(varargin{1})
     hh=varargin{1}; nh=numel(hh); varargin(1)=[];
     xz=cell(1,nh); yz=cell(1,nh);
     for i=1:nh, h=hh(i);
        [xz{i},yz{i},pp{i}]=findpeak_1(h,[],opts{:},varargin{:});
     end
     if nh==1, xz=xz{1}; yz=yz{1}; end
  else
     if numel(varargin)<2 || ~isnumeric(varargin{2})
          k=1; yd=varargin{k}; xd=1:length(yd);
     else k=2; xd=varargin{k-1}; yd=varargin{k};
     end
    
     [xz,yz,pp]=findpeak_1(xd,yd,opts{:},varargin{k+1:end});
  end

  if nargout
     if numel(pp)==1, pp=pp{1}; end
     Iout=add2struct('-',pp);
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [xz,yz,pp]=findpeak_1(xd,yd,xr,np,dk,ymin,wp,popts,varargin)

  if numel(varargin)
     wblog('WRN','ignoring additional options\N'); disp(varargin)
  end

  if isempty(yd) && ishandle(xd), h0=xd;
     if ~isline(h0), disp(h0), wbdie('got invalid line handle'); end
     xd=get(h0,'XData'); yd=get(h0,'YData');
  end

  n=numel(dk); xz=[]; yz=[]; pp=[];
  if n==1
     if dk>1 && dk==round(dk)
        if wp==0
           yd=avgdata(yd,dk,'-l');
        else
           nd=numel(yd); q=zeros(dk);
           q=sum(reshape(yd(1:nd-mod(nd,dk)),dk,[])',1);
           if wp>0
              % searching for maximum => keep interleaved set with largest values
                i=find(q==max(q),1);
           else i=find(q==min(q),1);
           end

           i=i:dk:nd; xd=xd(i); yd=yd(i);
        end
     elseif dk~=1, wblog('WRN','got dk=%g !?',dk); end
  elseif n==2
     if dk(1)<1 || dk(2)<1, dk, wbdie('invalid usage'); end
     i=dk(1):dk(2):numel(yd);
     xd=xd(i);
     yd=yd(i);
  elseif ~isempty(dk), dk, wbdie('invalid usage'); end

  pflag=~isempty(popts);
  if pflag
     ot={'tag',mfilename};
     om={{'Color',[1 0 1]} % color for maxima : 'm' (magenta)
         {'Color',[0 1 1]} % color for minima : 'c' (cyan)
     };

     mmark='*'; vflag=0; tfmt=''; bflag=0;
     to={'FontSize',8,'HorizontalAl','center','VerticalAl','bottom'};

     if ~isequal(popts,{'default'})
        getopt('INIT',popts);
           q=getopt('cm',[]);
           if ~isempty(q), om{1}{2}=q; om{2}{2}=q; 
           else
              om{1}{2}=getopt('cmax',om{1}{2});
              om{2}{2}=getopt('cmin',om{2}{2});
           end

           bflag=getopt('--bg');

           if getopt('-v'), vflag=vflag+1;
           elseif getopt('-V'), vflag=vflag+2; end

           if vflag, mmark='+'; end
           mmark=getopt('Marker',mmark);

           tfmt=getopt('tfmt',tfmt);

           q=getopt('topts',{});
           if ~isempty(q), to={to{:},q{:}}; end

        getopt('check_error');
     end

     if vflag
        if isempty(tfmt), f='%.3g';
           if vflag==1
                tfmt=['y=' f];
           else tfmt=['(' f ',' f ')'];
           end
        else
           i=numel(find(tfmt=='%'));
           if vflag>1 && ~isempty(regexp(tfmt,'^%[\w\.+-]*$'))
              tfmt=['(' tfmt ',' tfmt ')'];
           elseif vflag==1 && i~=1 || vflag>1 && i~=2
              wbdie('invalid tfmt=''%s'' having vflag=%g !?',tfmt,vflag); 
           end
        end
     end

     hold on
  end

  if ~isempty(xr)
     i=find(xd>=xr(1) & xd<=xr(2));
     if isempty(i)
        wblog('WRN','got empty search range [%g %g]',xr);
        return
     elseif numel(i)<10, wblog('WRN',...
       'got too few data points within search range [%g %g]',xr);
        return
     end
     xd=xd(i); yd=yd(i);
   % if pflag, plot(xd,yd,'m',ot{:}); end
  end

  nd=numel(xd);

  if size(xd,1)==1, xd=xd'; end
  if size(yd,1)==1, yd=yd'; end

  x1=avgdata(xd); n1=numel(x1);
  y1=diff(yd); s1=sign(y1);

% peak might be given by two points at exactly the same height
% (eg. when taking functional data!) // Wb,Mar17,18
  i=find(s1);
  if s1(1)==0, s1(1:i(1)-1)=s1(i(1)); end
  if s1(end)==0, s1(i(end)+1:end)=s1(i(end)); end
  while 1
     i=find(s1==0); i=i(find(s1(i-1))); s1(i)=s1(i-1);
     i=find(s1==0); i=i(find(s1(i+1))); s1(i)=s1(i+1);
     if isempty(i), break; end
  end

  if wp==0
       IZ=find(s1(1:end-1)~=s1(2:end)); % max or min
  elseif wp>0
       IZ=find(s1(1:end-1)> s1(2:end)); % max only
  else IZ=find(s1(1:end-1)< s1(2:end)); % min only
  end

  if ymin
     IZ=IZ(find(abs(yd(IZ))>=ymin | abs(yd(IZ+1))>=ymin));
  end

  iw = (sign(y1(IZ)) - sign(y1(IZ+1)))/4 + 1.5;

  if 1
   % take *more* data points than np+1
   % (otherwise data may become badly conditioned)
     np_=ceil((np+1)/2)*[1 1];
  else
     if mod(np,2)
          np_=(np+[ 1,+1])/2;
     else np_=(np+[ 0, 2])/2; end
  end

  nZ=numel(IZ); xz=zeros(1,nZ);

  for k=1:nZ
     i2=[IZ(k)-np_(1)+1, IZ(k)+np_(2)+1];
     if i2(1)<1, i2=[1, min(n1,1+np)];
     elseif i2(2)>n1, i2=[max(1,n1-np), n1]; end

     ir=i2(1):i2(2);

     x0=mean(x1(ir));
     p1=polyfit(x1(ir)-x0,y1(ir),np);
     z1=roots(p1)+x0;
     q=abs(z1-x1(IZ(k))); iz=find(q==min(q),1);
     z=real(z1(iz));

     x0=mean(xd(ir));
     p=polyfit(xd(ir)-x0,yd(ir),np);
     xz(k)=z;
     yz(k)=polyval(p,z-x0);

     if pflag
        xx=linspace(xd(ir(1)),xd(ir(end)),256);
        h1=plot(xx,polyval(p,xx-x0),'-',ot{:},om{iw(k)}{:});
        hm=plot(z,yz(k),mmark,ot{:},om{iw(k)}{:}); sms(hm,6);

        if bflag, mv2back([h1,hm]); end

        if vflag
           if vflag==1
                s=sprintf(tfmt,  yz(k));
           else s=sprintf(tfmt,z,yz(k));
           end
           text(z,yz(k)+0.05*diff(ylim),s,to{:},ot{:},om{iw(k)}{:});
        end
      % keyboard
     end
     p1(end)=p1(end)+x0;
     pp{end+1}=p1;
  end

% keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %


