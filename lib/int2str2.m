function ss=int2str2(dd,varargin)
% function ss=int2str2(dd)
%
%    convert integer to string,
%    e.g. introducing integer comma separator for thousand
%
% Options
%
%    '-th'     using 1st, 2nd, 3rd, etc.
%    '-TH'     same as '-th', but using superscript tex format: 1^{st},...
%
%    '-x'      convert integers to single digit chars
%              assuming range [0,35]
%    '-X'      same as -x, but case sensitive [0,61]
%    '-A'      similar to -X, but only use alphabet [0,51]
%   '--ext',.. extend beyond -X still appending the chars provided
%
% Examples
%
%    ss=int2str2(0:103,'--ext','_~'); % 103 = 26*4-1
%
% Wb,Aug09,10

  xflag=0;
  getopt('INIT',varargin);
     if     getopt('-th'), sflag=1;
     elseif getopt('-TH'), sflag=2;
     else sflag=0;
        xstr=getopt('--ext','');
        if     getopt('-x'), M=[ 10  26   0 ];
        elseif getopt('-X'), M=[ 10  26  26 ];
        elseif getopt('-A') || ~isempty(xstr), M=[  0  26  26 ];
        else M=[]; end

        xflag=sum(M)-1; xflag1=0;
        if ~isempty(xstr)
           if ~ischar(xstr), wbdie('invalid usage'); end
           xflag1=xflag+1; xflag=xflag1*numel(xstr)-1;
        end
     end

  getopt('check_error');

  if xflag>0
     e=norm(dd-round(dd));
     if e, wbdie('invalid usage (got non-integers !?)'); else e={}; end
     if any(dd<0), e{end+1}=sprintf('int out of bounds (%g)',min(dd)); end
     if any(dd>xflag)
        e{end+1}=sprintf('int out of bounds (%g/%g)',max(dd),xflag); end
     if ~isempty(e)
        e(2,:)={char(10)}; e(end)=[];
        if xflag1
             wbdie([e{:}]);
        else wbwrn([e{:}]); end
     end

     if xflag1, k=floor(dd/xflag1)+1; dd=mod(dd,xflag1); end
     q=[0, cumsum(M)]; c0={'0','a','A'};
     ss=dd;

     for l=2:numel(q)
        i=find(dd>=q(l-1) & dd<q(l));
        ss(i)=dd(i)-q(l-1)+c0{l-1};
     end
     ss=char(ss);

     if ~isempty(xstr)
        ss=reshape(ss,[],1);
        ss(:,2)=xstr(k);
     end

     return
  end

  nd=numel(dd);

  if sflag
     if nd~=1 || dd<0, wbdie('invalid usage'); end
     if     dd==1, ss='st';
     elseif dd==2, ss='nd';
     elseif dd==3, ss='rd'; else ss='th'; end

     if sflag>1
          ss=sprintf('%g^{%s}',dd,ss);
     else ss=sprintf('%g%s',dd,ss);
     end
     return
  end

  ss=cell(size(dd));
  for i=1:nd
      e=dd(i); if e, e=abs((e-round(e))/e); end
      if e>1E-12
         wblog('WRN','data appears to contain non-integer (%g)',dd(i));
         ss{i}=sprintf('%f',dd(i)); continue;
      end

      s=sprintf('%.0f',dd(i)); n=length(s);
      if n>3
         m=mod(n,3); if m, m=3-m; s=[repmat(' ',1,m),s]; end
         s=reshape(s,3,[]); s=[repmat(',',1,size(s,2)); s];
         s=s(m+2:end);
      end

      ss{i}=s;
  end
  if numel(ss)==1, ss=ss{1}; end

end

