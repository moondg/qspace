function s=num2str2(val,varargin)
% Usage #1: function s=num2str2(val, lim1,scale1,tag1, lim2,tag2,scale2, ...)
%
%    returns sprintf('%.3g %s',val/lim_i,tag1) if val<lim_i
%
% Usage #2: s=num2str2('-n',val1,str1,val2,str2, ...)
%
%    eg. num2str2('-n',1,'error',0,'warning',...)
%    returns: 1 error and no warnings (i.e. this checks plural s)
%    Wb,Mar10,11
%
% Usage #3: s=num2str2(val,eps)
%
%    e.g. num2str2(0.45093,0.009)
%    returns: 0.450(9)
%    where y=a(b) is to be read as follows: y=a with the uncertainty
%    b on the the last significant digit(s) shown with a comes
%    i.e., 0.450(9)
%        = 0.450 +/- '9' on the last significant digit here
%        = 0.450 +/- 0.009
%        = 0.4500(90) etc.
%    By default, (b) contains a single digit if b > ~ 3
%    and two digits if b is close to 1 for b>=1.
%
% Wb,Apr17,09

% https://chemistry.stackexchange.com/questions/58661/..
%    what-does-a-number-in-brackets-after-another-number-mean-i-e-211cm-1
% tags: bytes2str, byte2str, sizestr (meaning size=bytes)

  getopt('init',varargin);
     vflag=getopt('-v');
     fmt  =getopt('fmt',[]);

     if getopt('-b'), bflag=1;
     elseif getopt('--bytes'), bflag=2;
     else bflag=0; end

  varargin=getopt('get_remaining'); narg=length(varargin);

  if ~narg && numel(val)==1 && ~bflag
     if isempty(fmt), fmt='%.4g'; end
     if isreal(val)
        s=sprintf(fmt,val);
        if isempty(findstr(s,'.'))
           s=regexprep(s,'[eE][+0]*','E');
        else
           q=regexprep(s,'.*\.\d[eE][+0]*([1-9]+)','$1');
           if ~isequal(q,s), q=str2num(q)-1;
              s=sprintf([fmt 'E%g'],val*10^(-q),q);
           end
        end
     else
        r=real(val); i=imag(val);
        if r==0, s=sprintf([fmt 'i'],i);
        else
           s={ sprintf(fmt,r), sprintf(regexprep(fmt,'%\+*','%+'),i) };
           s{2}=regexprep(s{2},'^\s*([\+\-])',' $1 ');
           s=[s{:} 'i'];
        end
     end
     return
  end

% -------------------------------------------------------------------- %
  if isequal(val,'-n'), nset=1;
  elseif isequal(val,'-n0'), nset=2; else nset=0; end

  if nset
     if mod(narg,2), wberr('invalid usage (%g)',narg); end
     s={};
     for i=1:2:narg
        q=varargin{i}; t=varargin{i+1};
        if numel(q)~=1 || ~isnumeric(q), wberr(...
           'invalid usage (alternating numbers expected)'); end
        if ~ischar(t), wberr(...
           'invalid usage (alternating string expected)'); end
        if q~=0 || nset>1
           if q==0 s{end+1}=sprintf('no %ss',t);
           elseif q==1, s{end+1}=sprintf('%g %s',q,t);
           else s{end+1}=sprintf('%g %ss',q,t); end
        end
     end
     if numel(s)>1
        s(2,:)={', '}; s{2,end-1}=' and '; s{2,end}='';
     end

     s=cat(2,s{:}); return
  end

% -------------------------------------------------------------------- %
  if ~bflag && narg && isnumeric(varargin{1}) && isscalar(varargin{1})
     dval=varargin{1};

     if dval<=0
        if ~dval, s=sprintf('%g(0)',val);
        else wblog('ERR','got negative dval=%g',dval);
             dval=sprintf('%g',val);
        end; return
     elseif val==round(val)
        s=sprintf('%d(%d)',val,dval); return
     elseif dval>=abs(val)
        s=sprintf('%.3g +/- %.3g',val,dval); return
     end

     x=floor(log10(abs(val/dval)));

     if 1
        if isempty(fmt)
           if x>4, fmt=sprintf('%%.%dg',x+2);
           else fmt='%g'; end
        end

        s=sprintf(fmt,val);
        i=regexp(s,'[eE]');
        ss=strread(s,'%s','whitespace','eE'); lfac=0;

        if ~isempty(i), ss={s(1:i-1), s(i), s(i+1:end)};
             lfac=-str2num(s{3});
        else ss={s,'',''}; end

        i=regexp(ss{1},'\.'); n=length(ss{1});
        if ~isempty(i)
            n=n-i;
            lfac=lfac+n;
        else i=regexp(ss{1},'\d+'); n=n-i(1); end

        q=round(dval*(10^lfac));
        if q>999 && n>2, q=round(q/10); ss{1}=ss{1}(1:end-1); n=n-1; end
        if q> 99 && n>2, q=round(q/10); ss{1}=ss{1}(1:end-1); n=n-1; end
        if q> 29 && n>2, q=round(q/10); ss{1}=ss{1}(1:end-1); n=n-1; end

        ss{1}=[ss{1} '(' sprintf('%d',q) ')'];
        s=[ss{:}];
     else
        if isempty(fmt)
           if x>4, fmt=sprintf('%%.%dg',x+2);
           else fmt='%.16g'; end
        end
        ss={ sprintf(fmt,val)
             sprintf(fmt,val+dval)
             sprintf(fmt,val-dval) };
        s=sum(diff(double(strvcat(ss)),1).^2,1);
        i=regexp(ss{1},'[eE]'); n=numel(s);
        if ~isempty(i)
           if numel(i)~=1 || any(s(i:end))
              wblog('ERR','got no matching exponent !?');
              s=sprintf('[%g %g]',val,dval);
              return;
           end
        else i=n+1;
        end

        j=find(s,1);
        if isempty(j) || j>=i
           wblog('ERR','failed to represent accuracy of number');
           s=sprintf('<ERR %g!?>',val); return
        end

        if vflag
           b=repmat(' ',1,n); b(j)='|'; if i<=n, b(i)=':'; end
           fprintf(1,'\n'); fprintf(1,'  %s\n',ss{:},b);
           fprintf(1,'\n');
        end

        s=ss{1};
        s=[s(1:j-1) '(' s(j) ')' s(i:end)];
     end
     return
  end

% -------------------------------------------------------------------- %
  if bflag
     if ~narg
        varargin={ 5E3, 1,    ''
                   1E6, 2^10, 'k'
                   1E9, 2^20, 'M'
                  1E12, 2^30, 'G'
                  1E15, 2^40, 'T' };
        if bflag>1, varargin(:,3)={' bytes',' kB',' MB',' GB',' TB'}; end
        varargin=reshape(varargin',1,[]);
        narg=numel(varargin);
     else wberr('invalid usage'); end
  end

  if narg<3 || mod(narg,3) || ~isscalar(val) || ~isnumeric(val)
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  varargin=reshape(varargin,3,[]); n=size(varargin,2);
  for i=1:n
     if ~isscalar(varargin{1,i}) || ~isnumeric(varargin{1,i}) || ...
        ~isscalar(varargin{2,i}) || ~isnumeric(varargin{2,i}) || ...
        ~ischar(varargin{3,i}), wberr('invalid usage');
     end
  end

  [x,i]=sort(cat(2,varargin{1,:})); varargin=varargin(:,i);

  if isempty(fmt), fmt='%.4g%s';
  elseif isempty(regexp(fmt,'%.*s')), fmt=[ fmt ' %s']; end

  for i=1:n
     if val<varargin{1,i} || i==n
        s=sprintf(fmt,val/varargin{2,i},varargin{3,i});
        break;
     end
  end

end

