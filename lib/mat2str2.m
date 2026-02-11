function [sout,Iout] = mat2str2(M, varargin)
% function s = mat2str2(M [,opts])
%
%    Writes matrix M as string
%    also using simple representation of complex values where possible,
%    like 1i instead of matlab's 0.0000 + 1.0000i, etc.
%
% Options
%
%   'fmt',..   is the format string ('%8g')
%   'sep',..   is the separator string (' ')
%   'rowsep',..is the string to separate rows ('\n')
%   'istr',..  info/intro string
%   '-f'       no shortcuts (enforce full mode)
%   'eps',..   skip tiny numbers on the numerical noise level
%              (default: auto; use eps=0 to turn off)
%              
%   'phase'    abs|phase instead of real+imag
%   'nofac'    do not use overall factors pulled to the front (1 for scalar, 0 otherwise)
%
%   '-c'       return as cell array of strings, where the returned
%              cell array has exactly the same dimensions as M
%
%   '-p'       print result (default, if no output argument is requested)
%   '-p','..'  use string in '..' as name for the matrix shown
%
% See also existing MatLab routine mat2str().
% Wb,Jul11,03

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end
  sout='';

  getopt('init',varargin);
     fmt    = getopt('fmt',[]);
     sep    = getopt('sep',' ');
     rowsep = getopt('rowsep','\n');
     istr   = getopt('istr','');
     deps   = getopt('eps',-1); % former 'notiny' // Wb,Oct28,25
     phflag = getopt('phase');
     fflag  = getopt('-f');
     cflag  = getopt('-c');
     pflag  = getopt('-p');
     nofac  = getopt('nofac');
     vn=getopt('get_last','');
  getopt('check_error');

  if nargout>1, Iout=add2struct('-',fmt); end

  if cflag
     n=numel(M); if isempty(fmt), fmt='%g'; end
     M=reshape(mat2cell(M(:),ones(n,1),1),size(M));
     for i=1:n, M{i}=sprintf(fmt,M{i}); end
     sout=M; return
  elseif isempty(fmt), fmt='%8g';
  end

  s=size(M); r=numel(s); n1=s(1); n2=s(2); n=numel(M);
  if ~nofac && isscalar(M), if isreal(M), nofac=1; end; end

  if n==1
     sout=sprintf(fmt,M); return
  end
  if ~fflag
     if all(diff(M(:))==0)
        sout=sprintf([fmt ' (%s)'],M(1),vec2str(s,'sep','x')); return
     end
  end

  done=0;

  if r>2, wbdie('invalid usage (got rank-%g object)',r); end
  if ~fflag && n1==n2
     d=diag(M);
     if norm(M-diag(d))==0
        if all(diff(d)==0)
             sout=sprintf([fmt ' (eye; %gx%g)'],M(1),s);
        else sout=sprintf('diag([%s])',vec2str(d,'fmt',fmt)); end
        done=1;
     end
  end

if ~done
  if deps, a=abs(max(M(:)));
     if a>1e6*eps
        if deps<0, deps=a*1e6*eps; end
        if isreal(M), M(abs(M)<deps)=0;
        else
           i=find(abs(imag(M))<deps); M(i)=real(M(i));
           i=find(abs(real(M))<deps); M(i)=imag(M(i));
        end
     end
  end

  if nofac
     fac=1;
  else
     fac = max (abs(M(:)));
     if fac~=0, fac = floor(log10(fac)); end
     if strfind(fmt, 'd'), fac=1; end
     if abs(fac) > 3
        fac = 10^fac;
        M = M / fac;
     else
        fac = 1;
     end
  end

  if nargout>1
     Iout=add2struct(Iout,fmt,deps,fac,sep,rowsep);
  end

  sout=cell(n1,n2);

  if isreal(M)
     for i=1:n1
       for j=1:n2
          sout{i,j}=sprintf(fmt,M(i,j));
       end
     end
  else
      [fm1,r] =  strtok(fmt, '%.gefGEF');
      [fm2,r] =  strtok(r,   '%.gefGEF');
      [fmc,r] =  strtok(fmt, '%+-.0123456789');

      pdot = strfind(fmt,'.');
      if ~isempty(fm1) && ~isempty(pdot)
          if pdot<strfind(fmt,fm1)
             fm2 = fm1;
             fm1 = '';
          end
      end

      if  ~isempty(fm1), fm1 = str2num(fm1); else fm1=8; end

      if ~isempty(fm2)
         fm2 =str2num(fm2);
         fmtr=sprintf('%%.%dg', fm2);
         fmtc=sprintf('%%+.%dg', fm2);
      else
         fm2 = [];
         fmtr='%g';
         fmtc='%+g';
      end
	  fmts = sprintf('%%%ds', fm1);

      for i=1:n1
        for j=1:n2, mij=M(i,j);
          if real(mij)==0
             if imag(mij)==0
                  vstr = '0 '; % align with 1i etc. // previously '0.'
             else vstr = sprintf([fmtr 'i'], imag(mij));
             end
          else
             if imag(mij)==0
                vstr = sprintf(fmtr, real(mij));
             elseif ~phflag
                  vstr = sprintf([fmtr fmtc 'i'], real(mij), imag(mij));
             else vstr = sprintf([fmtr '|' fmtr], abs(mij), angle(mij)/pi);
             end
          end

          sout{i,j}=sprintf(fmts,vstr);
       end
     end
  end

  rowsep=sprintf(rowsep);
  sep=sprintf(sep);

  for i=1:n1
     for j=1:n2, s=sout{i,j};
        if abs(M(i,j))==0 && length(s)>=2 && isequal(s(end-1:end),'-0')
           s(end-1:end)=' 0';
        end
        if     j>1, s = [    sep s];
        elseif i>1, s = [ rowsep s];
        end
        sout{i,j}=s;
     end
  end

  sout=sout'; sout=[sout{:}];

  if ~isempty(istr) && isempty(find(istr=='='))
  istr=[istr ' = ']; end

  if fac~=1 | ~isempty(istr)
     if fac~=1, vstr = sprintf('%1.0E * ', fac);
     else       vstr = ''; end

     if size(M,1)>1 && (~isempty(strfind(rowsep,'\n')) | ~isempty(strfind(rowsep,10)))
         if     ~isempty(find(istr=='[')), bs='\n]';
         elseif ~isempty(find(istr=='{')), bs='\n}'; else bs=''; end
         sout = sprintf(['%s%s\n\n%s' bs], istr, vstr, sout);
     else
         sout = sprintf('%s%s[%s]',   istr, vstr, sout);
     end
  end
end

  if ~nargout || pflag
     if isempty(vn), vn=inputname(1); % if isempty(vn), vn='ans'; end
     elseif ~ischar(vn), vn, wbdie('invalid usage (variable name)');
     end

     if ~isempty(vn)
          fprintf(1,'\n   %s = \n\n',vn);
     else fprintf(1,'\n'); end
     if done==1, fprintf(1,'      '); end

     disp(sout); fprintf(1,'\n');

     if ~nargout, clear sout; end
  end

end

