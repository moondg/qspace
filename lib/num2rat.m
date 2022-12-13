function [s,q]=num2rat(xx,varargin)
% Function:[s,i] = num2rat(x [,fmt, sep])
%
%    convert number to rational string format using wbrat.mex
%    if rational approximation fails within default accuracy,
%    float number is printed using fmt (default: %.4g)
%
% Options
%
%    -r         include ..sqrt(..) formats
%   --rtex      include ..\surd(..) formats
%    'pmax',..  for x=p/q, abs([p,q]) <= pmax at least for p or q
%
% Wb,Dec29,09 ; Wb,Sep17,20

% for simple version based on matlab's rat
% see Archive/num2rat_110112.m

  getopt('init',varargin);
     if     getopt('-r'),     rflag=1;
     elseif getopt('--rtex'), rflag=2;
     else rflag=0; end
     pmax=getopt('pmax',99);
  args=getopt('get_remaining'); nargs=length(args);

  if nargs>2
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  if nargs && ~isempty(args{end}), sep=args{end};
  elseif isvector(xx), sep=', '; else sep=' '; end

  if nargs>1 && ~isempty(args{1}), fmt=args{1}; else fmt='%.4g'; end

  xx(find(abs(xx)<1E-12))=0;
  [z,q]=wbrat(xx,'-q');

  ss=cell(size(xx)); n=numel(xx);

  if rflag
     if rflag>1, sqrt_='\surd'; l=3;
     else sqrt_='sqrt'; l=6; end
  end

  for i=1:n, x=xx(i); sp=sprintf(fmt,x);
     if q.ee(i)==0 && any(abs([q.P(i),q.Q(i)])<pmax)
        if q.Q(i)==1
             ss{i}=sprintf(fmt,q.P(i));
        else ss{i}=sprintf([fmt '/%g'],q.P(i),q.Q(i));
        end
        if length(ss{i})>length(sp), ss{i}=sp; end
     elseif rflag
        [z(i),q1]=wbrat(x^2,'-q');
        if q1.ee==0
           if x<0, s='-'; elseif n>1, s=' '; else s=''; end
           if q1.Q==1
                ss{i}=[s sqrt_ sprintf('(%g)',q1.P)];
           else ss{i}=[s sqrt_ sprintf('(%g/%g)',q1.P,q1.Q)];
           end
           if length(ss{i})>l+length(sp), ss{i}=sp; end
        else ss{i}=sp; end
     else ss{i}=sp; end
   end

   for i=1:n, ss{i}=[sep ss{i}]; end

   if isvector(ss)
      s=cat(2,ss{:}); s=s(numel(sep)+1:end);
   else
      s=zeros(size(ss));
      for i=1:n, s(i)=numel(ss{i}); end;  n=max(s(:));
      for i=1:n, ss{i}(end+1:n)=' '; end
      for i=1:size(ss,1), ss{i,1}=cat(2,ss{i,:}); end
      s=cat(1,ss{:,1});
   end

end

