function I=showdeg(H,varargin)
% function [I=]showdeg(H [,opts])
%
%    show degeneracy of given rank-2 operator, either in 
%    diagonal representation, or to be diagonalized.
%
% Options
%
%   'dE',...  considere values degenerate within dE (1E-3)
%   'ng',...  number of degenerate groups to display
%   '-f'      force flag
%   '-v'      always print lowest part of spectrum
%   '-V'      more verbose still.
%
% Wb,Apr12,11

  if nargin<1 || numel(H.Q)~=2 || norm(H.Q{1}-H.Q{2})>1E-12
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('INIT',varargin);
     dE=getopt('dE',1E-3);
     ng =getopt('ng',[]);

     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end

     fflag =getopt('-f');
     ESflag=getopt('-ES');
  getopt('check_error');

  if ~isdiag(H), H_=H;
     [e,I]=eigQS(H);
     U=I.AK; H=QSpace(I.EK);

     if ESflag
        Emin=min(e(:,1)); n=numel(H.data);
        for i=1:n
           H.data{i}=-log(H.data{i});
        end
     end
  end

  n=numel(H.data); dd=cell(n,1);

  for i=1:n, m=length(H.data{i});
     dd{i}=[reshape(H.data{i},m,1), repmat(i,m,1),(1:m)',repmat(H.Q{1}(i,:),m,1)];
  end

  dd=cat(1,dd{:});
  [ee,i]=sort(dd(:,1)); iiq=dd(i,2:end);

  if nargout
     I=add2struct('-',H,'U?',ee,iiq);
     if isempty(ng) && ~vflag, return, end
  end

  n=length(ee); ii=[0; find(diff(ee)>dE); n]';
  if isempty(ng), ng=ii(min(end,4)); else ng=min(ng,n); end

  if ii(2)>12 && ~fflag, wblog('WRN',['got density matrix? ' ...
    '(%g-fold degeneracy!?)\nhint: use -ES flag or -f'],ii(2));
     return
  end

  if vflag>1
     for i=1:ng
        fprintf(1,'\n');
        for j=ii(i)+1:ii(i+1);
           fprintf(1,'   %2g  %2g  %8.5f [%s ]  @ %3g.%g\n',i,j,ee(j),...
           vec2str(iiq(j,3:end),'fmt','%3g'), iiq(j,1:2) );
        end
     end
   else
     for i=1:ng, jj=ii(i)+1:ii(i+1);
        if norm(diff(ee(jj)))>5*dE, ee(jj), norm(diff(ee(jj))),
           error('Wb:ERR','invalid degeneracy'); end
        fprintf('%s\n',repmat('-',50,1));

        [x,is]=sortrows(iiq(jj,3:end));
        iiq(jj,:)=iiq(jj(is),:);

        for j=jj
           if j==jj(1)
                e=mean(ee(jj));
                s={sprintf('  %2g.',i), sprintf('%10.6f  x%g',e,numel(jj)) };
           else s={'     ', ''}; end

           fprintf(1,'%s  [%s ]  @ %3g:%g %s\n',s{1},...
           vec2str(iiq(j,3:end),'fmt','%3g'), iiq(j,1:2),s{2} );
        end
     end
     fprintf('%s\n',repmat('-',50,1));
   end
end

