function dd=getSymXOM(sym)
% function dd=getSymXOM(varargin)
%    
%    generate statistics on XStore directory
%    with respect to largest outer multiplicity in x3 tensors.
%
% Wb,Jan30,15

  D0=get_dir(sym,'XStore');

  dd=dir(D0);
  dd=dd(find([dd.isdir])); nD=numel(dd);

  rmax=3;
  for id=1:nD, D=dd(id).name; if D(1)=='+'
     r=length(D); if rmax<r, rmax=r; end
  end, end
  rmax=rmax+1;
  SS={};

  for id=1:nD, D=dd(id).name; if D(1)~='.'
     fprintf(1,'\r   reading XStore/%s ... %10s\r',D,''); D=[D0 '/' D '/'];
     ff=dir(D);
     ff=ff(find([ff.isdir]==0));
     [~,is]=sort([ff.bytes]); ff=ff(is);

     nf=numel(ff);
     om=zeros(nf,1); r=zeros(nf,1);
     ss=zeros(nf,rmax);
     xs=zeros(nf,5);

     for k=1:nf, l=3*(k-1);
        q=load([D ff(k).name],'-mat');
        ss(l+1,:)=get_size_OM(q.x3m.a,rmax);
        ss(l+2,:)=get_size_OM(q.x3m.b,rmax);
        ss(l+3,:)=get_size_OM(q.x3m.c,rmax);

        q=get_size_x3(q.x3m.x3);
        xs(k,:)=[id,k,q];
     end

     SS{id}=uniquerows(ss);
     XS{id}=xs;
     FF{id}=ff;

  end, end

  ss=uniquerows(cat(1,SS{:}));
  [m,I,dm]=uniquerows(ss(:,1));

  mD=[m';dm]; if mD(1)==0, mD=mD(:,2:end); end
  fprintf(1,'\r%60s\n','');
  fprintf(1,'   statistics on outer multiplicity for symmetry %s:\n\n',sym);
  fprintf(1,'     OM    # entries (roughly, ignoring duals)\n',sym);
  fprintf(1,'     %2g  %8g\n',mD);
  fprintf(1,'\n');

  xs=uniquerows(cat(1,XS{:}));
  [m,I,dx]=uniquerows(sort(xs(:,3:end),2,'descend'));

  q=xs(I{end}(1),:);
  fprintf(1,'   largest x3 tensor obtained\n');
  fprintf(1,'     %s/%s @ %gx%gx%g\n\n', dd(q(1)).name, FF{q(1)}(q(2)).name,q(3:end));

end

% -------------------------------------------------------------------- %

function s=get_size_OM(a,rmax)
  if isempty(a), s=zeros(1,rmax); return; end

  r=length(a.qdir); s=double(a.cgd.S); l=length(s);
  if l<r || l>r+1, error('Wb:ERR','\n   ERR invalid cgd.size !?'); end
  if l==r, m=1; else m=s(end); end

  s=[m,zeros(1,rmax-r-1),s(1:r)];
end

% -------------------------------------------------------------------- %

function s=get_size_x3(x3)
  if isempty(x3.S), s=zeros(1,3); else s=x3.S; end
end

% -------------------------------------------------------------------- %

