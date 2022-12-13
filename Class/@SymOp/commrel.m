function IC=commrel(G,varargin)
% function IC=commrel(G [,opts])

%    get commutator relations of all elements in input symmetry G
%    this also checks / requires completeness of group G, and that
%    the commutator relations are canoncial. The resulting roots
%    of the symmetry are returned in the output structure IC.
%
% Options
%
%   -v   verbose flag
%
% Wb,Nov28,11

% adapted from $PHYS/commrel.m for input of cell-array
% see also $PHYS/SU3_commops.m

  getopt('init',varargin);
     vflag=getopt('-v');
  getopt('check_error');

  gg=getops(G);

  g2=getops(G,'-catV');

  [u,s,v]=svd(g2,'econ'); s=diag(s);
  k=find(s<1E-12); if ~isempty(k), error('Wb:ERR',['\n   ' ...
   'ERR input group is not linearly independent (%g)'],numel(k)); end

  n=numel(gg); xx=zeros(n);
  for i=1:n
  for j=i:n
     xx(i,j)=olap(gg{i},gg{j});
  end,end

  e=norm(triu(xx,1)); if e>1E-12
     wblog('WRN','got non-orthogonal operator space!'); end
  ng2=diag(xx);

  iz=idxzops(G);

  f3=zeros(n,n,n); ee=zeros(n,n); err=0; wrn=0; ss={}; io=0; nz=0;
  for i=1:n
  for j=i+1:n
     q=comm(gg{i},gg{j}); v=reshape(full(q),[],1);
     w=g2\v; r=g2*w-v; ee(i,j)=norm(r); l=find(abs(w)>1E-12); m=numel(l);

     if ee(i,j)>1E-12, err=err+1;
        wblog('WRN','CR leaves group (%g,%g @ %g)',i,j,ee(i,j));
     elseif m>1
        if ~isempty(setdiff(l,iz)) %|| norm(full(gg{i}-gg{j}'))>1E-12
           wrn=wrn+1; wblog('WRN',...
          'CR not in canonical form (%g,%g => %s)',i,j,vec2str(l));
        elseif vflag, ss{end+1}={i,j,l,w(l)}; end
     elseif m==1, if vflag, io=io+1;
        fprintf(1,['%5g.  [ G(%2g), G(%2g) ] = %4.6g G(%2g)  ' ...
          '[ %-12s, %-12s ] => %s\n'], io, i, j, w(l), l, ...
           G(i).istr, G(j).istr, G(l).istr);
        end
     else, nz=nz+1; if vflag,
        fprintf(1,'%7s [ G(%2g), G(%2g) ] = [ %s, %s ] = 0\n',...
        '',i,j,G(i).istr,G(j).istr); end
     end

     f3(i,j,:)= w;
     f3(j,i,:)=-w;
  end
  end

  if ~isempty(ss), fprintf(1,'\n');
  for is=1:numel(ss), [i,j,l,w]=deal(ss{is}{:}); io=io+1;

     fprintf(1,'%5g.  [ G(%2g), G(%2g) ] =',io,i,j);
     for k=1:numel(l), 
        s=sprintf('%+.6g', w(k));
        fprintf(1,' %s G(%g)',[s(1) ' ' s(2:end)],l(k));
     end
     fprintf(1,'\n');

     fprintf(1,'        [ %12s , %12s ] =>', G(i).istr, G(j).istr);
     for k=1:numel(l), fprintf(1,'  %s', G(l(k)).istr); end

     fprintf(1,'\n\n');
  end, end

  if vflag
     wblog(' * ','got %g non-trivial commutator relations',io);
     wblog(' * ','got %g commutators resulting in 0',nz);
  end

  f3(find(abs(f3)<1E-12))=0;
  gm=contract(f3,f3,[2 3],[3 2]);

     gm(find(abs(gm)<1E-12))=0;
     e=norm(imag(gm)); if e<1E-12, gm=real(gm);
     else wblog('WRN','got complex metric (@%g)',e); end

  g=eig(gm);
  if any(abs(g)<1E-12), err=err+1;
     wblog('ERR','metric is singular!');
  elseif err==0, wblog('SUC',...
    'got proper symmetry group with %g elements (@ %.1g)',n,norm(ee));
  end

  r=numel(iz);
  rr=cell(1,r); err=0;
  for l=1:r
     q=permute(f3(iz(l),:,:),[2 3 1]); rr{l}=diag(q);
     e=norm(q-diag(rr{l})); if e>1E-12, err=err+1;
       wblog('WRN','got invalid roots');
     end
  end

  if ~err
     rr=cat(2,rr{:}); e=norm(rr-round(rr)); if e<1E-12, rr=round(rr);
     else wblog('WRN','roots consist of non-integers!'); end

     is=find(sum(rr.^2,2)>1E-12);
     [ia,ib,i]=matchIndex(rr(is,:),-rr(is,:));
     if ~isempty(i.ix1) || ~isempty(i.ix2) || numel(ia)~=numel(is) || ...
        any(ia==ib), wblog('ERR','failed to match roots');
     else
        i=find(ia<ib); i2=[is(ia(i)), is(ib(i))]; n=size(i2,1);
        for i=1:n
           [x,j]=sortrows(fliplr(rr(i2(i,:),:)));
           if j(1)==1, i2(i,:)=i2(i,[2 1]); end
        end
        R=struct('iz',iz,'ip',i2(:,1),'im',i2(:,2),...
          'rz',rr(iz,:), 'rp',rr(i2(:,1),:), 'rm',rr(i2(:,2),:));
        if size(uniquerows(R.rp),1)==size(R.rp,1) ...
        && size(uniquerows(R.rm),1)==size(R.rm,1)
             wblog(' * ','got unique set of roots'); 
        else wblog('WRN','got degenerate roots !??'); end

        zz=zeros(n,r);
        for i=1:size(i2,1)
           f=squeeze(f3(i2(i,1),i2(i,2),:));
           zz(i,:)=f(iz); f(iz)=0; if norm(e)>1E-12
              wblog('WRN','[Sp,Sp''] leaves space of z-ops !??');
           end
        end
        R.zz=zz;
     end
  end

  IC=add2struct('-',f3,gm,g,rr,R);

end

