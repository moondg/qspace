function info2(A,varargin)
% function info2(A,[dim,opts])
% Options
% 
%   dim  if dimension is specified, plot will be generated
%   -ll  log-log scale on plot
% 
% show further size info // Wb,Sep17,11

  info(A)

  if nargin>1 && isnumber(varargin{1})
     dim=varargin{1}; varargin(1)=[];

     s=getDimQS(A);
     if size(s,1)>1
        fprintf(1,'  %-32s: %s => %s states (%.3g)\n','number of multiplets',...
        int2str2(s(1,dim)), int2str2(s(2,dim)), s(2,dim)/s(1,dim));
     else
        fprintf(1,'  %-32s: %g\n','number of states',s(1,dim));
     end
  else dim=0; end

  getopt('init',varargin);
     llflag=getopt('-ll');
  getopt('check_error');

  if dim
     [qq,dd,dc]=getQDimQS(A,dim); dd=dd';
     fprintf(1,'  %-32s: %g\n','number of symmetry sectors',size(qq,1));
  end

  cgflag=gotCGS(A);
  if dim && cgflag
     s='multiplet space dimension';
     d=sum(repmat(dd,1,size(dc,2)) .* dc, 1)/sum(dd); fprintf(1,'\n');
     fprintf(1,'  %-32s: %s\n',[s,' (avg)'],vec2str(d,'fmt','%6.1f'));
     fprintf(1,'  %-32s: %s\n',[s,' (max)'],vec2str(max(dc),'fmt','%6g'));
  end

  fprintf(1,'\n');
  s=whos('A'); sa=s.bytes;
     fprintf(1,'  %-32s: %s\n','actual QSpace size',num2str2(sa,'-b'));

if ~cgflag, fprintf(1,'\n'); return; end

  cc=A.info.cgr; s=whos('cc'); sc0=s.bytes;

  [n,m]=size(cc);
  for i=1:n
  for j=1:m
     sc(i,j,1)=max(1,cc(i,j).nnz);
     sc(i,j,2)=prod(cc(i,j).size);
  end
  end

  sz=getzdim(A);

  scf=8*sum(reshape(prod(sz,3),[],1));
  scF=prod(sz,2); r=length(size(scF)); scF=permute(scF,[1,3:r,2]);
  saf=sa-sc0+scf;

  fprintf(1,'  %-32s: %-12s (factor %s)\n',...
    'size if full cgr data', num2str2(saf,'-b'), int2str2(round(saf/sa)));

  sd=datasize(A);
  sA=8*sum(prod(sd.*scF,2));

  fprintf(1,'  %-32s: %-12s (factor %s; lower estimate)\n',...
    'size if full tensor(data,cgr)', num2str2(sA,'-b'), int2str2(round(sA/sa)));

  fprintf(1,'\n  %-32s: %-12s\n',...
    'size sparse cgr / size QSpace', sprintf('%.3g',sc0/sa));

  s=8*[ sum(sum(sc(:,:,1))),  sum(sum(sc(:,:,2))) ];
  fprintf(1,'  %-32s: %-.3g   (%s [%s] => %s)\n', 'CGS sparsity', ...
     s(1)/s(2), num2str2(s(1),'-b'), num2str2(sc0,'-b'), ...
     num2str2(s(2),'-b') ...
  );

  fprintf(1,'\n');

if dim<1, return; end

ah=smaxis(1,1,'tag',mfilename,'fpos',[1330 520 590 590]);
header('%M'); addt2fig Wb
setax(ah(1,1))

  if 0
     xx=prod([dd dc],2); yy=dd;

     plot(xx,yy,'o'); hold on
  else
     plot(prod(dc,2),dd,'o','Color',[.8 .8 .8],'Disp','prod(cgr)'); hold on
     for i=1:size(dc,2)
        plot(dc(:,i),dd,'o','Color',getcolor(i),'Disp',sprintf('dim=%g',i));
     end
     legdisp('erase');
  end

  if llflag
     set(gca,'YScale','log');
     x=xlim; if x(2)>20, set(gca,'XScale','log'); end

     xtight(1.1,'x1',1.00); ytight(1.1,'y1',1.00);
     x=[xlim; ylim]; x=[min(x(:,1)),max(x(:,2))];
     set(gca,'XLim',x,'YLim',x');

     x={get(gca,'XTick'), get(gca,'YTick')};
     if numel(x{1})>numel(x{2}), xx=x{1}; set(gca,'YTick',xx);
     elseif numel(x{1})<numel(x{2}), xx=x{2}; set(gca,'XTick',xx);
     else xx=x{1}; end
     xx=logspace(0,10,21);
     x=xlim; y=ylim; x1=x(1)*1.1; y1=y(1)*1.1;
     ph=linspace(0,0.5*pi); o={'Color',[.9 .9 .9],'LineW',2};
     for i=1:numel(xx)
        h=plot([x1, xx([i i])], [xx([i i]),y1],o{:});
        mv2back(h);
     end
  end
  sms(4);

  label('dim(cgr)','dim(multiplet space)')
  title2('dim=%g got %g multiplet spaces (%s)',dim,size(qq,1),A.info.qtype);

  i0=find(prod(dc,2)==1); m=numel(i0); if m
     fprintf(1,['  got %g multiplet%s with scalar cg-space:  ' ...
       '%% [ #, Q, size(data,dim) ]\n'], m, iff(m>0,'s','') );

     for i=1:m
        s=sprintf(' %3g',qq(i0(i),:));
        fprintf(1,'   %2g. [%s ]  %4g\n',i,s,dd(i0(i)));
     end
     fprintf(1,'\n');
  end

end

