function plot_op(S,H)
% function plot_op(S,H)
%
%    S operator (possibly IROP)
%    H Hamiltonian (in diagonal representation) => show energy range
%
% Wb,Aug20,12

  [i,i1,I]=matchIndex(S.Q{1},H.Q{1},'-s'); if ~isempty(I.ix1)
     wblog('ERR','got non-match Q-spaces (op/H; 1)'); end

  [i,i2,I]=matchIndex(S.Q{2},H.Q{2},'-s'); if ~isempty(I.ix1)
     wblog('ERR','got non-match Q-spaces (op/H; 2)'); end

  r=numel(S.Q);
  if r>3, wblog('ERR','got rank-%d object !??',r); end

  Sfac=0; nd=numel(S.data); 

  try
     Sfac=trace(S);
     if abs(Sfac)<1E-12, wblog('WRN',...
        'trace() too small (%.3g) => use max(data{:}) instead',Sfac);
        Sfac=max_data(S.data);
     end
  catch
     wblog('WRN','failed to calculate trace => use max(data{:}) instead');
     Sfac=max_data(S.data);
  end
  S=S*(1/Sfac);

ah=smaxis(2,2,'tag',mfilename,'fpos',[1210 450 710 670]);
addt2fig Wb

  w=[0.16, 0.4, 0.04];
  x=[0.14, 0.12+w(1)];
  set(ah(1,1),'Pos',[x+[1,-1]*(w(1)+w(3)), w([2 1])]);
  set(ah(1,2),'Pos',[x+w(1:2)+w(3), w([2 1])]);
  set(ah(2,1),'Pos',[x, w([1,2])]);
  set(ah(2,2),'Pos',[x+[w(1)+w(3), 0], w([2 2])]);
  set(ah(2:end),'FontSize',8);

  ee=cat(2,H.data{:}); er=[ min(ee), max(ee) ]; q=[];

  i=1:nd; 
  while keyiter(i), di=S.data{i};
     setax(ah(2,2))
        s=size(di); ri=numel(s);
          if ri>2, wblog('ERR','got rank-%d object !??',ri); end
        mp(di,'-gca','cmap');
        header('fright',get(get(gca,'title'),'String'));
        label('','','');
     setax(ah(2,1))
        ee=H.data{i1(i)}; n=numel(ee); plot(ee,0.5:n,'.-'); axis tight
        set(gca,'XDir','reverse','YDir','reverse');
        xlabel('E'); set(gca,'YLim',get(ah(2,2),'YLim')); xlim(er);
        if n~=size(di,1), ylabel(['\bf\color{red}', sprintf(...
          'length mismatch in S/H data (Y: %d/%d)',n,size(di,1)) ]);
        end
     setax(ah(1,2))
        ee=H.data{i2(i)}; n=numel(ee); plot(0.5:n,ee,'.-'); axis tight
        set(gca,'YAxisLoc','right');
        ylabel('E'); set(gca,'XLim',get(ah(2,2),'XLim')); ylim(er);
        if n~=size(di,2), title(['\bf\color{red}', sprintf(...
          'length mismatch in S/H data (X: %d/%d)',n,size(di,2))]);
        end

     for j=1:r, q(j,:)=S.Q{j}(i,:); end
     header('%M :: %2d/%d :: [%s] [%s] (Sfac=%.4g)',i,nd,...
       mat2str2(q,'rowsep','; ','fmt','%g'), ...
       mat2str2([H.Q{1}(i1(i),:), H.Q{2}(i2(i),:)],'rowsep','; ','fmt','%g'),...
       Sfac);

     setax(ah(1,1))
        if isequal(q(1,:),q(2,:))
           h=plot(diag(di),'o-'); sms(h,4); xtight('x1',0);
           ylabel(sprintf('diag(data\\{%d\\})',i))
        else cla; title(''); end

     set(ah(2:end),'FontSize',8);

  end

end

function q=max_data(data);

  n=numel(data); q=zeros(1,n);
  for i=1:n, q(i)=max(abs(data{i}(:))); end
  q=max(q);

end

