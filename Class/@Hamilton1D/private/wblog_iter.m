function Il=wblog_iter(k,kdir,L,r2,I,isw)
% function Il=wblog_iter(k,kdir,L,r2,I,isw)

% outsourced from sweepPsi.m // Wb,Jan20,16

  if nargin && ischar(k)
     if isequal(k,'header'), s={'',''};
        if nargin>1 && isequal(kdir,'-b')
             s{2}=' using bond updates';
        else s{1}=' 2-site'; end
        wblog(1,'@%s%s DMRG%s',mfilename('class'),s{:});
        fprintf(1,'\n%s\n%s\n%s\n',...
        '────────────┬─────────────────────────────────┬───────────────────',...
        '  systime   │sweep  Nkept  exp(SEnt)  disc.wt │ energy/site       ',...
        '────────────┴─────────────────────────────────┴───────────────────');
     elseif ~isempty(regexp(k,'sep')), fprintf(1,'%s\n',...
        '──────────────────────────────────────────────────────────────────');
     else disp(k), wbdie('invalid usage'); end

     return
  end

  if isempty(isw)
     if kdir>0, sw='-->'; else sw='<--'; end
  else
     if L<100, f='%02g'; else f='%03g'; end
     if kdir>0
          sw=sprintf(['│%g.' f '⟩'],isw,k);
     else sw=sprintf(['⟨%g.' f '│'],isw,k); end
  end

  if isfield(I,'E0')
       E0=I.E0(end,1)/L; s={sw,sprintf('%12.9f',E0),''};
  else E0=[]; s={sw,repmat(' ',1,15),'\r'};
  end

  if isfield(I,'se')
       se=I.se; s{4}=sprintf('%10.6f',exp(se));
  else se=[]; s{4}=''; end

  if r2 % avoid '0.000e+00' string // Wb,Jun01,19
       s{5}=sprintf('%9.2e',r2);
  else s{5}=sprintf('%9g',r2); end

  if isfield(I,'NK'), NK=I.NK;
  elseif isfield(I,'Nkeep'), NK=I.Nkeep;
  else NK=nan; end

  if isfield(I,'wrn') && ~isempty(I.wrn), tag='';
     q=regexprep(I.wrn,'(WRN|ERR)(?@tag=$1;)\s*','');
     if ~isempty(tag)
        [i,w]=wblog('--hl-check',tag);
        if i, q=[w{1} q w{2}]; end
     end
     s{2}(end+1:16)=' '; s{2}=[s{2} q];
  end

  fprintf(['  %s  %s %5g %s %s │ %s\n',s{3}],...
     datestr(now,'HH:MM:SS'),s{1},NK,s{4:5},s{2});

  if nargout
     Il.E0=E0;
     if ~isempty(E0) && isfield(I,'e0') && isfield(I,'se')
        Il.e0=I.e0;
        Il.se=se;
        Il.rr=I.rr; 
        Il.r2=r2;
     end
     Il.NK=NK;
     Il.Dk=reshape([I.Dk(:); I.DX(:)],size(I.Dk,1),[]);

     if isfield(I,'Eg')
        Il.Eg=I.Eg(end);
        Il.ex=I.ex;
     end
  end

end

