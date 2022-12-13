function display(HAM,varargin)
% function display(H [,opts])
% Wb,Apr08,14

  if isequal(HAM,Hamilton1D)
     fprintf(1,'\n   Got empty Hamilton1D object.\n\n');
     return
  end

  if nargin>1
     getopt('INIT',varargin);
        if getopt('-v'), vflag=1;
        elseif getopt('-V'), vflag=2; else vflag=0; end
     getopt('check_error');
  else vflag=0;
  end

  e=consistency_check(HAM,'-q');

  s=regexprep(HAM.info.istr,'[\n\s]*$','');
  s=regexprep(s,char(10),[char(10) '   ']);
  s=regexprep(s,'\\n',[char(10) '   ']);

  fprintf(1,'\n   %s\n',s);
  L=numel(HAM.mpo);

  if e, e=100+e; else
    AK=load_dmrg_data(HAM,1,'AK','-t'); vac=2;
    if ~isempty(AK), d=getDimQS(AK);
       if any(d(:,1)>1), vac=vac-1;
          if vac==1, fprintf(1,'\n'); end
          printfc('  \e[31m NB! got non-vacuum at left boundary \e[0m');
          if d(1,1)>1
             E=QSpace(getIdentityQS(AK,1));
             fprintf(1,'\n'); display(E); fprintf(1,'\n');
          else
             fprintf(1,'[%s] (d=%g)\n',vec2str(AK.Q{1}(1,:),'-f'),d(end,1));
          end
       end
    else
       printfc('  \e[31m NB! got AK(1) not yet initialized !? \e[0m\n\n');
       e=1;
    end
  end

  if ~e
    AK=load_dmrg_data(HAM,L,'AK','-t');
    if ~isempty(AK), d=getDimQS(AK);
       if any(d(:,2)>1), vac=vac-1;
          if vac==1, fprintf(1,'\n'); end
          printfc('  \e[31m NB! got non-vacuum at right boundary \e[0m');
          if d(1,2)>1
             E=QSpace(getIdentityQS(AK,2));
             fprintf(1,'\n'); display(E); fprintf(1,'\n');
          else
             fprintf(1,'[%s] (d=%g)\n',vec2str(AK.Q{2}(1,:),'-f'),d(end,2));
          end
       end
    else
       printfc('  \e[31m NB! got AK(L) not yet initialized !? \e[0m\n\n');
       e=2;
    end
  end

  fprintf(1,'\n       %-27s %s\n','symmetry',get_sym_info(HAM));
  fprintf(1,  '       %-27s %g\n','system length L',L);

  if isfield(HAM.info,'mpo'), IH=HAM.info.mpo;

  if isfield(IH,'lmax') && ~isempty(IH.lmax), l=IH.lmax;
     if isfield(IH,'ldist') && ~isempty(IH.ldist)
        dx=IH.ldist(:,1);
        if numel(dx)<4, s=sprintf('[ %s ]',vec2str(dx,'-f'));
        else s=sprintf('[ %g..%g ]',min(dx),max(dx)); end
     else s=' !?'; end

     fprintf(1,'       %-27s ','coupling range in H');
     if l<=max(2,L/2)
        fprintf(1,'%s\n',s);
     else
        HH=HAM.info.HH; dx=abs(HH(:,1)-HH(:,3));
        i=find(dx>max(2,L/2));
        if ~isempty(i)
           if numel(i)==1
                s2=sprintf('got 1 long-ranged bond (L-%g)',L-dx(i));
           else s2=sprintf('got %g long-ranged bonds',numel(i));
           end
           s=sprintf('%-19s %s\n',s,[' \e[31m ' s2 '\e[0m']);
           printfc(s);
        end
     end
  end
  end

  [dM,dMs]=getDimMPO(HAM);

  if ~e
    Il=load_dmrg_info(HAM,'Il','-t');
    if ~isempty(Il) && isfield(Il,'Dk')
       Dk=cat(3,Il.Dk); dS=permute(Dk(:,1,:),[3 1 2]);
       if size(dS,2)>1
            dSs=sprintf('%d (%d)',max(dS(:,1)),max(dS(:,end)));
       else dSs=sprintf('%d',max(dS)); end
    else
       [dS,dSs]=getPsiDim(HAM,'-a','-q');
       dS(find(dS<0))=0;
    end
  else
     dS=getopt(HAM.info,'sweep','Nkeep',nan);
     dSs='undef';
  end

  ntype=max([size(HAM.oez,2), size(HAM.ops,2)]);
  dloc=zeros(ntype,2);
  for j=1:ntype
     s=getDimQS(HAM.oez(1,j).op); s=s(:,1); s(end+1:2)=s(1);
     dloc(j,:)=s;
  end

  s=sprintf('D=%s, d=%s',dMs{:});
  fprintf(1,'       %-27s %-22s','mpo dimension',s);

  full_mpo=usingFullMPO(HAM);
  if full_mpo
       printfc('\e[34;1m full MPO\e[0m\n');
  else fprintf(1,' pseudo MPO\n'); end

  s={ ['D=' dSs], sprintf('d=%s', vec2str(dloc(:,1),'-f')), '' };
  if dloc(1,1)~=dloc(1,2)
     s{2}=[ s{2}, sprintf(' (%s)', vec2str(dloc(:,2),'-f')) ];
  end
  if ntype==1, s{3}='uniform sites'; end
  fprintf(1,'       %-27s %-22s %s\n','MPS dimension',[s{1} ', ' s{2}], s{3});

  if ntype>1, stype=[HAM.mpo.stype]; end

  for j=1:ntype
     fmt='       %-27s [ %s ]\n'; s1=''; s2='';

     d=getDimQS(HAM.oez(1,j).op);
     if d(1,1)~=d(end,1), s2=sprintf('(=%g)',d(end,1)); end

     if ntype>1
        i=find(stype==j); m=numel(i);
        s3=['at sites ' ind2str(i)];

        s1=sprintf('%g.',j);
        fprintf(1,'\n  %4s %-27s %3g %-13s ... %s\n',...
        s1,'local state space dim (d)',d(1,1),s2,s3);
     end

     fprintf(1,fmt,'nops (irop dim)',ivec2str([HAM.ops(:,j).dop]));
     fprintf(1,fmt,'hconj',ivec2str([HAM.ops(:,j).hconj]));

     ff=[HAM.ops(:,j).fermionic]; fmt='       %-27s [ %s ]\n';
     if any(ff) 
          fprintf(1,fmt,'fermionic operators',ivec2str(ff));
     end
  end

  if isstruct(HAM.mpo), kk=[]; 
     for k=1:L
        if isfield(HAM.mpo(k).user,'ops'), kk(end+1)=k; end
     end
     if ~isempty(kk), fprintf(1,...
        '\n   NB! got %g user specified site(s) through HAM.mpo([%s]).user\n',...
        numel(kk),vec2str(kk,'sep',','));
     end
  end

  if ~isempty(HAM.store)
       s=sprintf('global workspace (''%s'')',HAM.store);
  else s=sprintf('files %s_###.mat',repHome(HAM.mat));
  end

  if ~e
     fprintf(1,'\n   Data is stored in %s\n',s);
     getCurrentSite(HAM,'-v');
  else
     if e<100
       printfc('\n  \e[31m [data store in %s''\n   not yet initialized?] \e[0m \n',s);
     else
       printfc('\n  \e[38;5;235m [%s]  \e[0m \n',s);
     end
  end
  fprintf(1,'\n');

  if vflag
     disp(HAM.info.param); if vflag>1
     subsref(HAM,struct('type','.','subs','ops')), end
  end

end

% -------------------------------------------------------------------- %
function s=ivec2str(q)
  if numel(q)<8
       s=vec2str(q,'-f');
  else
     [~,I,D]=uniquerows(double(diff(q')==0));
     q=matcell(q);
     for i=1:numel(q), q{i}=num2str(q{i}); end
     for i=find(D>10)', ix=I{i}+1;
         q{ix(1)}= sprintf('..(%d)..',numel(ix)-2);
         q(ix(2:end-1))=[];
     end
     s=strjoin(q);
  end
end

% -------------------------------------------------------------------- %

