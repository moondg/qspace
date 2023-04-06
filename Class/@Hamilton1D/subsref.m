function varargout=subsref(HAM,S)
% function x=subsref(HAM,S)
%
%    Customizing index or field access into HAM.
%
%    HAM(k)  load data set for iteration k
%            where k=0 loads *_info.mat file, instedd
%
%    HAM.ops by itself and if no output is requested
%            displays fields in HAM.ops; of some letters
%            in .ops are capitalized to increases verbosity
%            e.g. also showing ops().op QSpace, etc.
%
% Wb,Apr08,14 ; Wb,Sep16,20

  if isequal(S(1).type,'()')
   % further arguments specify individual variables to load
   % e.g. HAM(k,'info') // Wb,Oct04,16
     k=S(1).subs{1};
     if ischar(k)
        switch k
           case {'end','last'}, k=numel(HAM.mpo);
           case {'kc'},         k=getCurrentSite(HAM);
           otherwise wberr('invalid usage');
        end
     end
     if ~isnumeric(k) || numel(k)~=1, k
        wberr('invalid data reference');
     end

     if numel(S(1).subs)>1 && isequal(S(1).subs{2},'-Tr')
          x=load_trotter_data(HAM,k,S(1).subs{3:end});
     else x=load_dmrg_data   (HAM,k,S(1).subs{2:end});
     end

     if numel(S)>1
        x=builtin('subsref',x,S(2:end)); end

     varargout={x};

  elseif numel(S)==1 && ...
     isequal(S.type,'.') && ~isempty(regexpi(S.subs,'^(ops|oez)$'))

     vflag=sum(S.subs<'a');
     S.subs=lower(S.subs);

     x=builtin('subsref',HAM,S);
     if nargout, varargout={x}; return; end

     dd=[]; for i=1:numel(x), dd(i)=length(x(i).info); end
     fmt=sprintf('%%-%ds',max(26,max(dd)));

     sx=size(x); mops=sx(2)>1; ds={};
     sh={'info'};
     if mops, fmt=[fmt ' %8s']; sh{2}='dims'; end

     if sx(2)>1, fprintf(1,'\n   got %g different site types\n',sx(2)); end
     for j=1:sx(2)
        if sx(2)>1, fprintf(1,'   site type %d\n',j); end
        fprintf(1,['\n%5s ' fmt ' dop  q-irop  flags\n'],'',sh{:});

        for i=1:sx(1), q=x(i,j);
           if mops
              E=q.op; if numel(E.Q)==3, E.info.otype='operator'; end
              dd=getDimQS(getIdentity(E));
              ds=sprintf('x%d',dd(2,:)); ds={ds(2:end)};
           end

           s={ num2str(q.dop), ['(' vec2str(q.qop,'-f') ')'] };
           if q.fermionic, s{end+1}='fermionic';    end
           if q.hconj,     s{end+1}='include H.c.'; end
           fprintf(1,['%4g. ' fmt ' %2s %8s  %s\n'],i,q.info,ds{:},s{1:2},...
           strjoin(s(3:end),', '));
        end
        fprintf(1,'\n');

        if vflag==1
           for i=1:sx(1), info(x(i,j).op); end
        elseif vflag
           q=[x(:,j).op];
           display(q(:),'nm',sprintf('ops(%%g,%g)',j),'-v');
        end
     end

  elseif numel(S)==2 && ...
     isequal(S(1).type,'.') && ~isempty(regexpi(S(1).subs,'^info$')) && ...
     isequal(S(2).type,'.') && ~isempty(regexpi(S(2).subs,'^HS$'))

     HS=HAM.info.HS; n=numel(HS);
     for i=1:n
        q=HS(i).iop; s={'',''};
        if numel(q)==4 && q(1,2)==q(2,2)
             s=sprintf('%4d %4d  %2d: %s',q(1:3),HAM.ops(q(3)).info);
        else s=sprintf(' %4d',q); s(1)=[];
        end
        fprintf(1,'%4d. %8s %s\n',i,num2str2(HS(i).hcpl),s);
     end

  else
     [varargout{1:nargout}] = builtin('subsref',HAM,S);
  end

end

