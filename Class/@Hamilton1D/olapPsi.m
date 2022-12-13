function X=olapPsi(HAM,varargin)
% function X=olapPsi(HAM [,opts])
%
%    '--x2'    exchange nearest-neighbor pairs (i.e. swap of rung for ladder)
%    '--RL'    compute overlap with left/right reflected wavefunction
%
% Wb,Dec15,17

% cd $LBIN/ts08.27003; load Wb161214o_tst_Hamilton1D_HBL51.mat % about 4min
% test for shorter ladder first!

% cd $LBIN/ts08.26650; load Wb161213v_tst_Hamilton1D_HBL51.mat; setRCStore rcl
% q=olapPsi(HAM,'-v');        % => x=1 @  6E-13
% q=olapPsi(HAM,'-v','--x2'); % => x=1 @ -2E-05
% q=olapPsi(HAM,'-v','--RL'); % => x=1 @ -6E-07

  nl=' ... \r';
  getopt('INIT',varargin);
     xflag=getopt('--x2');
     rflag=getopt('--RL');

     if getopt('-v'), vflag=1;
     elseif getopt('-V'); vflag=2; nl='\n'; fprintf(1,nl);
     else vflag=0; end

  getopt('check_error');

  L=length(HAM.mpo);

  if xflag && mod(L,2), wberr('invalid usage (L=%g !?)',L); end

  if rflag
     if xflag
        fmt=['\r   calculate overlap k=%2g/%g:  [ %s - %s X %s - %s ]' nl];
        for k=1:2:L/2
           A(1)=load_dmrg_data(HAM,  k,  'AK');
           A(2)=load_dmrg_data(HAM,  k+1,'AK');
           B(1)=load_dmrg_data(HAM,L-k,  'AK');
           B(2)=load_dmrg_data(HAM,L-k+1,'AK');

           if vflag, fprintf(1,fmt,k,L, ...
              A(1).info.itags{3}, A(2).info.itags{3},...
              B(1).info.itags{3}, B(2).info.itags{3});
           end

           for l=1:2
              t=B(l).info.itags; t{3}=A(l).info.itags{3};
              for i=1:2, t{i}=['R' t{i}(2:end)]; end
              B(l).info.itags=t;
           end

           if k>=L/2, Xlast=X; end

           if k>1, Q={X,A(1)}; else
              Q=A(1);
              B(2).info.itags{2}=A(1).info.itags{1};
           end
           X=contractQS(B(1),'*',{B(2),'*',{Q,A(2)}});

           if numel(X.Q)>3, wblog('WRN',''); keyboard, end
        end

        if ~mod(L,2), Xlast=X; end

        [Xlast.info.itags,ic]=adapt_itagsLR(Xlast.info.itags);

        X=contractQS(Xlast,ic,X);
     else
        fmt=['\r   RL calculate overlap k=%2g/%g:  [ %s -- %s ]' nl];

        for k=1:(L+1)/2
           A1=load_dmrg_data(HAM,  k,  'AK');
           A2=load_dmrg_data(HAM,L-k+1,'AK');

           if numel(A1.Q)==4
              A1.info.itags{end}=['X' A1.info.itags{end}];
           end

           if vflag, fprintf(1,...
              fmt, k,L, A1.info.itags{3}, A2.info.itags{3});
           end

           t=A2.info.itags;
           t{3}=A1.info.itags{3}; for i=1:2, t{i}=['R' t{i}(2:end)]; end
           A2.info.itags=t;

           if k>=L/2, Xlast=X; end

           if k>1, Q={X,A1}; else
              Q=A1;
              A2.info.itags{2}=A1.info.itags{1};
           end
           X=contractQS(A2,'!1*',Q);
           if numel(X.Q)>3, wblog('WRN',''); keyboard, end
        end

        fprintf(1,'\n');

        if ~mod(L,2), Xlast=X; end

        [Xlast.info.itags,ic]=adapt_itagsLR(Xlast.info.itags);

        X=contractQS(Xlast,ic,X);
     end
  else
     if xflag
        fmt=['\r   X2 calculate overlap k=%2g/%g:  [ %s X %s ]' nl];
        for k=1:2:L
           A1=load_dmrg_data(HAM,k,  'AK'); t1=A1.info.itags{3};
           A2=load_dmrg_data(HAM,k+1,'AK'); t2=A2.info.itags{3};

           if vflag, fprintf(1,fmt,k,L,t1,t2); end

           if k>1, Q={X,A1}; else Q=A1; end
           X=contractQS(Q,A2);
           i=findstrc(X.info.itags,'^s\d');
           if numel(i)==2
                X.info.itags(i)=X.info.itags(fliplr(i));
           else X.info.itags, wberr('unexpected itags !?');
           end

           if numel(A1.Q)==4
              A1.info.itags{end}=['X' A1.info.itags{end}];
           end
           X=contractQS(A2,'!2*',{A1,'*',X});
        end
     else
        fmt=['\r   calculate plain overlap k=%2g/%g:  [ %s ]' nl];
        for k=1:L
           Ak=load_dmrg_data(HAM,k,'AK');
           if vflag, fprintf(1,fmt,k,L,Ak.info.itags{3}); end
           if k>1, Q={X,Ak}; else Q=Ak; end
           if numel(Ak.Q)==4
              Ak.info.itags{end}=['X' Ak.info.itags{end}];
           end
           X=contractQS(Ak,'!2*',Q);
        end
     end
  end

  if vflag>1, fprintf(1,nl);
  elseif vflag, fprintf(1,'\r%80s\r',''); end

end

% -------------------------------------------------------------------- %

function [t,ic]=adapt_itagsLR(t)

  i=findstrc(t,'Psi');
  if ~isempty(i)
       j=1:3; j(i(1))=[]; ic=sprintf('!%g*',i(1));
  else j=[1 2]; ic='*'; end

  q=[ t{j(1)}(1), t{j(2)}(1) ];
  t{j(1)}(1)=q(2);
  t{j(2)}(1)=q(1);

end

% -------------------------------------------------------------------- %

