function [R,dd]=LoadRData(sym,varargin)
% function [R,dd]=LoadRData(sym,qpat)
%    
%    load specified RData structure.
%    if qpat is not specified, *all* RSets are loaded.
%
% Options
%
%   '--fix'   check if matrix elements of generators
 %            can be fixed to exact rationals or sqrt(rationals)
%
% Wb,Jan30,15

% adapted from LoadCData // Wb,Jun07,15

  getopt('init',varargin);
     if getopt('-v');, vflag=1;
     elseif getopt('-v'), vflag=2; else vflag=0; end

     rfix=getopt('--fix');
     eps =getopt('eps',1E-12);

  varargin=getopt('get_remaining'); narg=length(varargin);

  D0=get_dir(sym,'RStore');

  if narg
     q=regexprep([varargin{:}],'.*(',''); q=regexprep(q,').*','');
     q=regexprep(q,',(\d+)\*',';$1','once');
     q=regexprep(q,'\*','');

     f=[ D0 '/(' q ').rep'];
     if exist(f,'file')
        I=load2(f,'-mat'); R=I.RSet;
     else error('Wb:ERR','\n   ERR file not found "%s"',f); end

     if vflag || rfix
        fprintf(1,'\n   %s (%s) d=%d\n',R.type,qmat2char(R.J),size(R.Z,1));
        if isfield(I,'err' ), fprintf(1,'   @ accuracy = %.3g',I.err); end
        if isfield(I,'istr'), fprintf(1,' // %s\n',I.istr); 
        else fprintf(1,'\n'); end

        f_=dir(f);
        fprintf(1,'   @ generated %s\n',datestr(f_.datenum));

        if rfix
           [R,P,Q,ee,ok]=rat_fix(R,eps);
           R=add2struct(R,P,Q,ee,ok);

           ee=cat(1,R.ee{:}); e=max(ee(:,end));

           if all(R.ok(:)), fprintf(1,[
             '\n   ok. fixed rationals in RData @ P<=%g, Q<=%g, ' ... 
             'e=%.3g\n\n'],max(P(:)),max(Q(:)),e);
           else fprintf(1,[
             '\n   WRN failed to fix rationals in %g generators ' ... 
             '(e=%.3g/%g)\n\n'],numel(find(ok(:)==0)),e,eps);
           end
        end
     end
  else
     ff=dir([D0 '/*.rep']); nf=numel(ff); jj=cell(1,nf); dd=zeros(nf,1);
     fprintf(1,'\r   Loading all %g R-data ... \r',nf);
     for i=nf:-1:1
         I=load([D0 '/' ff(i).name],'-mat'); I=I.RSet; R(i)=I;
         jj{i}=QSet2str(I.J);
         dd(i)=size(R(i).Z,1);
     end
     fprintf(1,'\r%70s\r','');

     jj=cat(1,jj{:}); [~,is]=sortrows([dd,double(jj)]);
     R=jj(is,:); dd=dd(is);
  end

end

% -------------------------------------------------------------------- %
function [R,P,Q,ee,ok]=rat_fix(R,eps)

   rsym=numel(R.Sz);
   Sz=cell(1,rsym); Sp=cell(1,rsym); ee=cell(rsym,2);
   ok=zeros(rsym,2); P=ok; Q=ok;

   for k=1:rsym
      [Sz{k},ee{k,1},ok(k,1),P(k,1),Q(k,1)]=rat_fix_1(cgs2double(R.Sz(k)),eps);
      [Sp{k},ee{k,2},ok(k,2),P(k,2),Q(k,2)]=rat_fix_1(cgs2double(R.Sp(k)),eps);
   end

end

% -------------------------------------------------------------------- %

function [Sp,e3,ok,Pmax,Qmax]=rat_fix_1(Sp,eps)

   i=find(Sp); q0=full(Sp(i));

   [q,Ir]=wbrat(q0,'-r','-q');

   e3=[norm(Ir.ee), norm(Ir.relerr), norm(q0-q)/norm(q)];
   if e3(end)<eps
        ok=1; Sp(i)=q;
   else ok=0; end

   Pmax=max(abs(Ir.P));
   Qmax=max(abs(Ir.Q));

end

% -------------------------------------------------------------------- %

