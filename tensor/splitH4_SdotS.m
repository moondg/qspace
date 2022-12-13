function [X1,X2,Iout]=splitH4_SdotS(H4,varargin)
% function [X1,X2,I]=splitH4_SdotS(H4 [,opts])
%
%    Split input H4 with index order (i,j;i',j'), i.e. 2-in, 2-out (')
%    into an X1.X2' structure, i.e. X1 \cdot X2^\dagger.
%
% Options
%
%   'stol',.. cutoff (tolerance) on singular values
%   'istr',.. just info string to be used when -v is set
%   '--notr'  make traceless (i.e., subtract trace)
%   '-v'      verbose mode (e.g. report result on consistency check)
%
% Wb,Apr05,16

% adapted from Hamilton1D/private/setup_HeisenbergTriLadder.m

  getopt('init',varargin);
     vflag=getopt('-v');
     stol =getopt('stol',1E-12);
     istr =getopt('istr','H4' );
     notr =getopt('--notr');
  getopt('check_error');

  qdir=getqdir(H4(1),'-s');
  if ~isequal(qdir,'++--'), H4(1), qdir
     error('Wb:ERR','\n   ERR invalid input H4');
  end

  H41=sum(QSpace(H4));

  E1=QSpace(getIdentityQS(H41,[1 3])); Z1=QSpace(getIdentityQS(E1,1,'-0'));
  E2=QSpace(getIdentityQS(H41,[2 4])); Z2=QSpace(getIdentityQS(E2,1,'-0'));

  A1=QSpace(getIdentityQS(E1,Z1,2));
  P1=QSpace(permuteQS(contractQS(A1,2,Z1,'2*'),'132'));

  same12=isequal(E1,E2);
  if same12
     A2=A1; A4=A1; P2=P1;
  else
     A2=QSpace(getIdentityQS(E2,Z2,2));
     P2=QSpace(permuteQS(contractQS(A2,2,Z2,'2*'),'132'));

     A4=QSpace(getIdentityQS(E1,E2));
  end

  nH=numel(H4); tr4=[];
  Iout.stol=stol; Iout.Xb=QSpace; Iout.SS={};

  if notr
   % make traceless, since this only adds a trivial offset
   % NB! for Hamiltions with a single local multiplet, this just
   % removes the irop ontribution in the scalar sector
   % (since is the one that acts like an identity operator!)
   % => mostly leave as is, so to have consisteny in ground state
   % energies for benchmarks as is; it should not matter much
   % for the numerical efficiency! // Wb,Apr03,20
     U2=getIdentityQS(E1,E2); D=getDimQS(U2); D=D(end);
     E4=QSpace(contractQS(U2,3,U2,'3*'));
     for i=1:nH
         H2=QSpace(contractQS(U2,'12*',contractQS(H4(i),'34',U2,'12'),'12'));
         tr4(i)=trace(H2)/D;
         H4(i)=H4(i) - tr4(i)*E4;
         if 1 || vflag, wblog('NB!','subtracting trace %.8g',tr4(i)); end
     end
  end
  Iout.tr4=tr4;

  X1=QSpace(1,nH); X2=QSpace(1,nH);
  H2=QSpace(1,nH); H_=QSpace(nH,2);

  for i=1:nH
     X2(i)=contractQS(contractQS(P1,'12*',H4(i),'13'),'23',P2,'12*');
     X2(i)=skipzeros(X2(i));
     Iout.Xb(i)=X2(i); % `bond' term

     X1(i)=getIdentityQS(X2(i),1);
     if ~isequal(X1(i).Q{1},X2(i).Q{1})
        error('Wb:ERR','\n   ERR unexpected X2'); end

     nd=numel(X2(i).data); ss=cell(1,nd);
     for l=1:nd
        [u,s,v]=svd(X2(i).data{l});
           s=diag(s); ss{l}=s; j=find(s>stol);
           u=u(:,j); rs=diag(sqrt(s(j))); v=v(:,j);

        X1(i).data{l}=u*rs;
        X2(i).data{l}=rs*v';
     end

     Iout.SS{i}=sort(cat(1,ss{:}),'descend');

     X1(i)=contractQS(P1,3,X1(i),1);
     X2(i)=permuteQS(contractQS(P2,3,X2(i),2),'213*');

     X1(i).info.otype='operator';
     X2(i).info.otype='operator';

     H2(i)=contractQS(A4,'12*',contractQS(A4,'12',H4(i),'34'),'23');

     H_(i,1)=contractQS(A4,'12*',...
       contractQS(X1(i),'23',contractQS(A4,2,X2(i),'1*'),'14'),'13');
     H_(i,2)=contractQS(A4,'12*',...
       contractQS(X1(i),'13*',contractQS(A4,2,X2(i),'2'),'14'),'13');
     ee(i,1:3)=[normQS(H_(i,1)-H2(i)), normQS(H_(i,2)'-H2(i)), normQS(H2(i))];
  end

  if nH==1, Iout.SS=Iout.SS{1}; end

  e=sum(norm(ee(:,1:2)./ee(:,[3 3])),2); i=find(e>1E-12,1);
  if ~isempty(i), error('Wb:ERR',...
    '\n   ERR got %s inconsistency (@ %.3g; %.3g) !?',istr,e(i),ee(i,3));
  elseif vflag
     wblog(' * ','got consistent %s [@ %.2g]',istr,e);
  end

end

