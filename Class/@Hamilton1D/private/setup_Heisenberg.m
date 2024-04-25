function [HAM]=setup_Heisenberg(varargin)
% function [HAM]=setup_Heisenberg(varargin)
% Wb,Jan15,15

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     qloc = getopt('qloc',[]);
     qend = getopt('qend',[]);
     sym  = getopt('sym','SU2'); % value 'Spin' uses abelian
     L    = getopt('L',[]);

     B    = getopt('B',0.);

     J    = getopt('J', []);
     Jz   = getopt('Jz',[]);
     J2   = getopt('J2',[]);
     J3   = getopt('J3',[]);
     perBC= getopt('-perBC');
     tflag= getopt('-t');

     Ising= getopt('--Ising');

     if getopt('SfacQuella')
          Sfac=sqrt(2);
     else Sfac=getopt('Sfac',[]); end

     waklt=getopt('waklt','');

     use_mpo = getopt('--mpo',{'not_specified'});

  if isempty(L), L=getopt('get_last',[]); else
  getopt('check_error'); end

  if ~iscell(use_mpo)
     wbdie('invalid usage (--mpo requires cell)');
  end

% -------------------------------------------------------------------- %
% check L vs J specifications

  if length(L)>1
     if ~isempty(J), wbdie(['invalid usage ' ... 
        '(J specified twice, once implicitely through L)']); end
     J=adapt_Jformat(L); L=[];
  elseif isempty(L)
     J=adapt_Jformat(J);
  elseif isempty(J), J=1;
  end

  if isempty(L)
     if     size(J,1)>1, L=size(J,1)+1; if perBC, L=L-1; end
     elseif size(B,1)>1, L=size(B,1);
     else wbdie('missing length specification'); end
     wblog(' * ','got implicit length specification (L=%g)',L);
  end
  if L<=1 || L~=round(L)
     wbdie('invalid usage (L=%g) !?',length(J),L);
  end

  gotB=norm(B); gotDz=0; gotBz=0; gotBx=0;
  if ~gotB && Ising, gotB=-1; end
  if gotB
     if isequal(sym,'SU2'), sym='Spin'; end
     if size(B,1)==1 && size(B,2)>3, B=B.'; end
     if numel(B)>3 && size(B,1)~=L
        wbdie('invalid B (len=%g/%g) !?',numel(B),L);
     end
     gotBz=norm(B(:,1)); m=size(B,2);
     if m>1
        if m>3, wbdie('invalid usage'); end
        if m>1, gotDz=norm(B(:,2)); end
        if m>2, gotBx=norm(B(:,3)); end
     end
  end

  gotJz=0;
  if ~isempty(Jz) && ~isequal(J,Jz), gotJz=1; end

  if gotJz || gotBz || gotBx || gotDz, Aflag=1;
      if ~isequal(sym,'Spin')
         wblog('WRN','got sym=''%s'' -> ''%s'' (abelian)',sym,'Spin');
         sym='Spin';
      end
  else Aflag=isequal(sym,'Spin');
  end

  if perBC, L_=L; else L_=L-1; end
  if size(J,1)>1
     q=uniquerows(J); if size(q,1)==1, J=q; end
  end

  if isempty(J), J=1;
  elseif numel(J)~=1, s=size(J);
     if s(1)==1
        if s(2)>=L-1, J=J'; end
     elseif s(1)<s(2), wbdie('unexpected J !?');
     end
  end; s=size(J);

  if s(1)~=1 && s(1)~=L_
     wbdie('invalid usage (system length mismatch %g/%g)',s(1),L);
  end

  ncpl=0;
  if ~Aflag || s(2)>2
     while norm(J(:,end))==0, J(:,end)=[]; end
     ncpl=size(J,2)-1;
     if ncpl>8, wbdie('unexpected ncpl=%g !?',ncpl); end
     if ncpl && gotJz
        wbdie('invalid usage (got Jz with size(J,2)>1)');
     end
  end

  gotJ2=~isempty(J2);
  if gotJ2 && size(J2,1)~=1
     if size(J2,1)~=L_, wbdie('invalid size of J2'); end
     q=uniquerows(J2); if size(q,1)==1, J2=q; end
  end

  gotJ3=~isempty(J3);
  if gotJ3 && size(J3,1)~=1
     if size(J3,1)~=L_, wbdie('invalid size of J3'); end
     q=uniquerows(J3); if size(q,1)==1, J3=q; end
  end

  if Ising
     if gotDz, wbdie('invalid usage (got Dz with Ising)'); end
     if gotBz && gotBx, wbdie('invalid usage (got Bz and Bx with Ising)'); end
     if gotJ2 && gotJ3, wbdie('invalid usage (got J2 or J3 with Ising)'); end
     if ncpl, wbdie('invalid usage (got ncpl ith Ising)'); end
     if ~isempty(qend), wbdie('invalid usage (got ncpl ith Ising)'); end
  end

% -------------------------------------------------------------------- %

  if ~isempty(regexp(sym,'^SU\d+')) && str2num(sym(3:end))>4, o={'-V'};
  else o={'-v'}; end

  param=add2struct(L,J,J2,J3,B,perBC,Ising,gotJz,Aflag,qloc,ncpl,waklt);

  HAM=struct(Hamilton1D);

  if Aflag, q=unique(J(:,1));
     if Ising, sJ='Jz'; else sJ='J'; end
     if numel(q)>1
          jstr=sprintf('%s=[%.4g .. %.4g]',sJ,q([1,end]));
     else jstr=sprintf('%s=%g',sJ,q); end

     [s1,s2]=size(J);
     if s2>1
        if s2==2, q=unique(J(:,2)); gotJz=norm(diff(J,[],2));
           if Ising, sJ='Jx'; else sJ='Jz'; end
           if numel(q)>1
                jstr=[jstr, sprintf(', %s=[%.4g .. %.4g]',sJ,q([1,end]))];
           else jstr=[jstr, sprintf(', %s=%g',sJ,q)]; end
        elseif s1==1
           s=sprintf(',%.3g',J);
           jstr=sprintf('J=[%s]',s(2:end));
        else wbdie('invalid usage'); end
     end
     HAM.info.mpo.J=J;

  elseif size(J,1)==1
     if ncpl<=0 || norm(J(2:end))==0
        jstr=sprintf('J=%g',unique(J(:,1)));
     elseif numel(J)==2
        jstr=sprintf('BLQ J=[%g %g]',J);
     else
        s=sprintf(',%.3g',J);
        jstr=sprintf('J=[%s]',s(2:end));
     end
     HAM.info.mpo.J=J;
  else
     n=size(J,2); jstr=cell(1,n);
     for i=1:n
        q=[min(J(:,i)),max(J(:,i))];
        if diff(q)>1E-8
             jstr{i}=sprintf('%.3g~%.3g',q);
        else jstr{i}=sprintf('%.4g',q(1)); end
     end
     if n<=1 || norm(q(:,2:end))==0, jstr=['J=' jstr{1}];
     elseif n==2
        jstr=sprintf('BLQ J=[%s, %s]',jstr{1:2});
     else
        jstr(2,1:end-1)={', '};
        jstr=['J=[' jstr{:} ']'];
     end
     HAM.info.mpo.J=J;
  end

  HAM.info.mpo.J2=J2;
  if ~isempty(J2)
     if size(J2,2)>1
        wbdie('invalid J2 (biquadratic term not yet implemented)'); end
     jstr=regexprep(jstr,'\<J\>','J_1');
     if numel(J2)>1
          jstr=[jstr, sprintf(', J_2=[%.4g .. %.4g]',J2(1),J2(end))];
     else jstr=[jstr, sprintf(', J_2=%g',J2)];
     end
  end

  HAM.info.mpo.J3=J3;
  if ~isempty(J3)
     if size(J3,2)>1
        wbdie('invalid J3 (biquadratic term not yet implemented)'); end
     jstr=regexprep(jstr,'\<J\>','J_1');
     if numel(J3)>1
          jstr=[jstr, sprintf(', J_3=[%.4g .. %.4g]',J3(1),J3(end))];
     else jstr=[jstr, sprintf(', J_3=%g',J3)];
     end
  end

  Sflag=(~gotB && ~gotDz && ~gotJz);

  if isequal(sym,'SU2')
     if isempty(qloc), qloc=[1];
     elseif numel(qloc)~=1 || qloc<1, wbdie('invalid spin S=%g',S); end
     [S,IS]=getLocalSpace('Spin',qloc/2,o{:});
     if S.data{1}<0, S.data{1}=-S.data{1}; end

     if ~isfield(IS,'istr')
        IS.istr=sprintf('Sloc=%s',num2rat(qloc/2));
     end

     if Ising, s='Ising'; else s='Heisenberg'; end
     HAM.info.istr=sprintf('%s (%s, %s) having %s',s,jstr,sym,IS.istr);

  elseif isequal(sym,'Spin')
     if ncpl && (gotB || gotBx || gotDz || ~isempty(waklt))
        wbdie('invalid usage (got abelian U1 mode; no ncpl supported here)');
     end

     if isempty(qloc), qloc=[1];
     elseif numel(qloc)~=1 || qloc<1 || mod(qloc,1)
        wbdie('invalid spin S=%g',qloc/2);
     end

     if gotBx
          [S,IS]=getLocalSpace('Spin',qloc/2,'--nosym',o{:});
     else [S,IS]=getLocalSpace('Spin',qloc/2,'-A',o{:});
     end

     if gotB || isempty(S(1).info.qtype)
     if ncpl, S=sum(S);
     else 
        for i=1:numel(S)
           if isequal(S(i).Q{1},S(i).Q{2})
              S(i)=fixScalarOp(S(i));
              iz=i; break
           end
        end
        if gotDz
           iz2=numel(S)+1;
           q=S(iz)*S(iz); d=getDimQS(IS.E);
           q=q-(trace(q)/d(end,1))*IS.E;
           if norm(q)<1E-8
              wblog('WRN','ignoring irrelevant Sz^2 operator (Id)');
              B(:,2)=0;
           end
           S(iz2)=q;
        end
        if gotBx || gotJz
           if iz==1, ix=2;
           else wbdie('unexpected iz=%g',iz); end
        end
     end
     end

     if ~isfield(IS,'istr')
        IS.istr=sprintf('Sloc=%s',num2rat(qloc/2));
     end

     qq={gotBz,1,'B'; gotDz,2,'Dz'; gotBx,3,'Bx' }; bstr='';
     for i=1:size(qq,1)
        if qq{i,1}, q=unique(B(:,qq{i,2}));
           if norm(q)
              if numel(q)>1
                   bstr=sprintf(', %s=%g..%g',qq{i,3},q([1,end]));
              else bstr=sprintf(', %s=%g',qq{i,3},q); end
           end
        end
     end

     if Ising, s='Ising'; else s='Heisenberg'; end
     s=sprintf('%s (U1 %s, %s%s)',s,sym,jstr,bstr);
     HAM.info.istr=[s ' having ' IS.istr];

  elseif ~isempty(regexp(sym,'^S[Op]\d'))
     if gotB, wbdie('magnetic field not yet implemented (B=%g)',gotB); end

     [S,IS]=getLocalSpace(sym,qloc,o{:});
     HAM.info.istr=sprintf('Heisenberg %s, %s',jstr,IS.istr);
  else
     if gotB, wbdie('magnetic field not yet implemented (B=%g)',gotB); end

     [S,IS]=getLocalSpace('Spin',sym,qloc,o{:});
     S=normalize_spin_operator(S,Sfac);

     HAM.info.istr=sprintf(['Heisenberg (%s)\n' ... 
    'having %s with qloc=(%s)'],jstr,IS.istr,sprintf('%g',IS.qloc));
  end

  wblog('==>','%s',HAM.info.istr);

  HAM.info.param=param;
  HAM.info.IS=IS;

  HAM.oez=init_ops(IS.E,'local identity operator (E)');

  if numel(S)==1
     HAM.ops=init_ops(S,'spin operator S','~hconj');
  else
     if numel(S)~=3, wbdie('invalid S-op'); end
     if Ising
        HAM.ops=[
           init_ops(S(1),'spin-op (Ising Sz)','~hconj');
        ];
        if gotB>0 || gotJz
           m=size(B,3); if m>3 || m>1 && norm(B(:,1)), whos B
             wbdie('invalid usage (Ising @ L+%g)',L); end
          if gotBx || gotJz
             HAM.ops(end+1,1)=init_ops(...
             (1/sqrt(2))*(S(2)+S(2)'),'spin-op (Ising Sx)','~hconj');
          end
        end
     elseif Sflag
        HAM.ops=[
           init_ops(sqrt(1/2)*makeIrop(S(1))+S(2),'spin-op (sym)','-hconj');
        ];
     else
        HAM.ops=[
           init_ops(S(1),'spin-op Sz','~hconj');
           init_ops(S(2),'spin-op S_/sqrt(2)','-hconj');
        ];
     end
  end

  if isset('iz2')
     HAM.ops(iz2).info='Sz^2 (loc. anisotropy)';
  end

  if isequal(qloc,qend), qend=[]; end
  if ~isempty(qend) , qend_=qend;
     if gotB, wbdie('magnetic field not yet implemented (B=%g)',gotB); end
     if perBC, wbdie('invalid usage (perBC xor qend!)'); end

     m=size(qend,1);
     for i=1:m
        if isequal(sym,'SU2')
           [Se,Ie]=getLocalSpace('Spin',qend(i,:)/2,o{:});
           if Se.data{1}<0, Se.data{1}=-Se.data{1}; end
        else
           [Se,Ie]=getLocalSpace('Spin',sym,qend(i,:),o{:});
           Se=normalize_spin_operator(Se,Sfac);
        end

        if m==1, s=''; else if i==1, s='_L'; else s='_R'; end, end
        HAM.oez(1,1+i)=init_ops(Ie.E,['end-spin (E' s ')']);
        HAM.ops(1,1+i)=init_ops(Se,  ['end-spin (S' s ')'],'~hconj');
     end

     qend=m;
  else qend=0; end

  stype=[];

  if ~isempty(waklt)
     if numel(S.Q)~=3, wbdie('irop expected for S in case of waklt'); end
     if isequal(waklt,'Quella')
        q=getQuellaAKLT(waklt,qloc); n=numel(q);
        for i=1:n
           s=sprintf('%s [%s]',waklt,iff(mod(i,2),'a','b'));
           Sab(i)=init_ops(q(i),s,'~hconj');
        end
        if ncpl>0
           if n~=2*ncpl, wbdie(...
              'invalid usage (n_aklt = %g <> %g !?)',n,2*ncpl); end
           HAM.ops(2:n+1,1)=Sab;
        else
           HAM.ops(1:n,1)=Sab;
        end
     else wbdie('invalid waklt');
     end
  elseif ncpl
     if numel(S)~=1 || numel(S.Q)~=3
        wbdie('irop expected for S in case of (S.S)^k'); end
     if size(HAM.ops,1)~=1, wbdie('check setup'); end

     if size(J,1)==1
        [HAM.ops(1:2,1),Ix]=get_poly_ops(S,J); % '--notr'
        HAM.info.param.trOps=Ix.tr4;
        if qend
           [HAM.ops(1:2,2),Iy]=get_poly_ops(Se,J,'Send'); % ,'--notr'
           stype=ones(1,L); stype([1 end])=2;
        end
        J=1; ncpl=-1; % ncpl below indicates to include Sa.Sb' (!)
     else
        HAM.ops(2:2+2*ncpl-1,1)=get_poly_ops(S,ncpl+1);
        if qend
           HAM.ops(2:2+2*ncpl-1,2)=get_poly_ops(Se,ncpl+1,'Send');
           stype=ones(1,L); stype([1 end])=2;
        end
     end
  end

  HAM.info.lops=S;

  if ncpl<0 && size(HAM.ops,1)~=2, wbdie(...
    'got %g ops with ncpl=%g !?',size(HAM.ops,1),ncpl); end
  acpl=abs(ncpl);

  n=L-1+perBC;
  nops=size(HAM.ops,1)-acpl;

  if size(J, 1)==1, J =repmat(J, L_,1); end
  if size(J2,1)==1, J2=repmat(J2,L_,1); end
  if size(J3,1)==1, J3=repmat(J3,L_,1); end

  HH=zeros(n*max(1,1+ncpl),5); l=1;

  if perBC
     wblog('NB!','using periodic BC (interleaved setup)'); 
     L2=ceil(L/2);

     XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
     HAM.info.XY=XY;

     if gotJ2, wbdie('J2 not yet implemented for perBC'); end
     if gotJ3, wbdie('J3 not yet implemented for perBC'); end
     if qend, wbdie( ...
       'end-spin (open right perBC) not yet implemented'); end

     for i=1:n
        if mod(i,2)
             k=(i-1)/2+1;
        else k=n+1-i/2;
        end

        if ncpl>=0
           HH(l,:)=[ [i-1, 1  ], [i+1, 1    ], J(k,1  )]; l=l+1; for j=1:ncpl
           HH(l,:)=[ [i-1, 2*j], [i+1, 2*j+1], J(k,1+j)]; l=l+1; end
        else
           HH(l,:)=[ [i-1, 1  ], [i+1, 2    ], J(k,1  )]; l=l+1;
        end
     end

     HH(1:1+acpl,1)=1;
     HH(end-acpl:end,3)=L;

  else

     if gotJ3
        XY=[ (1:L)', zeros(L,1) ]; XY(2:3:end,2)=1; XY(3:3:end,2)=2;
        HAM.info.XY=XY;
     elseif gotJ2
        XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
        HAM.info.XY=XY;
     end

     for i=1:n
        if ncpl>=0
           HH(l,:)=[ [i, 1  ], [i+1, 1    ], J(i,1  )]; l=l+1; for j=1:ncpl
           HH(l,:)=[ [i, 2*j], [i+1, 2*j+1], J(i,1+j)]; l=l+1; end
        else
           HH(l,:)=[ [i, 1  ], [i+1, 2    ], J(i,1  )]; l=l+1;
        end

        if gotJ2 && i+1<L
        HH(l,:)=[ [i, 1  ], [i+2, 1    ],J2(i,1  )]; l=l+1; end

        if gotJ3 && i+2<L
        HH(l,:)=[ [i, 1  ], [i+3, 1    ],J3(i,1  )]; l=l+1; end
     end

     if ncpl<0, HH(1:l-1,4)=2; end

     if qend
        wblog('NB!','using end-spins (%s)',sprintf('%g',qend_')); 
        stype=ones(1,L); stype([1,end])=2;
        if qend>1, stype(end)=3; end

        i=[ find((HH(:,1)==1 | HH(:,1)==L) & HH(:,2)>1)
            find((HH(:,3)==1 | HH(:,3)==L) & HH(:,4)>1)
        ];
        HH(i,:)=[]; l=size(HH,1)+1;
     end
  end

  if isequal(sym,'Spin')
     if Sflag
        if size(HAM.ops,1)~=1, wbdie('invalid usage'); end
     elseif ~Ising || gotJz
        if size(HAM.ops,1)~=2
           wbdie('invalid usage');
        end
        i=find(sum(HH(:,[2 4])-1,2)==0);

        H3=repmat(HH(i,:),1,1,2); HH(i,:)=[];
        H3(:,[2 4],2)=2;

        if size(J,2)==2
           if Ising
              H3(:,end,ix)=J(:,2);
           else 
              if isempty(Jz), Jz=J(:,2); end
              H3(:,end,iz)=Jz;
           end
        end

        HH=[HH; reshape(permute(H3,[3 1 2]),[],size(H3,2)) ];
     end
  end

  if gotB, l=size(HH,1)+1;
     if gotBx
        if Ising
           if HAM.ops(ix).hconj, wbdie('invalid Sx'); end
        else
           B(:,3)=B(:,3)/sqrt(2);
           if ~HAM.ops(ix).hconj, wbdie('invalid Sx !?'); end
        end
     end

     if size(B,1)==1, B=repmat(B,L,1); end
     for i=1:L
        if gotBz && B(i,1), HH(l,:)=[ [i,iz ], [i,iz ], B(i,1)]; l=l+1; end
        if gotDz && B(i,2), HH(l,:)=[ [i,iz2], [i,iz2], B(i,2)]; l=l+1; end
        if gotBx && B(i,3), HH(l,:)=[ [i,ix ], [i,ix ], B(i,3)]; l=l+1; end
     end
  end

% if 0
%  % duplicate H (for testing purposes)
%    banner('creating duplicate of H!');
%    H2=HH; n=max([HH(:,2); HH(:,4)]); H2(:,[2 4]) = H2(:,[2 4]) + n;
%    HH=[HH;H2]; HH(:,end)=HH(:,end)/2;
%    HAM.ops=[HAM.ops; HAM.ops];
% end

  if isequal(use_mpo,{'not_specified'})
     HAM=setup_mpo(HAM,HH,stype);
  else
     HAM.info.HH=HH;
     if isequal(use_mpo,{1}), use_mpo={}; end
     HAM=setup_mpo_full(HAM,use_mpo{:});
  end

  HAM.store=['DMRG_HB_' sym];

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

function J=adapt_Jformat(J)
   s=size(J);
   if s(1)==1
      if s(2)>3, J=J.'; end
   elseif s(1)<s(2)
      wbdie('got unexpected J of size %gx%g !?',s);
   end
end

% -------------------------------------------------------------------- %
% Wb,Sep10,15

function S=normalize_spin_operator(S,Sfac)

   S=skipzeros(QSpace(S));

   switch S.info.qtype
     case 'SU3'

     case 'SU4'
      % S.S = (15/8) * Id in the defining irrep // Wb,Feb17,16
      % consistent with (N^2-1)/(2*N) // see PHYS notes.
      % NB! nearest neighbor S.S in the defining representation
      % yields two multiplets at energies [-5/8, 3/8] = [-0.625, 0.375]
      % Wb,Feb17,16

        if isempty(Sfac) || Sfac==1, return; end

        wbdie('use standard normalization!');

        if norm(Sfac-sqrt(2))>1E-12
           wbdie('unexpected Sfac=%g !?',Sfac); end

        wblog(' * ','applying Sfac=sqrt(2) [Thomas Quella]');
        S=Sfac*S;

        S2=QSpace(contractQS(S,'13*',S,'13'))
        if isequal(S.Q{1},[0 2 0]) && ~isIdentityQS((1/12)*S2)
           wbdie('unexpected normalization of spin operator !?');
        end

     case 'SU6'
        if isempty(Sfac) || Sfac==1, return; end

        if norm(Sfac-sqrt(2))>1E-12
           wbdie('unexpected Sfac=%g !?',Sfac); end

        wblog(' * ','applying Sfac=sqrt(2) [Thomas Quella]');
        S=Sfac*S;

        S2=QSpace(contractQS(S,'13*',S,'13'))
        if isequal(S.Q{1},[0 0 2 0 0]) && ~isIdentityQS((1/24)*S2)
           wbdie('unexpected normalization of spin operator !?');
        end

     otherwise
   end

end

% -------------------------------------------------------------------- %

