function varargout=getLocalSpace(model,varargin)
% function [FF,..,Iout]=getLocalSpace(model [,'sym1,sym2,...',varargin])
%
%   build local model state space as specified by means of the typical
%   associated operators such as spin (S), annihilation (F), or
%   charge parity operator (Z).
%
%   The residual info structure Iout contains further operators
%   (if applicable), such as the identity operator (E), the spinor for
%   particle-hole symmetry (C3), its equivalent to the Casimir operator
%   S^2 (C2), or Q2 := (N-1)^2, with N the total particle number summed
%   over all channels, as required for isotropic Coulomb interaction.
%
% Models available
%
%    [F,Z,S,I]=getLocalSpace('FermionS',sym [,opts]);   % spinful fermions
%    [F,Z,  I]=getLocalSpace('Fermion', sym [,opts]);   % spinless fermions
%    [S,    I]=getLocalSpace('Spin', S      [,opts]);   % spin-S operators [SU(2)]
%    [S,    I]=getLocalSpace('SUN',N        [,opts]);   % SU(N) site
%    [S,    I]=getLocalSpace('SU<N>',       [,opts]);   % same as previous
%    [S,    I]=getLocalSpace('SON',5,       [,opts]);   % SO(N) site
%    [S,    I]=getLocalSpace('SO<N>',       [,opts]);   % SO(N) site
%    [S,    I]=getLocalSpace('Sp<2N>',      [,opts]);   % Sp(SN) site
%
% Options
%
%   'NC',..  number of channels (fermionic systems only)
%   '-A',..  use abelian symmetry (in 'Spin' mode only)
%   '-v'     verbose mode
%
% Symmetries for sym (single string, separated by commas)
%
%   'Acharge'      abelian total charge;         'Acharge(:)'   ^1)
%   'SU2charge'    SU(2) total particle-hole;    'SU2charge(:)' ^1)
%   'Aspin'        abelian total spin;           'Aspin(:)'     ^1)
%   'SU2spin'      SU(2) total spin (S);         'SU2spin(:)'   ^1)
%   'SU2spinJ'     SU(2) total spin (J=L+S)
%   'AspinJ'       U(1) total spin (J=L+S)_z
%   'SUNchannel'   SU(N) channel symmetry
%   'SONchannel'   SO(N) channel symmetry
%   'SpNchannel'   Sp(N) particle/hole (charge) * channel symmetry
%
%    ^1) sym(:) indicates to use given symmetry <sym> in the charge
%        or spin sector for each of the NC channels individually.
%
% Examples
%
%  % three channels with particle-hole SU(2) in each channel and total spin SU(2)
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge(:),SU2spin','NC',3,'-v');
%
%  % three channels with abelian charge, SU(2) spin and SU(2) channel
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SUNchannel','NC',3,'-v');
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SU3channel','-v'); % same
%
%  % three channels with SO(N) symmetry
%  % so far NC=3 only since this reduces to SU(2) with integer spins
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SO3channel','-v');
%    --> IS.L3 contains L=1 `spin opertors' in orbital space
%
%  % three channels with SU2spinJ symmetry (total J=L+S)
%    [FF,Z,JJ,IS]=getLocalSpace('FermionS','Acharge,SU2spinJ','NC',3,'-v');
%
%  % three spinless channels with SU(3) channel symmetry
%    [FF,Z,IS]=getLocalSpace('Fermion','SUNchannel','NC',3,'-v');
%
%  % single spin-S site
%    [S,IS]=getLocalSpace('Spin',1,'-v'); 
%
%    [S,IS]=getLocalSpace('SUN',3,'-v');  % same as ...
%    [S,IS]=getLocalSpace('SU3','-v');
%
% Examples with flavor groups / split channels // Wb,Feb16,18
%
%    Flavor groups NC(i) preserve their respective number of particles.
%    Hence this also splits the 'Acharge' symmetry, accordingly.
%    An SUNchannel subsymmetry is added only for every NC(i)>1 then.
%
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2spin,Acharge,SUNchannel','NC',[2 1],'-v');
%    [FF,Z,IS]=getLocalSpace('Fermion','Acharge,SUNchannel','NC',[2 1],'-v');
%
% Comment on symmetry label for U(1) charge:
%
%    The total charge label is always taken relative to half-filling.
%    Therefore in the spinless case if NC (or a group therein) has
%    an odd number n of flavors, then the total charge labels are
%    given by q=(2*n-1) to maintain integer q-labels.
%    getLocalSpace issues a NB/WRN in that respect.
%    This is similar to SU(2) or U(1) spin labels.
%
% Wb,Jul09,11 ; Wb,Jul30,12

% ---------------------------------------------------------------------- %
% for testing purposes
% * two channels with particle-hole and spin SU(2) in each channel
%      [FF,Z,SS,IS]=getLocalSpace(...
%      'FermionS','SU2charge(:),SU2spin(:)','NC',2,'-v');
% ---------------------------------------------------------------------- %

% fprintf(1,'\n'); eval(['help ' mfilename]);
  if nargin<1 || ~ischar(model)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  if nargin>1 && isequal(varargin{end},'-fixA')
       fixA=1; varargin=varargin(1:end-1);
  else fixA=0; end

  if ~isempty(regexp(model,'^S[OU]\d$'))
     varargin=[{str2num(model(3:end))}, varargin];
     model=[model(1:2) 'N'];
  elseif ~isempty(regexp(model,'^Sp\d$'))
     q=str2num(model(3:end)); if mod(q,2)
        wbdie('invalid symmetry ''%s''',model); end
     varargin=[{q/2}, varargin];
     model=[model(1:2) '2N'];
  end

  switch model
     case 'FermionS'
        n=4; varargout=cell(1,n);
        [varargout{:}]=getLocalSpace_SpinfullFermions(varargin{:});

     case 'Fermion'
        n=3; varargout=cell(1,n);
        [varargout{:}]=getLocalSpace_SpinlessFermions(varargin{:});

     case 'Spin'
        n=2; varargout=cell(1,n);
        if nargin>1 && ischar(varargin{1}) && regexp(varargin{1},'^SU\d')
             [varargout{:}]=getLocalSpace_SpinSUN(varargin{:});
        else [varargout{:}]=getLocalSpace_Spin(varargin{:});
        end

     case 'SUN'
        n=2; varargout=cell(1,n);
        [varargout{:}]=getLocalSpace_SUN(varargin{:});

     case 'SON'
        n=2; varargout=cell(1,n);
        [varargout{:}]=getLocalSpace_SON(varargin{:});

     case 'Sp2N'
        n=2; varargout=cell(1,n);
        [varargout{:}]=getLocalSpace_Sp2N(varargin{:});

     otherwise wbdie('invalid model ''%s''',model);
  end

  if     nargout<2, varargout=varargout(1);
  elseif nargout<n, varargout=varargout([1:nargout-1,end]);
  elseif nargout>n, wbdie(['invalid usage\n' ... 
     'requesting too many output arguments (%g/%g)'],nargout,n); end

  if fixA
  for io=1:numel(varargout), q=varargout{io}; m=0;
     if isa(q,'QSpace')
        if any(isAbelian(q)==2)
           q=fixAbelian(q); m=m+1;
        end
     elseif isstruct(q), ff=fieldnames(q);
        for j=1:numel(ff), s=getfield(q,ff{j});
           if isa(s,'QSpace') && any(isAbelian(s)==2)
              q=setfield(q,ff{j},fixAbelian(s)); m=m+1;
           end
        end
     end
     if m, varargout{io}=q; end
  end
  end

end

% -------------------------------------------------------------------- %
% centralize initialization of qfac and jmap!
% Wb,Apr04,14

function q=init_qstruct(istr,sym,N)

   if nargin<2 || nargin>3 || ~ischar(sym) || ~ischar(istr)
      wbdie('invalid usage (init_qstruct)'),
   end

   q.info=istr;
   q.type=[];

   q.Sp=SymOp;
   q.Sz=SymOp;

   q.qfac=[];
   q.jmap=[];

   switch sym
    case 'SU'
       if nargin<3, wbdie('invalid usage (missing N)'); end
       if numel(N)~=1 || N<2 || N>10
          wbdie('invalid symmetry ''%s(%g)''',sym,r);
       end

       r=N-1; if r>1
          q.jmap=[ diag(1./(1:r)) - diag(1./(2:r),-1) ]';
       end

       q.type=sprintf('%s%g',sym,N);
       q.qfac=2;

    case 'SO'

       if nargin<3, wbdie('invalid usage (missing N)'); end
       if numel(N)~=1 || N~=3, wbdie('invalid symmetry ''%s(%g)''',sym,r);
       end

       q.type='SU2';
       q.qfac=2;

    case 'Sp'

       if nargin<3, wbdie('invalid usage (missing N)'); end
       if numel(N)~=1 || N<2 || N>20 || mod(N,2)
          wbdie('invalid symmetry ''%s(2*%g)''',sym,N/2); end
       if N==2, wblog('WRN',...
         'got Sp(%g) which is equivalent to SU(2)',N); end

       r=N/2;

     % if r==3 % tags: Sp6 SP6
     %  % see MATH//tex-notes // SPN_LIE
     %  % see also clebsch.cc::findMaxWeight()!
     %    q.jmap=[6 0 0; -3 3 0; 0 -2 2]'/6; % NB! same as for SU3 above!
     % end

       if r>1
          q.jmap=[ diag(1./(1:r)) - diag(1./(2:r),-1) ]';
       end

       q.type=sprintf('%s%g',sym,N);

    otherwise

       if nargin>2, wbdie(...
         'invalid usage (got 3rd argument, having ''%s'')',sym); end
       q.type=sym;
   end

end

% -------------------------------------------------------------------- %
% 'FermionS'
% -------------------------------------------------------------------- %

function [F,Z,S,Iout]=getLocalSpace_SpinfullFermions(Sym_,varargin)
% function [F,Z,S,Iout]=getLocalSpace_SpinfullFermions(Sym,varargin)
% supported symmetries
%
%   'Acharge'      abelian total charge
%   'SU2charge'    SU(2) total particle-hole
%   'Aspin'        abelian total spin
%   'SU2spin'      SU(2) total spin
%   'SU2spinJ'     SU(2) total spin (J=L+S)
%   'AspinJ'       U(1) total spin (J=L+S)_z
%   'SUNchannel'   SU(N) channel symmetry
%   'SONchannel'   SO(N) channel symmetry
%   'SpNchannel'   Sp(2N) channel symmetry with SU(2) particle/hole
%   'SpNchannelS'  Sp(2N) channel symmetry with SU(2) spin
%
% Example: 'Acharge,SU2spin'
% tags: symplectic symmetry, Sp2n, Sp4, Sp6, SpNchannel
% Wb,Jul09,11

  if nargin<1 || ~ischar(Sym_)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('INIT',varargin); oc={};
     NC   =getopt('NC',[]);
     Yops =getopt('-Y');

     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2;
     else vflag=0; end

  getopt('check_error');

  check_unique_syms(Sym_);

  NC=check_NC(NC,Sym_);

  [NC,NCx,NCxs,NCsplit]=check_NC_vec(NC,Sym_);

  istr=sprintf('fermionic %s-level system',NCxs);

  f2=sparse([0 1; 0 0]);
  z2=f2*f2'-f2'*f2; n2=f2'*f2; e2=speye(2);

  FF=SymOp(NC,2); NN=SymOp(NC,2); SS=SymOp(NC,3);
  CP=SymOp(NC,1); CZ=SymOp(NC,1); n2=2*NC;

  sx=spinmat(2,'-sp','-sym');
  sx=sx([1 3 2]);

  for i=1:NC
     for j=1:2, k=2*(i-1)+j;
        q=[ repmat({z2},1,k-1), {f2}, repmat({e2},1,n2-k) ];
        FF(i,j)=SymOp(getOpName('fops',i,j),mkron(q{:}));
        NN(i,j)=FF(i,j)'*FF(i,j); NN(i,j).istr=sprintf('N%g',i);
     end

     k=2*i-1;
        q=[ repmat({e2},1,k-1), {mkron(f2',f2')}, repmat({e2},1,n2-k-1) ];
        CP(i)=SymOp(sprintf('C%g(+)',i),mkron(q{:})); q{k}=0.5*comm(q{k},q{k}');
        CZ(i)=SymOp(sprintf('C%g(z)',i),mkron(q{:}));

     for j=1:3
        q=[ repmat({e2},1,k-1), {blkdiag(0,sx{j},0)}, repmat({e2},1,n2-k-1) ];
        SS(i,j)=SymOp(getOpName('spin',i,j),mkron(q{:}));
     end
  end

  SOP=[]; D=dim(FF);
  phflag=0; chflag=0; SONflag=0; Jflag=0; Sflag=0;

  Sym=strread(Sym_,'%s','whitespace',';, ')'; clear q

  for sym1=Sym, sym1=sym1{1};
    n=str2num(regexprep(sym1,'(SU|SO)(\d+)channel','$2'));
    if ~isempty(n)
       if ~isequal(n,NC), wbdie(...
         'invalid symmetry %s (having NC=%g)',Sym_,NC); end
       sym1=[sym1(1:2) 'Nchannel'];
    end

    n=0;
    sym1=regexprep(sym1,'Sp(\d+)channel(?@n=str2num($1);)','');
    if n
       if ~isequal(n,2*NC), wbdie(...
         'invalid symmetry %s (having NC=%g)',Sym{i},NC); end
       sym1=['SpNchannel' sym1];
    end

    n=str2num(regexprep(sym1,'Z(\d+)charge','$1'));
    if ~isempty(n)
       if n<2, wbdie(...
         'invalid symmetry %s (having NC=%g)',Sym{i},NC); end
       sym1='ZNcharge'; zs=sprintf('Z%g',n); nz=n;
    end

    if isequal(sym1(end-2:end),'(:)')
         cflag=1; if NC==1, sym1=sym1(1:end-3); end
    else cflag=0; end

    switch sym1

      case 'Pcharge'
          q=init_qstruct('total charge parity','P');
          q.Sz=2*sum(CZ); z=diag(full(q.Sz.op)); z=z-min(z);
          e=norm(diff(unique(z))-1); if e
            wbdie('unexpected charge labels !?'); end
          z=1-2*mod(z,2);
          q.Sz=SymOp(['P[' q.Sz.istr ']'],diag(z),'-disc');

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'ZNcharge'
          q=init_qstruct(sprintf('total charge %s',zs),zs);
          q.Sz=2*sum(CZ); z=diag(full(q.Sz.op)); z=z-min(z);
          e=norm(diff(unique(z))-1); if e
            wbdie('unexpected charge labels'); end
        % NB! the symmetry operation must be unitary!
        % e.g. see getSymmetryOps/get_symmetries_op()
        % => unitary symmetry operation is given by exp(i*(2*pi/nz)*Sz)
        % where Sz contains actual symmetry labels (as required
        % for getSymStates below) // Wb,Aug29,16
          z=mod(z,nz);
          q.Sz=SymOp(['P[' q.Sz.istr ']'],diag(z),'-disc');

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'Acharge'

          ic=[0, cumsum(NCx)];

          for k=0:NCsplit
             i0=ic(k+1); i1=ic(k+2)-1;
             nc=i1-i0+1; ir=(i0:i1)+1;

             if NCsplit, s=sprintf(...
               'U(1) charge group #%g (%s)',k+1,vsprintf(ir,','));
             else s='total charge U(1)'; end

             q=init_qstruct(s,'A');
             q.Sz=2*sum(CZ(ir));
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'Acharge(:)'
          for i=1:NC
             q=init_qstruct(sprintf('charge U(1) [%d/%d]',i,NC),'A');
             q.Sz=2*CZ(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'Aspin'

          q=init_qstruct('total spin U(1)','A');
          q.Sz=2*sum(SS(:,3));
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'Aspin(:)'

          for i=1:NC
             q=init_qstruct(sprintf('spin U(1) [%d/%d]',i,NC),'A');
             q.Sz=2*SS(i,3);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'SU2charge'

          q=init_qstruct('total particle-hole SU(2)','SU',2);
          q.Sp=sum(CP);
          q.Sz=sum(CZ);
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          phflag=1;

      case 'SU2charge(:)'

          for i=1:NC
             q=init_qstruct(sprintf('particle-hole SU(2) [%d/%d]',i,NC),'SU',2);
             q.Sp=CP(i);
             q.Sz=CZ(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end
          phflag=2;

      case 'SU2spin'

          q=init_qstruct('total spin SU(2)','SU',2);
          q.Sp=sum(SS(:,1));
          q.Sz=sum(SS(:,3));
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'SU2spin(:)'
        % suggested by Geng-Dong Zhou / Seung-Sup Lee [email 02/28/2025]
		  for i=1:NC
             q=init_qstruct(sprintf('spin SU(2) [%d/%d]',i,NC),'SU',2);
             q.Sp=SS(i,1);
             q.Sz=SS(i,3);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q;end
		  end

      case {'SU2spinJ','AspinJ'}

          lx=spinmat(NC,'-sp','-sym');
          lx=lx([1 3 2]);

          LL=SymOp(1,3); JJ=SymOp(1,3);
          for p=1:3, Q=SymOp; q=lx{p}; [I,J]=find(q);
             for l=1:numel(I), i=I(l); j=J(l);
                 Q=Q+full(q(i,j))*(FF(i,1)'*FF(j,1) + FF(i,2)'*FF(j,2));
             end
             LL(p)=Q;
             JJ(p)=Q+sum(SS(:,p));
          end

          if sym1(1)=='A'
             s=sprintf('total U(1) spin (J=L+S)_z having 2l+1=NC=%g',NC);
             q=init_qstruct(s,'A'); Jflag=1;
             q.Sz=2*JJ(3);
          else
             s=sprintf('total SU(2) spin J=L+S having 2l+1=NC=%g',NC);
             q=init_qstruct(s,'SU',2); Jflag=2;
             q.Sp=JJ(1);
             q.Sz=JJ(3);
          end
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'SUNchannel'

          ic=[0, cumsum(NCx)];

          for k=0:NCsplit
             i0=ic(k+1); i1=ic(k+2)-1;
             nc=i1-i0+1;
             if nc<2, continue; end

             Sp=SymOp(1,nc-1); Sz=SymOp(1,nc-1);
             for i=1:nc-1, l=i0+i;
                Sp(i)=FF(l,1)'*FF(l+1,1) + FF(l,2)'*FF(l+1,2);
                Sz(i)=sum(CZ(i0+1:l))-i*CZ(l+1);
             end
             if NCsplit
                  s=sprintf(' group #%g',k+1);
             else s=''; end

             q=init_qstruct(sprintf(['SU(%g) channel' s],nc),'SU',nc);
             q.Sp=Sp;
             q.Sz=Sz;

             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

          chflag=1; clear Sp Sz

      case {'SpNchannel','SpNchannelS'}
         if sym1(end)~='S', phflag=12;
         else
            Sflag=1;
            if ~isempty(regexp(Sym_,'SU2charge')), wbdie(...
              'invalid symmetry setting (%s requires abelian charge)',sym1);
            end
         end
         chflag=2;

          if NC<2, wbdie('Sp(2*NC) symmetry requires NC>=2 (%g)',NC); end
          if NCsplit, wbdie('NC=(%s) not yet implemented for %s',NCxs,sym1); end

          if phflag
             psi=FF; for i=1:size(FF,1), psi(i,2)=psi(i,2)'; end
             psi=psi(:);
          else
             psi=FF(:);
          end

          sz=SymOp(1,NC); sp=SymOp(1,NC);

          for i=1:NC
             is=sprintf('Sp%g:z%g',2*NC,i);
             z=ones(1,NC); if i<NC, z(i+1)=-i; z(i+2:end)=0; end
             sz(i)=SymOp(is,spdiags([z,-z]',0,n2,n2));

             is=sprintf('sp%g:p%g',2*NC,i);
             if i<NC, s=sparse(i,i+1,1,NC,NC);
                  sp(i)=SymOp(is,blkdiag(s,-s'));
             else sp(i)=SymOp(is,sparse(i,i+NC,1,n2,n2)); end
          end

          if 1
             e=ones(1,NC); e(2:2:end)=-1;
             J=blkdiag(speye(NC),fliplr(spdiags(e',0,NC,NC)));
             sz=transform(sz,J);
             sp=transform(sp,J);

             psi=psi2psi(J,psi);
          end

          q=init_qstruct(sprintf('channel Sp(2*%g)',NC),'Sp',2*NC);
          q.Sp=mat2ops(sp,psi);
          q.Sz=mat2ops(sz,psi);

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'SONchannel'
          if NC~=3, wbdie(['invalid usage ' ... 
            '(got %s with NC=%g; only implemented for NC=3)'],sym1,NC);
          end

          Sp=SymOp(1,2); Sz=SymOp(1,2);
          for l=1:2, i=2*l-1;
             for s=1:2
                 Sp(l)=Sp(l)+FF(l,s)'*FF(l+1,s);
                 Sz(l)=Sz(l)+FF(i,s)'*FF(i,  s);
             end
          end
          Sp=Sp(1)+Sp(2);
          Sz=Sz(1)-Sz(2);

          q=init_qstruct(sprintf('SO(%g) channel',NC),'SO',nc);
          q.Sp=Sp;
          q.Sz=Sz;

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          clear Sp Sz

          SONflag=1;

      otherwise helpthis, wbdie('invalid symmetry %s (%s)',sym1,Sym_);
    end
  end

  check_commrels(SOP);
  SOP=set_qzvac(SOP);

  gotsym=~isempty(SOP);
  if gotsym
       sym=strhcat({SOP.type},'-s',',');
  else sym='(none)'; end

  if phflag
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(FF,NN,SS,CP,CZ,P);
  else P=[]; end

  if gotsym
     [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

     if Jflag, transform(JJ,LL,U); end

     SOP=transform_all(SOP,U); transform(FF,NN,SS,CP,CZ,U);
     SOP=resparse(SOP,1E-12);

     if ~isempty(P), U=P*U; end
  else
     U=speye(D); Is=[];
  end

  if nargout>1
     Iout=add2struct('-',NC,SOP,'sym=Sym_',U,Is);
  end

% --------------------------------------------------------------------
% F-operators

  if Jflag>1
       oP={'-P'};
  else oP={}; end

  X=FF; F=QSpace;
  if phflag, X=[FF,FF']; end

  if gotsym, k=0; 
     while ~isempty(X), k=k+1;
        [x,Io,X]=getSymmetryOps(X(1),SOP,X,oP{:});

        for i=1:numel(x), x{i}=full(x{i}); end
        F(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));

        check_ferm_acomm(F(k),numel(x));
     end

     if numel(F)==numel(FF), F=reshape(F,size(FF)); end
  else
     x=cell(size(X)); for i=1:numel(x), x{i}=full(X(i).op); end
     F=QSpace([0;0;0],cat(3,x{:}),{'','*','*'});

     check_ferm_acomm(F,numel(X));
  end

  if vflag,
     wblog('<i>','%s',istr); wblog('==>','%s [%s]',Sym_,sym);
     [s,i]=dsize(F(1),'-d');
     s=[ prod(i.sd(:,:),2), permute(prod(i.sc,2),[1 3 2])];
     s=[ sum(s(:)), sum(prod(s,2)) ];
     n=numel(F); if vflag>1, inl(1); end
     wblog(' * ','%g F-op%s @ z=%.3g',n,iff(n~=1,'s',''),s(1)/s(2));
     if vflag>1, info(F(1)); end
  end

  F=setOpFlag(F);

% --------------------------------------------------------------------
% EZ-operators

  x=sum(NN(:));
  if norm(x,'-offdiag')>1E-12, wbdie('N-operator not diagonal'); end
  x=full(diag(x.op));
  if norm(x-round(x))>1E-12, wbdie('N-operator is non-integer'); end
  x=round(x); oo={'E','Z'};

  E=eye(D);
  Z=full(diag((-1).^x));
  if gotsym
     z=[SOP.qzvac];
     E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, E));
     Z=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, Z));
  else
     E=QSpace([0;0],E,{'','*'});
     Z=QSpace([0;0],Z,{'','*'});
  end

  if nargout>1
     Iout.E=E;
     Iout.Z=Z;
  end

% --------------------------------------------------------------------
% spin-operators (for every channel and total)

  S=QSpace;
  if SONflag || chflag || phflag>10
     SS=sum(SS,1);
  end

  nc=NC/size(SS,1);
  if norm(nc-round(nc))>1E-12, wbdie(['failed to ' ... 
    'determine effect number of channels for splitup spin-operator']);
  end

  if Jflag
     S=get_spin_ops(JJ,SOP,Is,nc,sym,Jflag);
     Iout.Sop=get_spin_ops(sum(SS,1),SOP,Is,nc,sym,Jflag);
     Iout.Lop=get_spin_ops(LL,SOP,Is,nc,sym,Jflag);
  elseif Sflag
     [x,Io]=getSymmetryOps(SS(1),SOP);

     for i=1:numel(x), x{i}=full(x{i}); end
     S=QSpace(compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:})));
  else
     for i=1:size(SS,1)
        q=get_spin_ops(SS(i,:),SOP,Is,nc,sym);
        if i>1, S(i,:)=q; else S=q; end
     end
  end

  if SONflag
     i=findstrc({SOP.info},'SO');
        if numel(i)~=1, wbdie('invalid SOP->SON data'); end
     LX=[SOP(i).Sp, -SOP(i).Sp', SOP(i).Sz];
     Iout.L3=get_spin_ops(LX,SOP,Is,nc,sym);
  end

  if vflag
     n=numel(S); if vflag>1, inl(1); end
     wblog(' * ','%g S-op%s',n,iff(n~=1,'s',''));
     if vflag>1
       wblog(' * ','%s => %s\N',Sym_,sym); info(S(1));
     end
  end

  if ~Jflag
     if size(S,1)>1
          Iout.S3=S; S=sum(S,1);
     else Iout.S3=[]; end
  end

% --------------------------------------------------------------------
% operators for particle-hole symmetry

  if phflag

     n=numel(CZ); if phflag>10, n=1; end
     C2=QSpace(1,n); C3=QSpace(1,n); Q=QSpace;
     for i=1:n
        if phflag>10
             [x,Io]=getSymmetryOps(sum(CZ),SOP);
        else [x,Io]=getSymmetryOps(CZ(i),SOP);
        end
        x=mfull(x); x=cat(3,x{:}); x(find(abs(x)<1E-12))=0;

        C3(i)=compactQS(oc{:},sym, Is.QZ,Is.QZ, Io.QZ, x);
        C2(i)=contractQS(C3(i),'13*',C3(i),'13');
     end

     Iout.C3=C3;
     Iout.C2=C2;

     oo{end+1}=sprintf('%gxC3',n);
     oo{end+1}=sprintf('%gxC2',n);

     for i=1:n
     for j=1:n, Q=Q+contractQS(C3(i),'13*',C3(j),'13');
     end; end

     Iout.Q2=(4/3)*Q;
     % NB! C3(i) are proper C=1/2 `spin' operators
     % factor 4   from (1/2)^2 in definition of Cz
     % factor 1/3 since there are three ops in the IROP C3
     % NB! the above only works for one(!) channel
     % since for N>1 channel, channel symmetry and p/h
     % needs to be combined into Sp(2N); therefore e.g.
     % for N=2 channels, Q2 has non(!)-integer eigenvalues
     % Wb,Feb02,21

     oo{end+1}='Q2';
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); info(Z); end
  end

% --------------------------------------------------------------------
% C'C operators (Schrieffer-Wolff) // Wb,Oct03,12
% see Archive/getLocalSpace_160213.m for explicit construction of C'C,
% --------------------------------------------------------------------
% NB! if hopping preserves flavor, for N flavors this *always*
% reflects an SU(N) symmetry in the number of symmetric flavors
% Wb,Mar05,20

  if Yops, F1=F;
     if numel(F)>1, s=size(F);
        if any(s==1)
             wblog('NB!','using sum(F) for Yops'); F1=sum(F);
        else wblog('NB!','using sum(F(:,1)) for Yops'); F1=sum(F(:,1));
        end
     end

     Eo=getIdentity(F1,3); Zo=getIdentity(Eo,'-0');
     E2=getIdentity(Zo,2,Zo,1); E3=contract(Zo,'2*',E2,1);

     E2=getIdentity(F1);

   % NB! contracting E3 introduces an additional unitary
   % transformation on the operrtor space, which affects
   % the sign structure the irops in Y! however, since
   % Schrieffer-Wolff always comes in pairs (e.g. S.S)
   % this sign is irrelevant // Wb,Feb13,16
     Y=contract({F1,'1*',F1,1},'24',E3,'12');

     [y1,I,D]=uniquerows(Y.Q{3}); Y2=QSpace; YY=QSpace(1,0); nY=numel(I);
     for i=nY:-1:1
        Q=getsub(Y,I{i}); if norm(Q)>1E-12
           YY(end+1)=Q;
           Y2=Y2+contract(Q,'12*',Q,'12');
        end
     end

     Iout.Y =YY; nY=numel(YY);
     Iout.Y2=Y2;

   % NB! consider traceless operators only! (for consistency
   % with older version; see Archive/getLocalSpace_160213.m)
   % => Cij=Fi'*Fj -> for i==j, must subtract 0.5*Id
   %    this also ensures that |Cij|^2=4 for all i and j
   % => this corresponds to making the scaler term in Y traceless!
     if norm(YY(end).Q{3}(:))>1E-12, wbdie(...
       'invalid usage (misplaced scalar operator in Y ?!)'); end
     a=fixScalarOp(YY(end)); q=trace(a)/trace(E2);
     Iout.Y(end)=makeIrop(skipzeros(a-q*E2));

     if vflag
     wblog(' * ','%g Y-op%s (SW)',nY,iff(nY~=1,'s','')); end

  end

  if ~gotsym && nargout>1
     Iout.SOP=init_qstruct('no symmetry used','');
  end

end

% -------------------------------------------------------------------- %
% make sure there is only at most one occurance for symmetry types
% Wb,Apr04,18 // outsourced Wb,Mar25,25

function check_unique_syms(sym)

  for q={'spin','charge','channel'}, i=regexpi(sym,q{1});
     if numel(i)>1, wbdie(1,...
        'invalid symmetry ''%s''\nhaving ''*%s'' specified %s',...
        sym,q{1},int2str2(numel(i),'-#'));
     end
  end

end

% -------------------------------------------------------------------- %
% Wb,Apr04,18

function NC=check_NC(NC,Sym)

   ss={'SU','Sp','SO'}; Sym_=Sym; q=[]; 
   for i=1:numel(ss), n=0; 
      Sym=regexprep(Sym,[ss{i} '(\d+)channel(?@n=str2num($1);)'],'');
      if n, if i==2, n=n/2; end
         q(end+1)=n;
      end
   end

   if numel(NC)<=1, q=[q NC]; e=0;
      if isempty(q), q=1;
      elseif norm(q-round(q)) || any(q<1), e=1;
      elseif any(diff(q)), e=2; end
      if e, wbdie('invalid NC = %s (e=%d)',vec2str(q,'sep',' / '),e); end
      NC=q(1);
   elseif ~isempty(q)
      wbdie(['invalid usage (ambiguous Nchannel symmetry)\n' ... 
      'having ''%s'' vs. NC=[ %s ]'],Sym_,vec2str(NC,'-f'));
   elseif norm(NC-round(NC)) || any(NC<1)
      wbdie('invalid NC=[ %s ]',vec2str(NC,'-f'));
   end
end

% -------------------------------------------------------------------- %
% allow split-up of channels into lower symmetry configuration
% follow up to discussion with Fabian Kugler // Wb,Feb09,18

function [NC,NCx,NCxs,NCsplit]=check_NC_vec(NC,sym)

  if isempty(NC) || ~isvector(NC) || any(NC<1) || any(NC~=round(NC))
     NC, wbdie('invalid NC');
  elseif NC>6, wbdie('failed plausiblity check (got NC=%d !?)',NC); 
  end

  NCx=NC; NCsplit=numel(NC)-1;

  if NCsplit
     NC=sum(NC);
     NCxs=['[' vsprintf(NCx,'+') ']'];

     if ~isempty(regexpi(sym,'SU2charge'))
        wbdie('got NC=%s together with SU2charge',NCxs);
     end

     if ~isempty(regexpi(sym,'(SU|Sp)[N\d]*channel')) && NC<2+NCsplit
        if NCsplit, wbdie(1,[
           'NC=%s inconsistent with channel symmetry ''%s''\n' ...
           'note that non-abelian would require sum(NC)>%g'],NCxs,sym,NC);
        else wbdie('channel symmetry requires NC>=2 (%g)',NC); end
     end
  end

  if ~NCsplit, NCxs=sprintf('%g',NC); end

end

% -------------------------------------------------------------------- %

function S=get_spin_ops(X,SOP,Is,NC,sym,Jflag)
  if nargin<6, Jflag=0; end

  e=norm(comm(X(1),X(2))+X(3)); if e>1E-12
     wbdie('invalid (sparse) S-ops (%.3g)',e); end

  S=QSpace; k=0; X=fliplr(X);

  gotsym=~isempty(SOP);
  if gotsym
     while ~isempty(X), k=k+1;
        [x,Io,X]=getSymmetryOps(X(1),SOP,X);

        for j=1:numel(x), x{j}=full(x{j}); end
        S(k)=compactQS(sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
     end
  else
     x=cell(size(X)); for i=1:numel(x), x{i}=full(X(i).op); end
     S=QSpace([0;0;0],cat(3,x{:}),{'','*','*'});
  end

  s=zeros(0,3); for i=1:numel(S)
     q=getDimQS(S(1)); q=q(end,:); q(end+1:3)=1;
     s(i,:)=q;
  end
  if norm(diff(s(:,3)))>1E-12, wbdie('invalid S-ops'); end
  if s(1,3)==1
     i=reshape(flipud(reshape(1:numel(S),3,[])),size(S));
     S=S(i);

     if numel(S(1).Q)>2
        e=normQS(QSpace(...
          contractQS(S(2),'13*',S(2),'13')) ...
        - contractQS(S(1),'13*',S(1),'13') - S(3));
     else
        e=normQS(QSpace(...
          contractQS(S(2),'1*',S(2),1)) ...
        - contractQS(S(1),'1*',S(1),1) - S(3));
     end

     if e>1E-12, wbdie('invalid spin-CR'); end

     if ~Jflag
        q=QSpace; for i=1:3
          if numel(S(i).Q)>2
               q=q+contractQS(S(i),'13*',S(i),'13');
          else q=q+contractQS(S(i),'1*', S(i),'1' ); end
        end
        e=eigQS(q);
        if norm(max(e(:,1))-(NC/2)*(NC/2+1))>1E-12
        wbdie('invalid S-ops'); end
     end
  end

end

% -------------------------------------------------------------------- %
% tags: _spinless

function [F,Z,Iout]=getLocalSpace_SpinlessFermions(Sym_,varargin)
% function [F,Z,Iout]=getLocalSpace_SpinlessFermions(Sym,varargin)
% supported symmetries
%
%   'Acharge'      abelian total charge
%   'SUNchannel'   SU(N) channel symmetry
%
% NB! NO SU(2) particle-hole symmetry in the absence of spin!
% consequently, also no Sp(2n) symmetry
%
% Example: 'Acharge,SUNchannel'
% Wb,Jun27,12

  if nargin<1 || ~ischar(Sym_)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('INIT',varargin);
     NC=getopt('NC',[]); oc={};

     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end

     if     getopt('-q1'), zflag=0;
     elseif getopt('-q2'), zflag=2; else zflag=-1; end

  getopt('check_error');

  check_unique_syms(Sym_);

  NC=check_NC(NC,Sym_);

  [NC,NCx,NCxs,NCsplit]=check_NC_vec(NC,Sym_);

  istr=sprintf('spinless fermionic %s-level system',NCxs);

  f2=sparse([0 1; 0 0]);
  z2=f2*f2'-f2'*f2; n2=f2'*f2; e2=speye(2);

  FF=SymOp(NC,1); NN=SymOp(NC,1);

  for k=1:NC
     q=[ repmat({z2},1,k-1), {f2}, repmat({e2},1,NC-k) ];
     FF(k)=SymOp(sprintf('F%g',k),mkron(q{:}));
     NN(k)=FF(k)'*FF(k); NN(k).istr=sprintf('N%g',k);
  end

  SOP=[]; D=dim(FF);

  Sym=strread(Sym_,'%s','whitespace',';, ')'; clear q

  for s=Sym, s=s{1};
    n=str2num(regexprep(s,'SU(\d+)channel','$1'));
    if ~isempty(n)
       if ~isequal(n,NC), wbdie(...
         'invalid symmetry %s (having NC=%g)',s,NC); end
       s='SUNchannel';
    end

    n=str2num(regexprep(s,'Z(\d+)charge','$1'));
    if ~isempty(n)
       if n<2, wbdie(...
         'invalid symmetry %s (having NC=%g)',s,NC); end
       s='ZNcharge'; zs=sprintf('Z%g',n); nz=n;
    end

    if isequal(s(end-2:end),'(:)')
         cflag=1; if NC==1, s=s(1:end-3); end
    else cflag=0; end

    switch s
      case 'Pcharge'
          q=init_qstruct('total charge parity','P');
          q.Sz=sum(NN); z=diag(full(q.Sz.op)); z=z-min(z);
          e=norm(diff(unique(z))-1); if e
             wbdie('unexpected charge labels !?'); end
          z=1-2*mod(z,2);
          q.Sz=SymOp(['P[' q.Sz.istr ']'],diag(z),'-disc');

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'ZNcharge'
          if NCsplit, wbdie(...
            'multiple NC not yet implemented with %s',s); end

          q=init_qstruct(sprintf('total charge %s',zs),zs);
          q.Sz=sum(NN); z=diag(full(q.Sz.op)); z=z-min(z);
          e=norm(diff(unique(z))-1); if e
             wbdie('unexpected charge labels !?'); end
          z=mod(z,nz);
          q.Sz=SymOp(['P[' q.Sz.istr ']'],diag(z),'-disc');

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'Acharge'
          ic=[0, cumsum(NCx)];

          if zflag<0, zflag=any(mod(NCx,2)); end
          if zflag==1 && vflag
             if NCsplit, s=' groups with'; else s=''; end
             wblog('NB!','got odd number of flavors [n=%s]',NCxs);
             wblog(' * ',[
               'using 2*n-1 as charge labels for' s, 10 ...
               'odd number of flavors to ensure symmetric' 10 ...
               'integer labels around half filling\N'
             ]);
          end

          for k=0:NCsplit
             i0=ic(k+1); i1=ic(k+2)-1;
             nc=i1-i0+1; ir=(i0:i1)+1;

             if NCsplit, s=sprintf(...
               'U(1) charge group #%g (%s)',k+1,vsprintf(ir,','));
             else s='total charge U(1)'; end

             q=init_qstruct(s,'A');

             q.Sz=sum(NN(ir));
             if mod(nc,2) || zflag>1
                if zflag
                   q.Sz=2*q.Sz-nc;
                end
             else q.Sz=q.Sz-nc/2;
             end

             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'Acharge(:)'
          for i=1:NC
             q=init_qstruct(sprintf('charge U(1) [%d/%d]',i,NC),'A');
             q.Sz=NN(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'SUNchannel'

          ic=[0, cumsum(NCx)];

          if NC<2, wbdie('channel symmetry requires NC>=2 (%g)',NC); end

          for k=0:NCsplit
             i0=ic(k+1); i1=ic(k+2)-1;
             nc=i1-i0+1;
             if nc<2, continue; end

             Sp=SymOp(1,nc-1); Sz=SymOp(1,nc-1);
             for i=1:nc-1, l=i0+i;
                Sp(i)=FF(l)'*FF(l+1);
                Sz(i)=0.5*(sum(NN(i0+1:l))-i*NN(l+1));
             end
             if NCsplit
                  s=sprintf(' group #%g',k+1);
             else s=''; end

             q=init_qstruct(sprintf(['SU(%g) channel' s],nc),'SU',nc);
             q.Sp=Sp;
             q.Sz=Sz;

             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

          clear Sp Sz

      otherwise helpthis, wbdie('invalid symmetry %s (%s)',s,Sym_);
    end
  end

  check_commrels(SOP);
  SOP=set_qzvac(SOP);

  sym=strhcat({SOP.type},'-s',',');

  if 1
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(FF,NN,P);
  else P=[]; end

  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); transform(FF,NN,U);
  SOP=resparse(SOP,1E-12);

  if ~isempty(P), U=P*U; end

  if nargout>1, Iout=add2struct('-',NC,SOP,'sym=Sym_',U,Is); end

% -----------------------------------------------------------
% F-operators
% -----------------------------------------------------------

  X=FF;

  k=0; F=QSpace;
  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for i=1:numel(x), x{i}=full(x{i}); end
     F(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));

     check_ferm_acomm(F(k),numel(x));
  end

  if vflag,
     wblog('<i>','%s',istr); wblog('==>','%s [%s]',Sym_,sym);
     [s,i]=dsize(F(1),'-d');
     s=[ prod(i.sd(:,:),2), permute(prod(i.sc,2),[1 3 2])];
     s=[ sum(s(:)), sum(prod(s,2)) ];
     n=numel(F); if vflag>1, inl(1); end
     wblog(' * ','%g F-op%s @ z=%.3g',n,iff(n~=1,'s',''),s(1)/s(2));
     if vflag>1, info(F(1)); end
  end

  F=setOpFlag(F);

% -----------------------------------------------------------
% E,Q,Z-operators
% -----------------------------------------------------------

  x=sum(NN(:));
  if norm(x,'-offdiag')>1E-12, wbdie('N-operator not diagonal'); end
  x=full(diag(x.op));
  if norm(x-round(x))>1E-12, wbdie('N-operator is non-integer'); end
  x=round(x); oo={};

  z=[SOP.qzvac];

  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';
  Z=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, full(diag((-1).^x)))); oo{end+1}='Z';

  if nargout>1
     Iout.E=E;
     Iout.Z=Z;
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); info(Z); end
  end

end

% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_Spin(qloc,varargin)
% function [S,Io]=getLocalSpace_Spin(qloc,varargin)
% Options
%
%   -A            use U(1) abelian symmetry
%   Z2spin|--Z2   use Z2 spin parity symmetry only
%   --nosym       use no symmetries at all
%
% Wb,Jul30,12

  if nargin<1 || ~isnumber(qloc)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('INIT',varargin); oc={}; 
     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end

     sym=getopt('get_last',''); Aflag=0; Z2=0; noSym=0; istr='';
     if isequal(sym,'-A'), Aflag=1; istr='abelian';
     elseif isequal(sym,'--Z2') || isequal(sym,'Z2spin')
        Z2=1; istr='Z2';
     elseif isequal(sym,'--nosym'), noSym=1;
     elseif ~isempty(sym), sym, wbdie('invalid usage'); end

  getopt('check_error');

  if norm(round(2*qloc)-2*qloc) || qloc<=0
     wbdie('invalid spin S=%g',qloc); end

  istr=sprintf([istr ' spin S=%g system'],qloc);
  Iout=add2struct('-',istr,'Sloc=qloc');

  sx=spinmat(2*qloc+1,'-sp','-sym');
  sx=sx([1 3 2]);

  if Z2 && qloc==1/2
   % e.g. used for Kitaev model on honeycomb lattice, where bonds have
   % either XX or YY or ZZ interactions only => S_ and S+ must have the
   % same irop q-labels so they can be added up! Therefore in Young tableaus
   % parity implies 4m boxes for even, and 4m+2 boxes for odd parity
   % with m an integer // Wb,Oct25,21

     [F,Z,IS]=getLocalSpace('Fermion','Z2charge');
     sym=regexprep(IS.sym,'charge','spin');
     SOP=IS.SOP;
     SOP.info=sprintf('spin S=%s with Z2 parity only',rat2(qloc,'-s'));

     E=getIdentity(Z);
     Z.info.otype='operator';

     S=QSpace(1,3);
     S(1)=makeIrop(Z/2); % = Sz = comm(S(2),S(2)')
     S(2)=F/sqrt(2);
     S(3)=-S(2)'; S(3).info.itags{end}='*';

     q=sum(S); e=normQS(contractQS(q,'13*',q,'13') - 0.75*E);
     if e>1E-12, wbdie('got spin-half operator inconsistency @ %.3g',e); end

     Iout=add2struct(Iout,SOP,sym,'U=[]','Is=struct',E,'S2=[]');
     return

  elseif noSym
     SOP.info=sprintf('spin S=%s without symmetries',rat2(qloc,'-s'));
     SOP.type='';

     E=QSpace([0;0],eye(size(sx{1},1)),{'','*'});
     S=[
        QSpace([0;0],full(sx{3}),{'','*'}) ...
        QSpace([0;0],full(sx{2}),{'','*'}) ...
        QSpace([0;0],full(sx{1}),{'','*'})
     ];

     sym='';
     Iout=add2struct(Iout,SOP,sym,'U=[]','Is=struct',E,'S2=[]');
     return
  end

  SS=SymOp(1,3); clear q
  for j=1:3
     s=getOpName('spin',1,j);
     SS(1,j)=SymOp(s([1,3:end]),sx{j});
  end

  if Z2
     q=init_qstruct('total spin Z2','Z2');
     q.Sz=2*sum(SS(:,3));
     wbstop
  elseif Aflag 
     q=init_qstruct('total spin U(1)','A');
     q.Sz=2*sum(SS(:,3));
  else
     q=init_qstruct('total spin SU(2)','SU',2);
     q.Sp=sum(SS(:,1));
     q.Sz=sum(SS(:,3));
  end
  SOP=q; D=dim(SS);

  check_commrels(SOP);
  SOP=set_qzvac(SOP);

  sym=strhcat({SOP.type},'-s',',');

  if 1
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(SS,P);
  else P=[]; end

  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); transform(SS,U);
  SOP=resparse(SOP,1E-12);

  if ~isempty(P), U=P*U; end

  if nargout>1
     Iout=add2struct('-','Sloc=qloc',SOP,'sym=''SpinS''',U,Is);
  end

% -----------------------------------------------------------
% further operators: E
% -----------------------------------------------------------

  oo={}; z=[SOP.qzvac];
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';

  if nargout>1
     Iout.E=E;
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); end
  end

% -----------------------------------------------------------
% spin-operators
% -----------------------------------------------------------

  k=0; S=QSpace; X=sum(SS,1); oo={};

  e=norm(comm(X(1),X(2))+X(3)); if e>1E-12
     wbdie('invalid (sparse) S-ops (%.3g)',e);
  end

  X=fliplr(X);

  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for i=1:numel(x), x{i}=full(x{i}); end
     S(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  end

% NB! in case of SU(2), sign of S.data{1} may change depending
% on whether qloc <= 2 (=2S) or qloc > 2, since the underlying CGT
% is sorted, and hence may acquire a permutation; keep the sign
% from the above! // SUN_SIGN_SPIN_OP // Wb,Jul07,22
% in particular, in the defining irep, S.data{1}<0 (!)
% since S = [-S_+/sqrt(2),  S_z, S_-/sqrt(2) ] // SU2_SIGN_SPIN_OP
% this start with a *negative* coefficient, whereas the sign
% convention with CGTs is to start with a positive non-zero value

  S=fixMixedScalar(S);

  if qloc==1, oo{end+1}='and T-op';
     [x,Io]=getSymmetryOps(SS(1)*SS(1),SOP);
     for i=1:numel(x), x{i}=full(x{i}); end
     Iout.T=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  end

  if numel(S)==1
     q=contractQS(S,'13*',S,'13');
  else
     X=QSpace;
     for i=1:numel(S)
        if numel(S(i).Q)>2
             q=contractQS(S(i),'13*',S(i),'13');
        else q=contractQS(S(i),'1*', S(i),'1' ); end
        X=X+q;
     end
     q=X;
  end

  e=eigQS(q);
  if norm(max(e(:,1))-qloc*(qloc+1))>1E-12
  wbdie('invalid S-ops'); end

  if vflag
     Iout.S2=QSpace(q);
     n=numel(S); if vflag>1, inl(1); end
     wblog(' * ',['%g S-op%s ' oo{:}],n,iff(n~=1,'s',''));
     if vflag>1
       wblog(' * ','%s => %s\N',Sym_,sym); info(S(1));
     end
  end

end

% -------------------------------------------------------------------- %
% define spin-operator in defining representation through
% "Schrieffer-Wolff transformation" from fermionic setup
% via getLocalSpace('Fermion','SUNchannel','NC',N,...);
% then build TOTAL spin-operator until qloc is found
% which also ensures proper (=standard) normalization
% of the spin operator // NORM_CASIMIR

function [S,Iout]=getLocalSpace_SpinSUN(sym,qloc,varargin)
% function [S,Iout]=getLocalSpace_SpinSUN(sym,qloc [,opts])
%
%    qloc is the symmetry of the local multiplet to choose
%    (default: defining representation).
%
%    NB! qloc needs to be specified in terms of general SU(N)
%    notation based on the underlying Young diagram, hence
%    qloc=2*S in the case of SU(2) [note that this is different
%    from getLocalSpace_Spin which interprets qloc in terms
%    of standard half-spins].
%
% Wb,Apr07,14

   if nargin<2, helpthis, wbdie('invalid usage'); end
   if isempty(regexp(sym,'^SU\d+$'))
      helpthis, wbdie('invalid symmetry ''%s''',sym);
   end

   N=str2num(sym(3:end)); r=N-1;
   if ischar(qloc)
      q=regexprep(qloc,'^\s*\(*\s*(\d+)\s*\)*\s*','$1');
      if numel(q)~=r || isempty(regexp(q,'^\d+$'))
         wbdie('invalid qloc ''%s''',qloc);  end
      qloc=double(q-'0');
   end

   if ~isnumeric(qloc), helpthis, qloc
      wbdie('invalid qloc for symmetry ''%s''',sym);
   end

   getopt('INIT',varargin); oc={};
      if getopt('-v'), vflag=1;
      elseif getopt('-V'), vflag=2;
      elseif getopt('-q') || 1, vflag=0; end
   getopt('check_error');

   q0=[1, zeros(1,r-1)];

   if isempty(qloc), qloc=q0;
   elseif numel(qloc)~=r || size(qloc,1)~=1 || ...
     ~isequal(round(qloc),qloc) || any(qloc<0)
      qloc, wbdie('invalid qloc');
   end

   [~,~,I]=getLocalSpace('Fermion','SUNchannel','NC',N,oc{:});

   I.SOP(1).info=regexprep(I.SOP(1).info,'channel','Spin');

   if numel(q0)>1, q0(2,:)=flip(q0); end
   i=matchIndex(I.E.Q{1},q0); nS=size(q0,1);
   if numel(i)~=nS, wbdie('failed to identify def irep or its dual'); end
   q=getsub(I.E,i);

   Z=getIdentityQS(q,'-0');
   X=getIdentity(q,Z);
   X=contract(X,2,Z,'2*',[1 3 2]);
   i=find(sum((X.Q{1}-X.Q{2}).^2,2)==0 & sum(X.Q{3}.^2,2));
   S0=getsub(X,i);
   qs=S0.Q{3}(1,:); % qlabels for `spin'

   q2=q0(1,:); q2(end)=q2(end)+1;
   if size(S0.Q{3},1)~=nS || norm(diff(S0.Q{3},[],1)) || ~isequal(qs,q2)
      wbdie('failed to select spin operators'); 
   end

   istr=sprintf('%s spin-operator (q=%s)',sym,sprintf('%g',qs));

   S0=QSpace(norm_SUNspin_def(S0));

   n=6; Sk=S0;
   for k=1:n
      i=matchIndex(Sk.Q{1},qloc);
      if ~isempty(i)
         Sk=getsub(Sk,i); i=matchIndex(Sk.Q{2},qloc);
         Sk=getsub(Sk,i);
         break;
      end

      if k<n && (k==1 || dim(Ak,2)<1E6)
         if k>1
            Ak=QSpace(permuteQS(getIdentityQS(Ak,2,S0,1),[1 3 2]));
            Q=contractQS(Ak,1,Sk,2);  % RsL'o
         else
            Ak=QSpace(permuteQS(getIdentityQS(S0,1,S0,1),[1 3 2]));
            Q=contractQS(Ak,1,S0,2);  % RsL'o
         end

         Sk=contract(Ak,'13*',Q,'32') + ...  % RR'o
            contract(Ak,'13*',contractQS(Ak,3,S0,2),'13');
      else
         wbdie('failed to find qloc=(%s) within %g iterations',...
         sprintf('%g',qloc),k);
      end
   end
   S=skipzeros(Sk); S.info.otype='operator';

   for i=1:numel(S.data)
      q=S.data{i}; s=size(q); s(end+1:3)=1; e2=eye(s(1:2)); e=0;
      for j=1:prod(s(3:end)), qj=q(:,:,j);
          e=e+norm(qj-qj(1)*e2);
      end
      if e>1E-12
         wbdie('got unexpected spin operator (e=%.3g)',e);
      else
         s(1:2)=1;
         S.data{i}=reshape(q(1,1,:),s);
      end
   end

 % ------------------------------------------------------------------- %
 % if numel(S.data{1})==1 && S.data{1}<0
 %    S=-S; % Wb,Jul08,16  //  ** DON'T! // Wb,Jul06,22
 % end
 % ------------------------------------------------------------------- %
 % NB! once the sign in S0 is fixed [cf. norm_SUNspin_def(S0) above]
 % this fixes the sign of the spin operator for all other multiplets!
 % E.g. consider the fusion of two single multiples, say SU2 multiplets
 % S1 and S2; then the scalar operator S1.S2 in the combined state
 % space via identity A-tensor will have a value that depends on
 % the fused multiplet Stot considered since
 %    S1.S2 = [ (S1+S2)^2 - S1^2 - S2^2 ] / 2    (*)
 % here an independent sign convention on S1 and S2 individually
 % will result in an arbitrary sign of S1.S2 which quickly leads
 % to inconsistencies, [e.g. for SU2, S=3/2 the above would flip
 % the sign, since the underlying CGT picks up a permutation
 % with S=3/2 > 1 = spin of adjoint]. Yet, the sign is well-defined
 % on physical grounds! Conversely, building the total spin operator
 % by combining Stot = S1 + S2, a wrong sign in (*) can give
 % wrong quadratic Casimirs!
 %
 % When having local state space comprised of several multiplets q,
 % the spin operator for that combined state space may be built via
 %     S0 = \sum_q getLocalSpace('Spin','SUN',q);
 % which must have consistent sign across the multiplets to avoid
 % inconsistencies // SUN_SIGN_SPIN_OP // Wb,Jul06,22
 % ------------------------------------------------------------------- %

   if nargout<2, return; end

   Iout.istr=istr;
   q=add2struct(rmfield(I,{'Z','NC'}),qloc,qs);

   Iout=structmerge(Iout,q);
   Iout.sym=sprintf('Spin%s',sym);

   d=getDimQS(S); d=d(end,:);

   if vflag
      wblog('<i>','setup for %s d=%d',istr,d(end));
      wblog(' * ','having qloc=(%s) d=%d', sprintf('%g',qloc),d(1));
   end

   Iout.Ef=Iout.E;

   Iout.E=QSpace(contractQS(S,'13*',S,'13'));
   Iout.E.data{1}=1;

end

% -------------------------------------------------------------------- %
% get spin operator in the defining representation only
% for a more generic setting, see getLocalSpace_SpinSUN.m above
% -------------------------------------------------------------------- %

function [S0,Iout]=getLocalSpace_SUN(N,varargin)
% function [S0,Io]=getLocalSpace_SUN(N,varargin)
%
%    spin operator for plain SU(N) site in the defining
%    representation with single local IROP transforming
%    according to the adjoint representation.
%
% Wb,Nov03,12

  if nargin<1 || ~isnumber(N)
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('INIT',varargin); oc={};
     if getopt('-v'), vflag=1;
     elseif getopt('-V'), vflag=2; else vflag=0; end
  getopt('check_error');

  istr=sprintf('SU(%g) site',N);

  if N<2 || N>8, wbdie('SU(%g) not yet implemented !?',N); end

  Sp=SymOp(1,N-1); Sz=SymOp(1,N-1);
  for i=1:N-1
     Sp(i)=SymOp(sprintf('Sp%g',i),sparse(i,i+1,1,N,N));
     z=ones(N,1); z(i+1)=-i; z(i+2:end)=0;
     Sz(i)=SymOp(sprintf('Sz%g',i),spdiags(z/2,0,N,N));
  end

  q=init_qstruct(sprintf('plain SU(%g) site',N),'SU',N);
  q.Sp=Sp;
  q.Sz=Sz;

  SOP=q; D=dim(Sp); clear Sp Sz

  check_commrels(SOP);
  SOP=set_qzvac(SOP);

  sym=strhcat({SOP.type},'-s',',');

  if 1
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P);
  else P=[]; end

  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U);
  SOP=resparse(SOP,1E-12);

  if ~isempty(P), U=P*U; end

  if nargout>1
     Iout=add2struct('-',N,SOP,sprintf('sym=''SU%g''',N),U,Is);
  end

  if vflag
  wblog('<i>','setting up %s',istr); end

% -----------------------------------------------------------
% get operators
% -----------------------------------------------------------

  oo={}; z=[SOP.qzvac];
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';

  if nargout>1
     Iout.E=E;
  end

  if vflag
     s=iff(numel(oo)>1,'s','');
     wblog(' * ','got { %s } op%s',strhcat(oo,'sep',', '),s);
  end

  [x,Io]=getSymmetryOps(SymOp('MW',sparse(1,N,1,N,N)),SOP);
  for i=1:numel(x), x{i}=full(x{i}); end

  S0=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  S0=QSpace(norm_SUNspin_def(S0));

  if vflag
     wblog(' * ','got IROP for adjoint representation');
  end

end

% -------------------------------------------------------------------- %
% get spin operator in given representation // Wb,Mar24,18
% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_SON(N,qloc,varargin)
% function [S,Io]=getLocalSpace_SON(N,qloc,varargin)
%
%    plain SO(N) site with single local IROP,
%    corresponding to the adjoint represenation.
%
% Wb,Mar24,18

  if mod(N,2)
     if N==5, qadj='02';
     elseif N>5, qadj=['01' repmat('0',1,(N-5)/2)];
     wbdie('invalid usage [SO(%g) ?]',N); end
  else
     if N==6, qadj='011';
     elseif N>6, qadj=['01' repmat('0',1,(N-4)/2)];
     else wbdie('invalid usage [SO(%g) ?]',N); end
  end

  sym=sprintf('SO%g',N);
  if isnumeric(qloc), qloc=sprintf('%g',qloc); end

  [q2,i2]=sort({qloc,qadj});
  c=SymStore(sym,'-C',[q2{1} ',' q2{2} ';' qloc]);
  S=QSpace(c);
  S=conj(QSpace(permuteQS(S,[3 i2])));
  S.info.otype='operator';

% normalize spin operator
% make sure to (also) have spin operator in defining irep (which defines f below)
% NB! by convention, each generator (spin operator) in the defining
% representation (10*) has |S_alpha|^2 = 2
% which is consistent e.g. with SO(3): S=1 SU(2) spin irep

  qdef=['1' repmat('0',1,numel(qloc)-1)];
  if isequal(qloc,qdef)
     Sdef=S;
  else
     [q2,i2]=sort({qdef,qadj});
     c=SymStore(sym,'-C',[q2{1} ',' q2{2} ';' qloc]);
     Sdef=QSpace(c);
     Sdef=conj(QSpace(permuteQS(Sdef,[3 i2])));
     Sdef.info.otype='operator';
  end

  d=getDimQS(Sdef); da=d(end);
  Sdef=(sqrt(da/2)/normQS(Sdef))*Sdef; % standard `spin' normalization!

  ar=QSpace(contractQS(contractQS(S,2,S,1),'13',S,'12*'));
  ad=QSpace(contractQS(contractQS(Sdef,2,Sdef,1),'13',Sdef,'12*'));

% making use of (adj x adj -> adj) which must always exist
% NB! the constant f below represents the norm of the structure constants
% and hence must be independent(!) of the irep (sic!)
% (cf. PHYS notes -> ref. Peskin) // Wb,Apr20,18
% => e.g. determine f (once and for all) in the defining irep
  if 1
     f=normQS(ad)/normQS(Sdef)^2;
  else
     c=SymStore(sym,'-C',[qadj ',' qadj ';' qadj]);
     a=QSpace(c);
     f=getscalar(QSpace(contractQS(a,'123',ad,'123')));
  end

  x=f*normQS(S)^2/normQS(ar);
  S=x*S;

  Iout.istr=sprintf('%s in (%s) [a=%s]',sym,qloc,qadj);
  Iout.qloc=S.Q{1};
  Iout.sym=sym;
  Iout.E=QSpace(getIdentityQS(S,1));
  Iout.SOP.info='';
  Iout.SOP.type='';

end

% -------------------------------------------------------------------- %
% get spin operator in given representation // Wb,Mar24,18
% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_Sp2N(N,qloc,varargin)
% function [S,Io]=getLocalSpace_Sp2N(N,qloc,varargin)
%
%    plain Sp(2N) site with single local IROP,
%    in the adjoint represenation.
%
% Wb,May14,18

  if N<2 || N~=round(N), wbdie('invalid Sp(2*%g) symmetry',N); end

  qadj=['2' repmat('0',1,(N-1))];

  sym=sprintf('Sp%g',2*N); tflag=0;
  if isnumeric(qloc), qloc=sprintf('%g',qloc);
  elseif ischar(qloc)
     if regexp(qloc,'--*test-*(.*)(?@tflag=$1;)')
        qloc=['1' repmat('0',1,N-1)];
        if isempty(tflag), tflag=1; end
     end
  end
  if isempty(regexp(qloc,'^[0-9A-Z]+$'))
     qloc, wbdie('invalid usage (qloc)');
  end

  [q2,i2]=sort({qloc,qadj});
  c=SymStore(sym,'-C',[q2{1} ',' q2{2} ';' qloc]);
  S=QSpace(c);
  S=conj(QSpace(permuteQS(S,[3 i2])));
  S.info.otype='operator';

% normalize spin operator
% make sure to (also) have spin operator in defining irep (which defines f below)
% NB! by convention, each generator (spin operator) in the defining
% representation (10*) has |S_alpha|^2 = 2
% which is consistent e.g. with SO(3): S=1 SU(2) spin irep

  qdef=['1' repmat('0',1,numel(qloc)-1)];
  if isequal(qloc,qdef)
     Sdef=S;
  else
     [q2,i2]=sort({qdef,qadj});
     c=SymStore(sym,'-C',[q2{1} ',' q2{2} ';' qloc]);
     Sdef=QSpace(c);
     Sdef=conj(QSpace(permuteQS(Sdef,[3 i2])));
     Sdef.info.otype='operator';
  end

  d=getDimQS(Sdef); da=d(end);
  Sdef=(sqrt(da/2)/normQS(Sdef))*Sdef; % standard `spin' normalization!

  ar=QSpace(contractQS(contractQS(S,2,S,1),'13',S,'12*'));
  ad=QSpace(contractQS(contractQS(Sdef,2,Sdef,1),'13',Sdef,'12*'));

% making use of (adj x adj -> adj) which must always exist
% NB! the constant f below represents the norm of the structure constants
% and hence must be independent(!) of the irep (sic!)
% (cf. PHYS notes -> ref. Peskin) // Wb,Apr20,18
% => e.g. determine f (once and for all) in the defining irep
  if 1
     f=normQS(ad)/normQS(Sdef)^2;
  else
     c=SymStore(sym,'-C',[qadj ',' qadj ';' qadj]);
     a=QSpace(c);
     f=getscalar(QSpace(contractQS(a,'123',ad,'123')));
  end

  x=f*normQS(S)^2/normQS(ar);
  S=x*S;

  Iout.istr=sprintf('%s in (%s) [qadj=%s]',sym,qloc,qadj);
  Iout.qloc=S.Q{1};
  Iout.sym=sym;
  Iout.E=QSpace(getIdentityQS(S,1));
  Iout.SOP.info='';
  Iout.SOP.type='';

  if ~tflag, return; end

  sp=[0 1; 0 0];   s1=diag([1 0]); ir2=1/sqrt(2);
  sz=diag([1 -1]); s2=diag([0 1]);

  SOP=init_qstruct('TST symplectic root space','Sp',2*N);
  SOP.Sz=SymOp(1,N);
  SOP.Sp=SymOp(1,N);

  for i=1:N
     z=zeros(1,N); z(i)=ir2;
     SOP.Sz(i)=SymOp(getOpName('Sz',i),mkron(diag(z),sz));

     if i>1
        q = mkron(sparse(i,i-1,ir2,N,N),s1) ...
          - mkron(sparse(i-1,i,ir2,N,N),s2);
     else
        q = mkron(diag([1 repmat(0,1,N-1)]),sp);
     end
     SOP.Sp(i)=SymOp(getOpName('Sp',i),q);
  end

  check_commrels(SOP);
  Iout.SOP=set_qzvac(SOP);

  [x,Io]=getSymmetryOps(SOP.Sp(1),SOP,'-t');
  Iout.TT=x;
  Iout.IT=Io;

  if isequal(tflag,'fig') && (N==2 || N==3)
     ah=smaxis(1,1,'tag',mfilename); addt2fig Wb
     header('%M');

     qz=Io.qz{1}; s=size(qz);
     qz(:,:,2)=0;
     qz(:,:,3)=nan;
     dd=reshape(permute(qz,[3 1 2]),3*s(1),s(2));

     if N==2
          plot (dd(:,1),dd(:,2),        'o-');
     else plot3(dd(:,1),dd(:,2),dd(:,3),'o-');
     end
     sms(6); axis equal
  end
end

% -------------------------------------------------------------------- %
% keep scalar operator dimension (3rd index)
% if non-scalar operators exist within the same QSpace array

function F=fixMixedScalar(F)

   if numel(F)<2, return; end

   sF=size(F);
   mark=zeros(sF); n=prod(sF); otype=cell(sF);
   for i=1:n, r=numel(F(i).Q); otype{i}=F(i).info.otype;
      if r==3, q=norm(F(i).Q{3});
         if q>1E-12, mark(i)=1; else mark(i)=-1; end
      elseif r~=2
         wbdie('invalid usage (got r=%g !?)',r);
      end
   end

   if any(mark(:)>0)
      I=find(mark==0); n=numel(I);
      if n &&~isempty(otype), i=cellfun(@isempty,otype);
         otype=unique(otype(find(~i))); i=numel(otype);
         if     i==0, otype='';
         elseif i==1, otype=otype{1};
         else otype, wbdie('got mixed otype setting'); end
      end
      for k=1:n, i=I(k);
          F(i)=makeIrop(F(i));
          if otype, F(i).info.otype=otype; end
      end
   else
      I=find(mark<0); n=numel(I);
      for i=1:n, F(i)=squeeze(F(i)); end
   end
end

% -------------------------------------------------------------------- %

function s=getOpName(tag,varargin)

   switch tag
      case 'fops'

         if numel(varargin)~=2, wbdie('invalid usage (%s)',tag); end
         [i,s]=deal(varargin{:});

         switch s
            case 1, s='u';
            case 2, s='d';
            otherwise s, wbdie('invalid spin');
         end
         s=sprintf('F%g%s',i,s);

      case 'spin'
         if numel(varargin)~=2, wbdie('invalid usage (%s)',tag); end
         [i,s]=deal(varargin{:});

         switch s
            case 1, s='+';
            case 2, s='-';
            case 3, s='z';
            otherwise s, wbdie('invalid spin');
         end
         s=sprintf('S%g(%s)',i,s);

      case { 'Sp', 'Sz' }
         if numel(varargin)~=1, wbdie('invalid usage (%s)',tag); end
         s=[tag '_' num2str(varargin{1})];

      otherwise, tag, wbdie('invalid tag'); 
   end
end

% -------------------------------------------------------------------- %

function F=setOpFlag(F)

   for i=1:numel(F)
      if numel(F(i).Q)==3,F(i).info.otype='operator'; end
   end

end

% -------------------------------------------------------------------- %

function SOP=transform_all(SOP,U)

   s=inputname(2); if ~isempty(s), o={'-n',s}; else o={}; end
   for i=1:numel(SOP)
      SOP(i).Sp=transform(SOP(i).Sp,o{:},U);
      SOP(i).Sz=transform(SOP(i).Sz,o{:},U);
   end

end

function SOP=resparse(SOP,varargin)

   for i=1:numel(SOP)
      SOP(i).Sp=sparse(SOP(i).Sp, varargin{:});
      SOP(i).Sz=sparse(SOP(i).Sz, varargin{:});
   end

end

function SOP=strip4mex(SOP)

   for i=1:numel(SOP)
      SOP(i).Sp=getops(SOP(i).Sp,'-f');
      SOP(i).Sz=getops(SOP(i).Sz,'-f');
   end

end

% -------------------------------------------------------------------- %
% Wb,Dec03,11

function G=mat2ops(G,psi)

  s=size(psi(1)); q0=sparse(s(1),s(2));

  for k=1:numel(G)
     [i,j,g]=find(G(k).op); q=SymOp('',q0);
     for l=1:numel(g), q=q+g(l)*psi(i(l))'*psi(j(l)); end
     G(k)=q;
  end

end

function Upsi=psi2psi(U,psi)

  if numel(psi)~=size(U,1), wbdie(...
    'size mismatch (%g/%g)',size(U,1),numel(psi)); end
  e=norm(U*U'-speye(size(U)),'fro'); if e>1E-12, wbdie(...
    'unitary matrix U expected! (%.3g)',e); end

  s=dim(psi); Upsi=SymOp(size(psi)); q0=sparse(s,s);

  for i=1:numel(psi)
     [x,j,g]=find(U(i,:)); q=SymOp('',q0);
     for l=1:numel(g), q=q+g(l)*psi(j(l)); end
     Upsi(i)=q;
  end
end

% -------------------------------------------------------------------- %
% check all commutators for symmetry operations (within sparse still)
% outsourced from getLocalSpace_SpinfullFermions()
% Wb,Jul31,12

function check_commrels(SOP)

  for i=1:numel(SOP)
     S1=SOP(i).Sz; n1=numel(S1);
     for i1=1:n1, for i2=i1+1:n1
        e=norm(comm(S1(i1),S1(i2)));
        if e>1E-12, wbdie('invalid z-operator'); end
     end, end

     S1=[ SOP(i).Sz, SOP(i).Sp ]; n1=numel(S1);

     for j=i+1:numel(SOP)
        S2=[ SOP(j).Sz, SOP(j).Sp ]; n2=numel(S2);
        for i1=1:n1, for i2=1:n2
           e=norm(comm(S1(i1),S2(i2)));
           if e>1E-12, inl(1); disp(SOP(i)), disp(SOP(j))
              wbdie(['got non-commuting symmetries ' ...
             'SOP(%d).%d and SOP(%d).%d (e=%g)'],i,i1,j,i2,e);
           end
        end, end
     end
  end
end

% -------------------------------------------------------------------- %
% safeguard: check fermionic commutator relations
% outsourced // Wb,Nov04,16

function check_ferm_acomm(F,nops)

  for k=1:numel(F), rk=numel(F(k).Q);
     if     rk==3, ic={'13*','13'; '23','23*'};
     elseif rk==2, ic={'1*', '1' ; '2', '2*' };
          wblog(1,'WRN','got fermionic F-operator of rank %d',rk);
     else wbdie('unexpected fermionic F-operator of rank %d',rk);
     end

     a=contract(F(k),ic{1,1},F(k),ic{1,2});
     b=contract(F(k),ic{2,1},F(k),ic{2,2});
     Q=a+b; n=Q.data{1}(1); n(2)=max(1,round(n));

     e=abs(diff(n));
     if nargin<2, nops=n(2); else e(2)=abs(n(1)-nops); end

     if any(e>1E-12), wbdie(1,...
        'invalid fermionic creation/annihilation ops (n=%g)',n(1)); end

     Q=(1/nops)*Q;
     if ~isIdentityQS(Q), wbdie(1,...
        'invalid fermionic creation/annihilation ops');
     end
  end
end

% -------------------------------------------------------------------- %
% Wb,Nov04,16

function [SOP,qzvac]=set_qzvac(SOP)

  for i=1:numel(SOP)
     if isequal(SOP(i).type,'P')
          SOP(i).qzvac=+1;
     else SOP(i).qzvac=zeros(1,numel(SOP(i).Sz)+numel(SOP(i).Sp));
     end
  end

end

% -------------------------------------------------------------------- %
% Enforce standard normalization of Casimir in the DEFINING irep
% e.g. see Peskin & Schroeder, and also PHYS notes Eq.(25.16)
%
%    Casimir for SU(N) in the defining irep: S2=(N^2-1)/2N
%       SU(2)   spin-1/2   S2 = 3/4  ok. Casimir(sym,q)
%    Casimir for SU(N) in the adjoint irep: S2=N
%       SU(2)   spin-1     S2 = 2    ok. Casimir(sym,q)
%    for a check on the defining and adjoint reps for other groups
%
% see MEX//Casimir.m -> SymStore(sym,'Casimir',q);
% tags: NORM_CASIMIR // Wb,Dec03,15

function S=norm_SUNspin_def(S)

 % this routine may be called for the spin operator in the defining
 % AND its dual representation; hence S.data can have multiple entries,
 % even if there cannot be any OM in the defining or its dual irep.
 % Enforce that the entries must be identical in size and norm though
 % Wb,Jul07,22

   nS=numel(S.data); sym=S.info.qtype;
   r=size(S.Q{1},2); q0=[1, zeros(1,r-1)];
   if nS>1
      [~,is]=sortrows(S.Q{1},'descend');
      S=getsub(S,is);
      q0(2,:)=flip(q0);
   end
   if ~isequal(S.Q{1},q0) || ~isequal(S.Q{1},S.Q{2}) || ...
       isempty(regexp(sym,'^SU\d+$')), S
       wbdie('unexpected Sop input');
   end
   if nS>1
      q=S.data; q=permute(cat(4,q{:}),[4 3 1 2]);
      if abs(q(1))<1E-12, q=fliplr(q); end
      if norm(diff(size(q))) || norm(q-q(1)*eye(size(q)))
         wbdie('unexpected input spin operators'); 
      end
   end

   cgr=S.info.cgr; cgw=[cgr.cgw];
   if size(cgr,2)~=1 || norm(diff(cgw(:)))>1E-12
      wbdie('invalid CGR data'); end
   cgw=cgw(1);

   for k=1:nS, S.data{k}=1; end

   d=getDimQS(S); d=d(end,:);
   da=d(3);

   q=sqrt(da/2)/cgw;

   if isequal(sym,'SU2'), if nS>1, wbdie('invalid usage'); end
    % standard conventions imply that the reduced matrix element of the
    % SU(2) operator in the defining irep is *negative*, having
    % S = [-S_+/sqrt(2), S_z, S_-/sqrt(2) ]. Note that this *is*
    % already the case with getLocalSpace_Spin()! Thus for consistency,
    % change sign also to negative here // SU2_SIGN_SPIN_OP // Wb,Jul07,22
      q=-q;
   end

   for k=1:nS, S.data{k}=q; end

   if nS>1
    % WRN! the spin operator for the dual of the defining irep may have
    % a first non-zero first value that is negative! (it actually does
    % for SU3!). Note that the spin operator in the dual irep may be
    % obtained from the total spin operator by tensor product of N-1
    % defining ireps; this spin operator is well-defined with no! freedom
    % to choose a sign! (while the sign of the input to this routine was
    % chosen rather arbitrarily!). // Wb,Jul08,22

      S1=getsub(S,1); S2=getsub(S,2);
      A=getIdentity(S1,S2,[1 3 2]);
      S12=contract(A,'13*',contractQS(S1,'13*',contractQS(A,3,S2,2),'14'),'13');
      i=find(sum(S12.Q{1}.^2,2)==0);
      if S12.data{i}>0
         S.data{1}=-S.data{1};
      end
   end

end

% -------------------------------------------------------------------- %

