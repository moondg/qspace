function [X3,Iout]=get_Perm3(E0,varargin)
% function [X3,Iout]=get_Perm3(E0 [,idx])
%
%    Get 3-site permuation operator set, where the returned X3
%    contains the three operators for the equivalent operation of
%    acting with these operators X3(1..3) on sites 1,2,3 in that order.
%
%    X3 is the center operator that describes a 3-site permutation.
%    It is rank-4 with two trailing irop indices, to be contracted
%    across 3-steis also using Iout.EX as in EX' * X3 * EX
%    (see also discussion for chi below). The inverse permutation
%    is returned in Iout.X3i (see below).
%
% The output structure Iout contains various QSpaces
%
%    .E0   identity operator of a given site (based on the input)
%    .EX   complete local operator basis
%    .X2   (X2' * X2) fully describes to the 2-site permuation,
%          i.e., the swap operator
%    .X3i  inverse operation to X3 (first argument returned),
%          having X3i=X3' (but the dagger is hidden due to the
%          rank-4 structure of X3.
%
%    .chi  chiral operator built from clockwise and counter-clockwise
%          permutation: chi = (X3-X3')/4i (for --chi-real: 1/4i -> 1/4)
%          For a spin-half Hamiltonian, this is the regular S.(SxS).
%          However, e.g. for a Hubbard site (d=4), this also includes
%          other terms, as chi corresponds to the full chiral term
%          as for an SU(d) state spaces, which acquires substructure
%          if a smaller symmetry is used. This then splits in various
%          contributions, e.g. the chiral operator in the charge channel
%          (pseudo-spin, etc!). See option wsys to constrain chi, in
%          which case the full chiral operator is store in Iout.CHI.
%
%          Eventually, the chiral operator needs to be applied 
%          across three different sites as follows: EX' * chi * EX 
%
% Options
%
%    -t          test mode
%    -k          stop with keyboard before returning
%
%   --chi-real   do not include 1i factor in the chiral operator Iout.chi
%                i.e., keep its matrix elements real (since e.g. for chiral
%                correlation, <chi'_i * chi_j>, this phase won't matter).
%
%   'wchi',{..}  reduce chiral operator to subspace in case of larger 
%                local state space (by, default, full chiral operator
%                is returned). This also reduces EX to the relevant
%                irop space for chi, returned in Iout.Ex.
%
% Examples:
%
%   'wchi',{'spin',IS.sym }
%          just pick S.(SxS) term where location of spin q-label
%          is derived from IS.sym as returned by getLocalSpace.m
%
% Wb,Mar02,21

% outsourced from setup_HeisenbergSL_chiral.m

   if nargin<1
      helpthis, if nargin || nargout, wbdie('invalid usage'), end
      return
   end

   getopt('init',varargin);
     kflag=getopt('-k');
     tflag=getopt('-t');

     rchi =getopt('--chi-real');
     wchi =getopt('wchi',{});

   args=getopt('get_remaining');

   E0=QSpace(getIdentityQS(E0,args{:}));

   Z0=getIdentity(E0,'-0');
   dloc=getDimQS(E0); dloc=dloc(end,1);

   A2=getIdentity(E0,E0); % setitags(A2,{'s1','s2','K2'});

   EX=contract(getIdentity(E0,E0),'2',Z0,'1*',[1 3 2]);

   if 0
      Z=getIdentity(EX,3,'-0');
      EX_=contract(EX,'3*',Z,'1*','213');
   end

   if nargout>1 || tflag
    % build 2-site permutation operator
    % ----------------------------------------------------------------- %
    % NB! by construction, EX fuses the local state spaces of bra and ket
    % without truncation; contracted with EX' on the neighboring site,
    % then just crosslinks, i.e., exchanges them, which is equivalent
    % to the swap operator! // EX_IS_SWAP // liegroups.tex // Wb,Dec08,20
      X2=contract(A2,'12*',A2,'21');
      X2=contract(A2,3,contract(X2,2,A2,'3*'),1);

      [Sp,Sp_,Ix]=splitH4_SdotS(X2,'-v');
      e=norm(Sp-Sp_); if e>1E-12
         wbdie('unexpected Sp''.Sp structure (e=%g)',e); end

      if kflag || tflag, H2=QSpace(1,2);
         for k=1:2, if k==1, Q=Sp; else Q=EX; end
            Q=contract(contract(Q,'1*',A2,1),'23',Q,'32');
            H2(k)=contract(A2,'12*',Q,'13');
         end
         q=[norm(H2(1)-H2(2)), norm(H2(1))];

         e=q(1)/q(2); if e>1E-12, wbdie(['unexpected ' ...
             'operator setting (EX inequivalent from Sp @ %.3g)'],e);
         else wblog('ok.','got proper swap operator @ %.3g',e); 
         end
      end
   end

 % -------------------------------------------------------------------- %
   A3=contract(A2,3,getIdentity(A2,3,E0,1),1);
   x3=contract(A3,'231*',A3,'123'); % s1,s2,s3,s123* '+++-'

   X3=contract(A3,4,contract(x3,2,A3,'4*'),1); % rank-6 '+++---'

   if nargout>1 || tflag
      X3i=contract(A3,4,contract(x3',2,A3,'4*'),1);
      n3=(normQS(X3-X3i)/4)^2;
   end

   X3=contract(EX,'21',contract(X3,'36',EX,'12*'),'13',[2 3 1 4]);
   X3=skipzeros(X3);      % oss'o' -> ss'oo' with btlr index order

   if nargout<2, return; end

 % if X3 is clockwise, then X3' is counter clockwise
 % X3i=permute(X3,'2143*'); % same as X3! // norm(X3-X3i) < 1E-14
 %
 %               t                 %                b(!)
 %               |                 %                |
 %               ^                 %         Z*Z'   v     Z*Z'
 %         v     |      v          %          v     |      v
 % EX' ----<---- X3 ----<---- EX   %  EX' ---->---- X3'---->---- EX
 %         l     |      r          %          l     |      r
 %               ^                 %                v
 %               |                 %                |
 %               b                 %                t(!)
 %
 % X3i=contractQS(contractQS(X3,'3*',Z,'1*'),3,Z,'1',[2 1 3 4]);
 % WRN! not working either! // WRN_EX_DAGGER
 % => also need to contract X3i with EX from scratch to ensure consistency
 % => see X3i computed from x3' above

   X3i=contract(EX,'21',contract(X3i,'36',EX,'12*'),'13',[2 3 1 4]);
   X3i=skipzeros(X3i);

   if rchi, q=1/4; else q=1/4i; end
   chi=skipzeros(q*(X3-X3i));

   n2=norm(chi)^2;

   if n2<1E-6, wbdie(...
     'failed to obtain chiral operator (|chi|^2=%.3g',n2); end

   n=dloc;
   q=[ n2, (n/8)*(n^2-1) ];

   err=norm(diff(q)); if e>1E-12
        wblog('WRN','got unexpected chiral operator (e=%.3g)',err);
   else wblog('ok.','got proper chiral operator @ %.3g',err); 
   end

   Iout=add2struct('-',E0,EX,'X2=Sp',X3i,chi,dloc);
   Iout.real_chi=rchi;
   Iout.norm_chi2 = [n3 q];
   Iout.opts_chi=wchi;

   % if 1
   %    Chi=contract(A3,4,contract((x3-x3')/4,2,A3,'4*'),1);
   %    Chi=contract(EX,'21',contract(Chi,'36',EX,'12*'),'13',[2 3 1 4]);
   %    Chi=skipzeros(Chi);
   %    wbstop -> same as chi, indeed (@ 1E-15)
   % end

   if ~isempty(wchi), done=0;
      if iscell(wchi) && numel(wchi)==2 
         if ischar(wchi{1}) && ischar(wchi{1})
            if ~isempty(regexpi(wchi{1},'spin'))
               Iout.CHI=Iout.chi;
               Iout.opts_chi=wchi;
               Iout.chi=get_sub_SU2spin(chi,wchi{2}); done=1;
            end
         end
      end
      if ~done, wchi, wbdie('invalid usage'); end

      E1=getIdentity(Iout.chi,3);
      E2=getIdentity(Iout.chi,4);
      if ~isequal(E1,E2), wbdie('unexpected chiral operator'); end
      q={getDimQS(Iout.EX), getDimQS(E1)}; q=[ q{1}(end,3), q{2}(end,1) ];
      wblog(' * ','reducing dimension of EX based on chi (%d->%d)',q);

      Iout.Ex=contract(Iout.EX,3,E1,1);
   end
   Iout.err=err;

   if kflag, wbstop; end
end

% -------------------------------------------------------------------- %

function chi=get_sub_SU2spin(chi,sym)

   q=textscan(sym,'%s','whitespace',',; '); q=q{1};
   isym=findstrc(q,'SU2spin','-i');

   if numel(isym)~=1
      wbdie('unexpected symmetry (SU2spin required)'); end

   q=zeros(1,size(chi.Q{1},2)); q(isym)=2;
   Q=chi.Q; i=matchIndex([Q{3:4}],[q q]);

   if isempty(i)
        wbdie('failed to match symmetries for chiral operator');
   else chi=getsub(chi,i);
   end
end

% -------------------------------------------------------------------- %

