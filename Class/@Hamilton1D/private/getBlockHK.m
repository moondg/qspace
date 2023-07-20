function HK=getBlockHK(Ik,varargin)
% function HK=getBlockHK(Ik [,use_hconj])
%
%    Extract block Hamiltonian from Ik with is a structure
%    with two expected fields: HK and AK. For pseudo-MPO,
%    this just returns Ik.HK, while also checking H.c.
%
%    When using full MPO for Hamiltonian, however, HK needs
%    to be tweaked to obtain proper part of HK that represents
%    *completed* block Hamiltonian.
%
%    Default use_hconj uses Ik.info.hconj, or auto if not set.
%
% See also setupStorage(), use_Hconj().
% Wb,Dec04,20

% merged with former initNRG/fix_HK.m // Wb,Dec06,20

   getopt('init',varargin);
      zfix=getopt('--fix'); % fix by adding `zero' blocks
   hconj=getopt('get_last',-1);

   HK=Ik.HK; r=numel(HK.Q);
   if ~r, if ~zfix, return; end
   elseif r==3
      HK=getsub(HK,'-0');
      k=get_kidx(HK);

      i=itags2odir(Ik.AK); n=size(HK.data{1},3); 
      if i==1, ic=1; elseif i==2, ic='1*';
      else wbdie('invalid setting (odir=%g)',i); end

      v=zeros(n,1); if i<=n, v(i)=1; end
	  V=getvac(HK,'-1d'); V.data{1}=v;

      HK=contract(HK,3,V,ic); % `squeeze'

      if hconj<0 && isfield(Ik.info,'hconj') && Ik.info.hconj
         hconj=10;
      end

   elseif r~=2
      wbdie('invalid usage (got rank-%g tensor)',r);
   end

   q=normQS(HK);
   i=itags2odir(Ik.AK); if ~i, i=1; end

   if HK
        k=get_kidx(   HK.info.itags(1));
   else k=get_kidx(Ik.AK.info.itags(i)); end
   k=abs(k);

   if hconj>0, HK=HK+HK';
   elseif q
      q(2)=normQS(HK-HK'); e=q(2)/q(1);
      if e>1E-6 && hconj, HK=HK+HK';
      elseif e>1E-12
         wblog('WRN','got non-Hermitian H(%g) @ %.3g',k,e);
      end
   else q=1; end

   if ~zfix, return; end

   da=getDimQS(Ik.AK); da=da(1,:);
   dh=getDimQS(HK);

   if ~isempty(dh), dh=dh(1,1:2);
   else
      wblog('WRN','got empty H(%g)',k); dh=0;
   end

   if any(dh~=da(i))
      wblog(' * ','adding zero-blocks to H(%g): D=%g->%g',k,dh(1),da(i));
      E=getIdentity(Ik.AK,i);
      HK=HK+(1E-24*q(1))*E;
   end

end

% -------------------------------------------------------------------- %

