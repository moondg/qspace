function [HPsi,E]=get_HPsi(Psi,X1,X2,varargin)
% function [HPsi,E]=get_HPsi(Psi,X1,X2 [,J,hconj])
% compute H|psi> within the bond update setting
% Wb,Jan20,16

% outsourced from sweepPsi.m // Wb,Jan20,16

  persistent init_DMRG_mpo_fix

% allow global state index 'g' with Psi // Wb,Jul28,16
% index order convention: EE'[g]
  qdir=getqdir(Psi,'-s'); r=rank(Psi);

  if     isequal(qdir,'++' ), NPsi=0; p132=[];
  elseif isequal(qdir,'++-'), NPsi=1; p132=[1 3 2];
  else wberr('unexpected rank-%g Psi (%s)',r,qdir); end

  k=get_kidx(Psi);
  if diff(k(1:2)) || (r==3 && ...
    ~isnan(k(3)) && isempty(regexpi(Psi.info.itags{3},'^(dav|Psi)')))
     Psi, wberr('unexpected Psi for bond-update');
  end

  if nargin==4 && ischar(varargin{1}) && ...
     regexp(varargin{1},'^--full-mpo')

     q={ X1.info.hconj, X2.info.hconj };
        if isequal(q{:}), q=q{1};
        else disp(q), wberr('inconsistent hconj'); end
        if numel(q)~=1 || q<0, disp(q), wberr('invalid hconj'); end
     hconj=q;

     k1=get_kidx(X1.HK);
     k2=get_kidx(X2.HK); X2_HK=X2.HK;

     if ~isequal(k1(1:2),k2(1:2)) || sum(k1(1:2)) || sum(k2(1:2))
        wberr('invalid MPO setting');
     elseif k1(3)+k2(3)
        if regexp(varargin{1},'INIT')
           if isempty(init_DMRG_mpo_fix), init_DMRG_mpo_fix=1; 
              wblog('WRN',['fusing mpo-bonds k=%g~%g ' ...
               '(assuming initDMRG)'],k1(3),-k2(3));
           else init_DMRG_mpo_fix=init_DMRG_mpo_fix+1;
           end
           setitags(X2_HK,X1.HK.info.itags{3},3);
        else
            wbdie('got Hamilton MPO inconsistency (%g/%g)',k1(3),k2(3));
        end
     end

     HPsi=QSpace(contractQS({X1.HK,Psi},X2_HK,p132));
     if hconj
        HPsi=HPsi+contractQS({X1.HK,'*',Psi},X2_HK,'*',p132);
     end

  elseif numel(varargin)==2
     J=varargin{1}; hconj=varargin{2};

     HPsi=QSpace(contractQS(X1.HK,Psi)) + contractQS(Psi,X2.HK,p132);

     J=reshape(J,1,[]);
     for j=J

        if hconj(j), m=2; else m=1; end
        for ic=1:m
           if ic==1
              o1={'*'}; o2={};
           else 
              o1={}; o2={'*'};
           end

           try
              Q=contractQS({X1.OP(j),o1{:},Psi},X2.OP(j),o2{:},p132);
              HPsi=HPsi+Q;
           catch me
              save2tmp
              save2(sprintf('%s_%s_dbg.mat',wbstamp,mfilename));
              rethrow(me)
           end
        end
     end
  else disp(varargin), wberr('invalid usage'); end

  if nargout>1
     E=getscalar(QSpace(contractQS(Psi,'*',HPsi)));
  end

end

