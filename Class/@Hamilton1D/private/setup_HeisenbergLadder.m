function [HAM]=setup_HeisenbergLadder(varargin)
% function [HAM]=setup_HeisenbergLadder([L, opts])
% Wb,Apr07,14 ; Wb,Aug28,15

% Wb,Sep25,15 :: added missing ~hconj for spin operators!
% (up to now, by adding hconj, thus simply resulted in twice the ground
% state energy; no other effect expected otherwise so far)

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     qloc = getopt('qloc', []);
     qend = getopt('qend' ,[]);
     qproj= getopt('qproj',[]);
     sym  = getopt('sym','SU2');
     B    = getopt('B', 0.);
     J    = getopt('J', [1 1]);
     L    = getopt('L', []);
     smoRB= getopt('smoRB',[]);
     perBC= getopt('-perBC');
     addS2= getopt('-addS2');
     tflag= getopt('-t');

  if isempty(L) L=getopt('get_last',[]); else
  getopt('check_error'); end

  naklt=0;

  if isempty(L)
     L=size(J,1);
     if L<=1
        wbdie('invalid usage (system length not specified)'); end
     wblog(' * ','got implicit length specification (L=%g)',L);
  elseif numel(J)==1,  J=repmat(J,L,2);
  elseif size(J,1)==1, J=repmat(J,L,1);
  end

  if ~perBC
     if     size(J,1)==L,   J(end,  1)=J(end-1,1);
     elseif size(J,1)==L-1, J(end+1,1)=J(end,  1); end
  end
  if L~=size(J,1) || size(J,2)~=2
     wbdie('invalid usage (L=%g; J: [%s]) !?',L,vec2str(size(J),'-f'));
  end

  if B, wbdie('finite B (%.3g) not yet implemented',B); end

  o={'-v'};
  if regexp(sym,'^SU\d+')
     N=str2num(sym(3:end)); if N>4, o={'-V'}; end
  end

  param=add2struct('-',L,sym,J,B,perBC,qloc,qend,qproj);
  q=uniquerows(J); if size(q,1)==1, param.J=J; end

  HAM=struct(Hamilton1D);

  q=uniquerows(J);
  if size(q,1)==1
     jstr=sprintf('J=[%g, %g]',q);
     HAM.info.mpo.J=q;
  else
     jstr=sprintf('J=[%.4g %.4g]..[%g %g]',min(J,[],1),max(J,[],1));
     HAM.info.mpo.J=J;
  end

  if ~isempty(smoRB)
     if numel(smoRB)<2, smoRB(2)=L/20; end
     if smoRB(1)<=1
        i0=L-3.6*smoRB(2);
     else
        i0=smoRB(1);
     end
     if i0<0.40*L
        wblog('WRN','got smoothing range of %.1f%% !?',100*i0/L); 
     end

     a=smoRB(2); f=0;
     for i=ceil(i0+1E-8):size(J,1)
        x=(i-i0)/a; f=1/(1+exp(((x^2-5)/x)));
        J(i,:)=J(i,:)*f;
     end

     s=sprintf('smoothening RB (i0=%g, a=%g) down to f=%.3g',i0,a,f);
     if f>1E-4, wblog('NB!','%s',s); else wblog('WRN','%s !?',s); end

     HAM.info.mpo.J=J;
  end

  if isequal(sym,'SU2')
     if isempty(qloc), qloc=[1];
     elseif numel(qloc)~=1 || qloc<1
        wbdie('invalid spin S=%g',qloc);
     end
     [S,IS]=getLocalSpace('Spin',qloc/2,o{:});
     if S.data{1}<0, S.data{1}=-S.data{1}; end

     if ~isfield(IS,'istr')
        IS.istr=sprintf('Sloc=%s',num2rat(qloc/2));
     end

     istr=sprintf(...
    'Heisenberg ladder (%s, %s) having %s',jstr,sym,IS.istr);
  else
     [S,IS]=getLocalSpace('Spin',sym,qloc,o{:});

     istr=sprintf(['Heisenberg ladder (%s)\n' ... 
    'having %s with qloc=(%s)'],jstr,IS.istr,sprintf('%g',IS.qloc));
  end

  IS_=IS; IS.S=S;
  for i=1:size(qend,1)
     if isequal(qend(i,:),qloc), IS(i+1)=IS(1);
     else 
        [s,q]=getLocalSpace('Spin',sym,qend(i,:),o{:});
        q.S=s; IS(i+1)=q;
     end
  end

  HAM=struct(Hamilton1D);
  HAM.info.istr=istr;
  HAM.info.param=param;
  HAM.info.IS=IS_;

  HAM.oez=init_ops(IS(1).E,'local identity operator (E)');
  HAM.ops=init_ops(IS(1).S,'local spin operator (S)','~hconj');

  if ~isempty(qend), s=size(qend);
     if s(1)>2 || s(2)~=length(qloc), qend
        wbdie('invalid usage (unexpected end-spin)');
     end
  end

  if ~isempty(qproj)
     A0=QSpace(getIdentityQS(S,1,S,1));
     [i1,i2,I]=matchIndex(A0.Q{3},qproj);

     if size(qproj,1)~=1 || isempty(i1), wbdie('invalid qproj'); end
     wblog('NB!','using rung-projection with q=(%s)',sprintf('%g',qproj)); 

     A1=getsub(A0,i1);

     S1=QSpace(contractQS(A1,'12*',contractQS(A1,1,S,2),'31'));
     S2=QSpace(contractQS(A1,'12*',contractQS(A1,2,S,2),'13'));

     HAM.oez(1)=init_ops(getIdentityQS(A1,3),'E_rung (projected)');
     HAM.info.S=S;

     e=[norm(S1-S2), norm(S1+S2)];

     if norm(e(1))<1E-12
        wblog(' * ','projection leads to same spin operator');
        Stot=S1+S2;
        HAM.ops(1)=init_ops(Stot,'S_rung (projected)','~hconj');

        if addS2, naklt=1;
        HAM.ops(2:3,1)=get_aklt_ops(Stot,naklt); end
     else
        wblog('WRN','projection leads to 2 distinct operators (%.3g,%.3g)',e);
        HAM.ops(1)=init_ops(S1,'Sa_rung (projected)','~hconj');
        HAM.ops(2)=init_ops(S2,'Sb_rung (projected)','~hconj');
     end

     HAM.info.istr=[istr sprintf(' and qproj=(%s)',sprintf('%g',qproj))];
     HAM.info.XY=[]; if ~addS2, J=J(:,1); end

     if ~isempty(qend)
        wblog('NB!','using end-spins with q=(%s)',sprintf('%g',qend')); 
        HAM.info.istr=[istr ' with end-spins'];

        m=size(qend,1);
        for i=1:m
           if m>1, if i==1, s='_L'; else s='_R'; end, else s=''; end
           HAM.oez(1,i+1)=init_ops(IS(i+1).E,['end-spin (E' s ')']);
           HAM.ops(:,i+1)=init_ops(IS(i+1).S,['end-spin (S' s ')'],'~hconj');
           if naklt
           HAM.ops(2:3,i+1)=get_aklt_ops(IS(i+1).S,naklt); end
        end

        if size(qend,1)==1
           L=L+1;
           J(end+1,:)=J(1,:);
        else
           L=L+2;
           J(end+(1:2),:)=J([1 1],:);
        end

        stype=ones(1,L); stype([1 end])=2;
        if size(qend,1)>1, stype(end)=3; end
     else stype=[]; end

     nops=size(HAM.ops,1)-naklt;
     HH=zeros(nops*(L-1),5); l=1;

     if perBC
        wblog('NB!','using periodic BC (interleaved setup)'); 
        L2=ceil(L/2);

        XY=[ (1:L)', zeros(L,1) ]; XY(2:2:end,2)=1;
        HAM.info.XY=XY;;

        for j=1:nops
           HH(l,:)=[ [1, j], [2, j], J(1,1)]; l=l+1;
        end
        for k=1:2:L-2
           for j=1:nops
              HH(l,:)=[ [k  , j], [k+2, j], J(k+1,1)]; l=l+1;
           end
           for j=1:nops
              HH(l,:)=[ [k+1, j], [k+3, j], J(k+2,1)]; l=l+1;
           end
        end

        if isempty(qend)
           if mod(L,2)==0, k=L-1;
              for j=1:nops
                 HH(l,:)=[ [k, j], [k+1, j], J(k,1)]; l=l+1;
              end
           else
              i=size(HH,1)-(0:(nops-1));
              HH(i,3)=L;
           end
        else
           if size(qend,1)~=1, wbdie('invalid qend'); end
           wblog('NB!','using end-spin (open right boundary!)'); 
           if mod(L,2)
              i=size(HH,1)-(0:(nops-1));
              HH(i,:)=[];
           end

           stype=ones(1,L); stype([end-1:end])=2;
        end
     else
        for k=1:L
           for j=1:nops
              HH(l,:)=[ [k, j], [k+1, j], J(k,1)]; l=l+1;
           end
        end
        l=size(HH,1)-(0:(nops-1));
        HH(l,:)=[];
     end

     if naklt
        if nops~=2, wbdie('unexpected nops=%g !?',nops); end
        i=find(HH(:,4)==2);
        HH(i,4)=3;
        HH(i,5)=J(HH(i,1),2);
     end

     HAM=setup_mpo(HAM,HH,stype);

     HAM.store='DMRG_HBL';

     return
  end

  if ~isempty(qend)
     wblog('NB!','using end-spins with q=(%s)',sprintf('%g',qend')); 
     HAM.info.istr=[istr ' with end-spins'];

     L=L+1;
     J(end+1,:)=J(1,:);
  end

% DMRG site order -------------------------------------------------- %
%
%   1 --- 3 --- 5 --- 7 ---  ...
%   |     |     |     |
%   |     |     |     |
%   2 --- 4 --- 6 --- 8 ---  ...
%
% ------------------------------------------------------------------ %

% XY=[
%   0 1    \   0 0 1 0
%   0 0    /  
%   1 1    \   1 1 1 0
%   1 0    /  
%   ...        ...
% ];

  XY=[ repmat((0:(L-1))',1,2), repmat([1 0],L,1) ];
  XY=reshape(permute(reshape(XY,[L,2,2]),[3 2 1]),2,[])';
  HAM.info.XY=XY;

  HH=zeros(2*(L-1)+L,5); l=1;

  if perBC
     wblog('NB!','using periodic BC (interleaved setup)'); 

     i=[3:4:2*L, 4:4:2*L];
     XY(i,2)=XY(i,2)+2;

     HH(l,:)=[ [1, 1], [3, 1], J(1,1)]; l=l+1;
     HH(l,:)=[ [2, 1], [4, 1], J(1,1)]; l=l+1;

     for k=1:L, i=2*k;
        HH(l,:)=[ [i-1, 1], [i  , 1], J(k,2)]; l=l+1;
        HH(l,:)=[ [i-1, 1], [i+3, 1], J(k,1)]; l=l+1;
        HH(l,:)=[ [i  , 1], [i+4, 1], J(k,1)]; l=l+1;
     end

     HH(end-1:end,:)=[];
     HH(end+(-2:-1),3)=HH(end+(-2:-1),3)-2;

     if ~isempty(qend)
        wbdie('invalid usage (got end-spin with perBC !?)');
     end
  else
     for k=1:L, i=2*k;
        HH(l,:)=[ [i-1, 1], [i  , 1], J(k,2)]; l=l+1;
        HH(l,:)=[ [i-1, 1], [i+1, 1], J(k,1)]; l=l+1;
        HH(l,:)=[ [i  , 1], [i+2, 1], J(k,1)]; l=l+1;
     end

     l=size(HH,1)+(-1:0);
     HH(l,:)=[];
     if ~isempty(qend)
        HH([3 end],:)=[];
        HH(:,[1 3])=HH(:,[1 3])-1;
        HH(1:2,1)=1;
        HH(end-1,3)=HH(end-2,3);

        XY([1 end],:)=[];
        XY([1 end],2)=0.5;
     end
  end

  HAM.info.XY=XY;
  HAM=setup_mpo(HAM,HH);

  HAM.store='DMRG_HBL';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

