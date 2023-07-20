function [HAM]=setup_AKLT(varargin)
% function [HAM]=setup_AKLT(varargin)

% outsourced from Hamilton1D()
% deprecated / old version // Wb,Aug28,15

  if nargin<1
     helpthis, if nargin || nargout, wbdie('invalid usage'), end
     return
  end

  getopt('init',varargin);
     qloc = getopt('qloc',[]);
     sym  = getopt('sym','SU2');
     B    = getopt('B', 0.);
     J    = getopt('J', 1.);
     L    = getopt('L', []);
     perBC= getopt('-perBC');
     tflag= getopt('-t');
  if isempty(L)
       L=getopt('get_last',[]);
  else getopt('check_error'); end

  if isempty(L), wbdie('length L not specified'); end
  if B, wbdie('finite B (%.3g) not yet implemented',B); end

  o={'-v'};
  if regexp(sym,'^SU\d+')
     N=str2num(sym(3:end)); if N>4, o={'-V'}; end
  end

  HAM=struct(Hamilton1D);

  if isequal(sym,'SU2')
     if isempty(qloc), qloc=[1];
     elseif numel(qloc)~=1 || qloc<1, wbdie('invalid spin S=%g',qloc);
     end
     [S,IS]=getLocalSpace('Spin',qloc/2,o{:});
     if S.data{1}<0, S.data{1}=-S.data{1}; end

     [~,q]=getLocalSpace('FermionS','SU2spin','NC',2); q=q.E;
     for i=1:numel(q.data), q.data{i}=1; end
     IS.E=q;

     if ~isfield(IS,'istr')
        IS.istr=sprintf('Sloc=%s',num2rat(qloc/2));
     end

     HAM.info.istr=sprintf(...
       'Heisenberg Hamiltonian (J=%g, %s) having %s',...
        J,sym,IS.istr);

  else wbdie('invalid usage'); end

  qz=zeros(1,size(S.Q{3},2));

  E=contractQS(S,[1 3],S,[1 3],'conjA');
  for i=1:numel(E.data), E.data{i}=eye(size(E.data{i})); end
  E=QSpace(E);

  HAM.info.IS=IS;

  HAM.ops=[ E, sqrt(abs(J))*S ]; J=sign(J);
  HAM.qop=[ qz; S.Q{3} ];

  HAM.hconj=[0 0];

  HAM.info=[];

  HAM.dop=zeros(size(HAM.ops));
  for i=1:numel(HAM.ops)
     s=dim(HAM.ops(i),3,'-a'); HAM.dop(i)=s(end);
  end

  if ~perBC
     M=cat(3, ...
         [
           1 0 0
           0 1 0
           0 0 0
         ], ...  
         [
           0 0 1
           0 0 0
           0 J 0
         ] ...
     );

     HAM.mpo=repmat(struct('M',M,'iop',[0 0 1],'user',[]),1,L);
     HAM.mpo(1).M=M(1,:,:);
     HAM.mpo(L).M=M(:,2,:);

  else
     M=cat(3, ...
         [
           1 0 0 0
           0 1 0 0
           0 0 0 0
           0 0 0 1
         ], ...
         [
           0 0 1 0
           0 0 0 0
           0 J 0 0
           0 0 0 0
         ] ...
     );

     HAM.mpo=repmat(struct('M',M,'iop',[0 0 2 2],'user',[]),1,L);
     q=M; q(1,4,2)=1; HAM.mpo(1).M=q(1,:,:);
     q=M; q(4,2,2)=J; HAM.mpo(L).M=q(:,2,:);
  end

  HAM.store='DMRG';

  if tflag, keyboard, end

end

% -------------------------------------------------------------------- %

