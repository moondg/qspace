function [A1,A2,r2,I]=ortho2site(A1,A2,kdir,varargin)
% function [A1,A2,r2,I]=ortho2site(A1,A2,kdir[,rtol[,Nkeep]][,opts]);
%
%    Orthonormalize the two consecutive A-tensor
%    (..-A1-A2-..) as specified by kdir \in {'>>','<<'}.
%    Truncation is performed via a threshold on the
%    eigenvalues of the explicitly computed reduced
%    density matrix R (default: rtol=1E-14).
%
% Options
%  
%   '-l' lenient (perform orthonormalization even if non-standard)
%   '-b' switch to bond picture (e.g. for bond update)
%
% Further output
%
%    r2 returns the discarded weight of the reduced
%    density matrix R. I contains further information
%    on eig(R).
%
% NB! this is a private function (only for member functions)
% Wb,Apr08,14 ; Wb,Jul06,16

% Wb,Jul06,16: allow global state index either on A1 or A2
% for former version, see Archive/ortho2site_160728.m

  if nargin<3 || numel(A1)~=1 || numel(A2)~=1
     helpthis, if nargin || nargout, wberr(['invalid usage ' ... 
         '(input A1 and A2 must contain to consecutive A-tensors)']);
     end, return
  end

  rtol=[];
  if ~isempty(varargin) && isnumber(varargin{1})
     rtol=varargin{1}; varargin(1)=[]; end
  if isempty(rtol) || rtol<0, rtol=1E-14;
  elseif rtol>0.1, wberr('got likely invalid rtol=%g / 0.1',rtol); end
  os={'stol',sqrt(rtol)};

  Nkeep=[];
  if ~isempty(varargin) && isnumber(varargin{1})
     Nkeep=varargin{1}; varargin(1)=[]; end
  if ~isempty(Nkeep)
     if ~Nkeep, wberr('invalid Nkeep=%g',Nkeep); end
     os={os{:},'Nkeep',Nkeep};
  end

  getopt('init',varargin);
     lflag=getopt('-l');
     bflag=getopt('-b');
  getopt('check_error');

  if ~isa(A1,'QSpace'), A1=QSpace(A1); end
  if ~isa(A2,'QSpace'), A2=QSpace(A2); end

  kdir=check_dir(kdir); i=1:3;

  if lflag
       s='invalid orthonormalziation';
  else s='not current site'; end

  if kdir>0
     if lflag, i=[1 3]; end
     q=getqdir(A1); if numel(q)<3 || any(q(i)<1), A1, wberr('A1 %s',s); end
     t=regexprep(getitags(A1,2),'^[A-Z]*','E');
     E=getIdentity(A1,1,A1,3, t);
     X=contract(E,'*',A1); % R'R;g // A_LRsg_ORDER
  else
     if lflag, i=[2 3]; end
     q=getqdir(A2); if numel(q)<3 || any(q(i)<1), A2, wberr('A2 %s',s); end
     t=regexprep(getitags(A2,1),'^[A-Z]*','E');
     E=getIdentity(A2,2,A2,3, t);
     X=contract(E,'*',A2); % L'L;g // A_LRsg_ORDER
  end

  if bflag
     [B,S,U,I]=svdQS(X,1,os{:});
  else
     ot={'^[A-Z]+','X'};
     t=getitags(X,2); itagrep(X,2,ot{:}); 
     if kdir>0
          itagrep(A2,1,ot{:});
     else itagrep(A1,2,ot{:}); end

     [U,B,I]=orthoQS(X,1,os{:},'itag',t);
  end

  I.DX=getDimQS(X);
  Dk=getDimQS(U); I.Dk=Dk(:,2);

  rb=numel(B.Q);
  if rb<2 || rb>3, wberr('got rank(B)=%g !?',rb); end

  if kdir>0
     p2=[]; if rb==3
        p2=initperm(numel(A2.Q)+1,'--2end',2);
     end

     A1=contract(E,U,[1 3 2]);
     A2=contract(B,A2,p2);
  else
     p1=[]; if rb==3, r1=numel(A1.Q);
        p1=initperm(r1+1,'--2pos',{r1,2});
     elseif rb==2, p1=[1 3 2];
     else wberr('got rank(B)=%g !?',rb); end

     A1=contract(A1,B,p1);
     A2=contract(U,E);
  end

  ot={'^[A-Z]+','K'};
  itagrep(A1,2,ot{:}); itagrep(A2,1,ot{:});

  if bflag
     if kdir>0, S=permuteQS(S,'21'); end
     S=itagrep(QSpace(S),ot{:});

     ot={'^[A-Z]','X'}; % ot={'(\d+)','$1~'};
     if     kdir>0, S=itagrep(S,2,ot{:}); A2=itagrep(A2,1,ot{:});
     elseif kdir<0, S=itagrep(S,1,ot{:}); A1=itagrep(A1,2,ot{:}); end

     I.S=diag(S);
  end

  if nargout>3, I=add2struct(I,'dir=kdir'); end
  r2=I.svd2tr;

% AB2=contractQS(A1,2,A2,1); e=norm(QSpace(AB1)-AB2)^2;
% if e>1E2*rtol
%    wblog(iff(e<1E-10,'WRN','ERR'),...
%       'invalid ortho2site (%.3g/%g @ r2=%.3g)',e,rtol,r2);
%    if e>1E-6, keyboard, end
% elseif e>1E-20, wblog('TST','e=%.3g',e);
% end

end

