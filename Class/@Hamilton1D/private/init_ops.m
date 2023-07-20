function S=init_ops(x_,varargin)
% function S=init_ops(x [,'-ferm'])
%
%   NB! fermionic and hconj needs to be specified in the caller
%   routine; the remaining field qop is for convenience as it
%   can be derived from ops.
%
% Options
%
%  '-ferm'  whether got fermionic operators in that Jordan Wigner
%           string is required when computing matrix elements.
%
% Wb,Aug25,15

  nx=numel(x_);
  if nx<1, wbdie('got empty operator'); end

  getopt('init',varargin);
     qflag=getopt('-q');
     f=getopt('-ferm');

     if     getopt('-hconj'), hconj=1;
     elseif getopt('~hconj'), hconj=0; elseif f, hconj=1;
     else
        hconj=0;
        for i=1:nx, xi=x_(i);
           if ~isequal(xi.Q{1},xi.Q{2})
              hconj=2; break;
           elseif numel(xi.Q)>2 && any(reshape([xi.Q{3:end}],1,[]))
              hconj=3; break;
           elseif numel(xi.data)==1
              H=xi.data{1}; s=size(H);
              if diff(s), hconj=4; else q=norm(H,'fro');
                 if q && norm(H-H','fro')/q>1E-14;, hconj=5; end
              end
           end
        end
     end

  istr=getopt('get_last','');
  if ~ischar(istr), wbdie('invalid usage (istr must be string)'); end

  S=repmat(struct(...
     'info',istr, ...
     'op', [], ...
     'dop',[], ...
     'qop',[], ...
     'fermionic',f, ...
     'hconj',hconj ...
    ),nx,1);

  % ------------------------------------------------------------ %
  % *) NB! qop is unique for an IROP; however, in order to deal  %
  % with  former 'hcpow', allow composite operators that         %
  % combine several multiplet qop's!                             %
  % e.g. consider (S.S')^2 for spin-1/2 AKLT, which results in   %
  %                                                              %
  %       |                       |           |         |        %
  %       ^                       ^           ^         ^        %
  %       |    s=1                |           |   s=1   |        %
  %      [S]->-               ->-[S]         [S]--->---[S]       %
  %       |    \  s={0,1,2}  /    |           |         |        %
  %       |    [X]---->----[X]    |     =     |         |        %
  %       |    /             \    |           |   s=1   |        %
  %      [S]->-               ->-[S]         [S]--->---[S]       %
  %       |    s=1                |           |         |        %
  %       ^                       ^           ^         ^        %
  %       |                       |           |         |        %
  %                                                              %
  % which represents a *mix* of scalar with a non-scalar IROPs!  %
  % tags: QOP_ALLOW_MULT // Wb,Aug26,15                          %
  % ------------------------------------------------------------ %

  for i=1:nx, xi=x_(i);
     if numel(xi.Q)==3 && isempty(xi.info.otype)
     xi.info.otype='operator'; end

     S(i).op=QSpace(xi); d=getDimQS(xi); d=d(end,:); d(end+1:3)=1;
     S(i).dop=d(3);

     if nx>1
        if regexp(istr,'^[a-zA-Z]+$')
             S(i).info=sprintf('%s%g',   istr,i);
        else S(i).info=sprintf('%s [%g]',istr,i); end
     end

     Q=xi.Q; r=numel(Q);
     if r==2
        if isequal(Q{1},Q{2})
             q=zeros(1,size(Q{1},2));
        else q=uniquerows(Q{1}-Q{2});
        end
     elseif r==3, q=uniquerows(Q{3});
        if ~isequal(getqdir(xi),[1 -1 -1])
           disp(xi); wbdie('unexpected rank-%g operator',r);
        end
     elseif r==4
        if ~isequal(getqdir(xi),[1 -1 -1  1]) ...
        || ~isequal(getQDimQS(xi,3),getQDimQS(xi,4)), disp(xi)
           wbdie('unexpected rank-%g operator',r); end
        z=contractQS(xi,'1234',xi,'2143');
        q=uniquerows(cat(1,Q{3:end}));
     else disp(xi)
        wbdie('invalid rank-%g operator',r);
     end

     if size(q,1)~=1
        if qflag, s={'-->',''}; else s={'NB!','got '}; end
        wblog(s{1},[s{2} 'composite `%s'''],S(i).info);
     end
     S(i).qop=q;

     if hconj<0
      % auto determine hconj // QOP_ALLOW_MULT ok!
      % Note that e.g. the Heiseberg interaction Si.Sj'
      % in between sites i and j (i!=j) is already
      % perfectly fine *without* h.conj.(!) even though
      % S does not transform as a scalar operator!
      % => such cases need to be specified externally.
        if norm(S(i).qop), S(i).hconj=1; end
     elseif hconj && norm(S(i).qop)==0
        if numel(xi.Q)>2, xi=fixScalarOp(xi); end
        e=normQS(xi-xi');
        if S(i).fermionic
           wbdie('invalid usage (local operator cannot be ferminoic)');
        elseif e>1E-6
             wblog('NB!','using hconj=%g for local operator',  hconj);
        else wblog('WRN','got hconj=%g for scalar operator !?',hconj);
        end
     end
  end

end

