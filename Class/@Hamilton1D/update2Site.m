function varargout=update2Site(HAM,k1,k2,dir,varargin)
% function update2Site(HAM,...)
% Usage #1: [r2,I]=update2Site(HAM,k1,k2,dir [,rtol,Nkeep])
%    this loads and stores data
% Usage #2: [X1,X2,I]=update2Site(HAM,X1,X2,dir [,rtol,Nkeep])
%    this just transformes given nearest neighbor pair of sites
% Wb,Apr20,14

  if nargin<4
     helpthis, if nargin || nargout, wberr('invalid usage'), end
     return
  end

  getopt('init',varargin);
     iflag=getopt('-i');
  args=getopt('get_remaining');

  rtol=[]; Nkeep=[];
  for i=1:numel(args), q=args{i};
     if isnumber(q) && q>0
        if q<1 && isempty(rtol), rtol=q; continue;
        elseif q>1 && q==round(q) && isempty(Nkeep), Nkeep=q; continue;
        end 
     end
     wberr('invalid usage (i=%g: %g)',i,q);
  end

  i=[ isstruct(k1), isstruct(k2) ];
  if xor(i(1),i(2)), wberr('invalid usage'); end
  gotX=i(1);

  if gotX
     X1=k1; k1=str2num(regexprep(X1.AK.info.itags{3},'^\w(\d+)$','$1'));
     X2=k2; k2=str2num(regexprep(X2.AK.info.itags{3},'^\w(\d+)$','$1'));
  end

  if ischar(dir), dir=check_dir(dir); end
  if iflag, if k1>=k2, wberr(...
    'invalid usage (to unsorted k''s [%g %g]',k1,k2); end
  elseif k1+1~=k2, wberr(...
    'invalid usage (expecting two consecutive sites; got [%g %g]',k1,k2);
  end

  if ~gotX
     X1=load_dmrg_data(HAM,k1);
     X2=load_dmrg_data(HAM,k2);
  end

  oo={};
  if isempty(rtol)
     rtol=getopt(HAM.info,'sweep','rtol',[]); end
  if ~isempty(rtol) && rtol>0, oo{1}=rtol; end

  if isempty(Nkeep)
     Nkeep=getopt(HAM.info,'sweep','Nkeep',[]);
     d1=getDimQS(X1.AK); 
     d2=getDimQS(X2.AK); q=min([ d1(1,1:2), d2(1,1:2) ]);
     if Nkeep<q, Nkeep=q; end
  end
  if ~isempty(Nkeep), oo{2}=Nkeep; end

  [X1.AK,X2.AK,r2,I]=ortho2site(X1.AK,X2.AK,dir,oo{:});

  if dir>0
     X1=updateHK(HAM,k1,'>>',X1);
     X2.info.odir=0;
     X2.Psi=X2.AK;
  else
     X2=updateHK(HAM,k2,'<<',X2);
     X1.info.odir=0;
     X1.Psi=X1.AK;
  end

  X1.info.itags=X1.AK.info.itags;
  X2.info.itags=X2.AK.info.itags;

  if ~gotX
     save_dmrg_data(HAM,X1,k1);
     save_dmrg_data(HAM,X2,k2);
     if nargout, varargout={r2,I}; end
  else
     varargout={X1,X2,I};
  end

end

