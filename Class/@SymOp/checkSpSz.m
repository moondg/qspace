function cr=checkSpSz(Sp,Sz,sym)
% function cr=checkSpSz(Sp,Sz,sym)
% Wb,Jun28,10

  m=numel(Sz); n=numel(Sp); cr=zeros(m,n);

  for i=1:m
  for j=1:n
      q=comm(Sz(i),Sp(j)); if norm(q)<1E-12, continue; end
      cr(i,j)=olap(q,Sp(j)) / olap(Sp(j),Sp(j));

      e=norm(q-cr(i,j)*Sp(j)); if e>1E-12
      wberr('invalid Sp/Sz pair (%d/%d)',i,j); end
  end
  end

  if nargin<3, return; end

  N=-1; sym_=sym;
  if regexp(sym,'^[SZ]');
     i=regexp(sym,'\d+$'); if ~isempty(i)
     N=str2num(sym(i:end)); sym=[sym(1:i-1) 'N']; end
  end

  switch sym

    case {'A','P','ZN'}
       if ~isempty(Sp) || numel(Sz)~=1 || ~isempty(cr)
       wberr('invalid Ablian Sz/Sp setting'); end

    case 'SUN'
       if isempty(N) || N<2 || N>12
          wberr('invalid symmetry (%s)',sym_);
       elseif N>2
          check_SU(N,Sp,Sz,cr,sym_);
       else
          if numel(Sz)~=1 || numel(Sp)~=1 || ...
           ~isnumber(cr) || abs(cr-1)>1E-12
          wberr('invalid SU(2) Sz/Sp pair (%.3g)',abs(cr-1)); end
       end

    case 'SpN'

       if isempty(N) || mod(N,2) || N>12
            wberr('invalid symmetry (%s)',sym);
       else check_Sp(N,Sp,Sz,cr,sym_);
       end

    otherwise wberr('invalid symmetry (%s)',sym_);
  end
end

% -------------------------------------------------------------------- %
% outsourced from above // Wb,May07,21

function check_SU(N,Sp,Sz,cr,sym)

   r=N-1; e=-1;
   if r<2 || numel(Sz)~=r || numel(Sp)~=r
      wberr('invalid %s operators',sym); end

   if r==2
      e=norm(cr - [1 -0.5; 0 1.5],'fro');
   elseif r<=9
      e=(0:.5:r); q=diag(e(3:2+r)) + diag(-e(2:r),1);
      e=norm(cr - q,'fro');
   end
   if abs(e)>1E-12, wberr('invalid %s operators (e=%g)',sym,e); end

   for i=1:r
   for j=2:r
      Sij=comm(Sp(i),Sp(j)); e=trace(Sij);
      if abs(e)>1E-12
         wberr('got trace-full %s generator! (%g,%g; %.3g)',sym,i,j,e);
      end
   end
   end
end

% -------------------------------------------------------------------- %
% Wb,Nov29,11 % Sp4: Wb,Apr13,12
% outsourced from above // Wb,May07,21

function check_Sp(N,Sp,Sz,cr,sym)

   r=N/2;
   if numel(Sz)~=r || numel(Sp)~=r
      wberr('invalid %s operators',sym); end

   if r==2
      e=norm(cr - [2 -2; 0 2],'fro');
   elseif r==3
      e=norm(cr - [2 -1 0; 0 3 -4; 0 0 2],'fro');
   elseif r==4
      e=norm(cr - [2 -1 0 0; 0 3 -2 0; 0 0 4 -6; 0 0 0 2],'fro');
   else e=-1; end
   if abs(e)>1E-12, wberr('invalid %s operators (e=%g)',sym,e); end

   S3=comm(Sp(1),Sp(2)); e=trace(S3);
   if abs(e)>1E-12
      wberr('got trace-full 3rd %s generator! (%.3g)',sym,e);
   end
end

% -------------------------------------------------------------------- %

