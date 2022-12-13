function [I,xy] = intxy(xy,yd, varargin)
% FUNCTION : intxy - integrate xy data
% USAGE 1  : I = intxy(xy_data [,meth])
% USAGE 2  : I = intxy(xdata, ydata [,meth])
%
%    meth = {'linear'}, 'spline', 'cubic'
%    returns rough estimate of integral of discrete [x,y] data
%
% Wb,Mar22,05  May25,06

  I=0;

% put all data into [ X, Y ] data block
% NB! also accept complex integrals => no conjugation, just transpose!
  if nargin==1 || nargin>1 && ischar(yd)
     [m n]=size(xy); if m<n, xy=xy.'; end
     [m n]=size(xy); if n<2, wblog('ERR Invalid XY data'); return; end
  elseif nargin>1
     if ~isvector(xy), wblog('ERR X data is not a vector'); return; end
     m=length(xy);
     if     m==size(yd,1), xy=[xy(:), yd  ];
     elseif m==size(yd,2), xy=[xy(:), yd.'];
     else wblog('ERR XY data does not match'); return; end
  else
     eval(['help ' mfilename]); return
  end

  if nargin==2 && ischar(yd), varargin={yd}; end

  dxy=diff(xy(:,1));
  if any(dxy<0) && any(dxy>0)
  wblog('WRN input data not sorted!! (%g)', min(diff(xy(:,1))) ); end

  if ~isempty(varargin), meth=varargin{1}; else meth='linear'; end

  switch lower(meth)

   case 'linear'

     for m=2:size(xy,2)
     I(m-1) = trapz(xy(:,1), xy(:,m)); end

   case {'spline','cubic'}

     xx=linspace(min(xy(:,1)), max(xy(:,1)), min([10*size(xy,1), 1E4]));
     dx=xx(2)-xx(1);

     for m=2:size(xy,2)
        pp=spline(xy(:,1), xy(:,m));

        yy=ppval(pp, xx); YY{m-1}=yy(:);
        y2=interp1(xy(:,1), xy(:,m), xx, 'cubic');

        q=norm(yy-y2)/sqrt(length(yy));
        if q>0.1, wblog('WRN Interpolation quality = %g', q); end

        if strcmp(lower(meth),'spline')
             I(m-1)=dx*sum(yy);
        else I(m-1)=dx*sum(y2); end
     end

   otherwise, wblog('ERR - invalid method >%s<', meth);
  end

  if nargout>1
  xy = [xx(:), cat(2,YY{:})]; end

end

