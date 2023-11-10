function [varargout]=eig(H,varargin)
% function [varargout]=eig(H,varargin)
%    ee=eig(H,varargin);
%    [A,E]=eig(H,varargin);
% Wb,Aug30,12

  getopt('init',varargin);
     qsflag=getopt('-QS'); if ~qsflag
     xflag =getopt('-x' ); end
  varargin=getopt('get_remaining');

  if ~isscalarop(H)
     wbdie('invalid usage (scalar operator required)'); end

  if qsflag && nargout<2
     for i=1:numel(H.data)
        if 1
           q=eig(H.data{i});
           if ~isreal(q) || any(diff(q)<0)
              e=norm(imag(q))/norm(q);
              if e>1E-10, wblog('WRN','got non-symmetric H !? (%.3g)',e);
              else q=sort(real(q)); end
           end
           H.data{i}=q;
        else
           q=H.data{i};
             e=norm(q-q','fro')/norm(q,'fro');
             if e>1E-10, wblog('WRN','got non-symmetric H (%.3g)',e);
             else q=0.5*(q+q'); end
           H.data{i}=eig(q);
        end
     end
     varargout{1}=H;
     return
  end

  [ee,I]=eigQS(H,varargin{:});

  if nargout<=1
     if qsflag, varargout={ QSpace(I.EK) };
     elseif xflag
        if isfield(I.EK.info,'cgs') && ~isempty(I.EK.info.cgs)
           EK=I.EK; dd=getzdim(QSpace(EK),1,'-p');
           for k=1:numel(EK.data)
              EK.data{k}=reshape(repmat(EK.data{k},dd(k),1),1,[]);
           end
           varargout={sort(cat(2,EK.data{:}))};
        else varargout={ ee }; end
     else varargout={ ee }; end
  elseif nargout==2, varargout={ QSpace(I.AK), QSpace(I.EK) };
  else
     I.ee=ee;
     varargout={QSpace(I.AK), QSpace(I.EK), I};
  end

end

