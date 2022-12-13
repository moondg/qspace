function  h=addt2fig(varargin)
% Function h=addt2fig([tstamp, ] sflag)
% add time stamp at the right top corner of gcf
% default time is now(), but can be altered by numeric value tstamp.
%
%     default: Dec11,07 11:43:58
%     'M'      Muenchen, Mon dd, yyyy
%     'H'      hostname, Mon dd, yyyy
%     'wb'     AW, mm/dd/yy
%     'Wb'     AW, Mon dd, yyyy
%     'MWb'    AW, Muenchen, Mon dd, yyyy
%
% Wb,Dec01,00

  if nargin && isnumeric(varargin{1})
  t=varargin(1); varargin(1)=[]; else t={}; end
  narg=length(varargin);

  if narg>0, sflag=varargin{1}; tstr=time(t{:},'-l');
     switch(sflag)
        case 'wb',  tstr=['AW ' datestr(now,'mm/dd/yy')];
        case 'Wb',  tstr=['AW, ', tstr(1:12)];
        case 'MWb', tstr=['AW, Munich, ', tstr(1:12)];
        case 'WB',  tstr=['AW, ', tstr      ];
        case 'M',   tstr=['Munich, ',     tstr(1:12)];
        case 'H',   tstr=[getenv('HOST'), ', ', tstr(1:12)];
     otherwise
        wblog('ERR','\NInvalid flag %s\N', sflag);
        error(' ');
     end
  else
     tstr=time(t{:});
  end

  h = header('hright', tstr);

  s=get(0,'DefaultAxesFontSize');
  set (h(1),'FontSize',0.7*s);

  if nargout==0, clear h, end

end

