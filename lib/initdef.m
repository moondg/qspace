function val=initdef(vname,val)
% Function: val=initdef(vname,val)
%
%    If variable <vname> is not set yet in workspace
%    it is set to default <val>
%
% DEPRECATED - see setdef(). Wb,Mar01,08
% Wb,Jun15,07

  global initdef_val__ initdef_flag__

  if nargin~=2
     eval(['help ' mfilename]);
     if nargin, wbdie('invalid usage'); else return; end
  end

  if nargout
     initdef_flag__=0;
     evalin('caller', sprintf(['global initdef_val__ initdef_flag__; ' ...
     'if exist(''%s'',''var''), initdef_val__=%s; initdef_flag__=1; end'],...
     vname,vname));
     if initdef_flag__, val=initdef_val__; end
  else
     initdef_val__=val;
     evalin('caller', sprintf(...
     'global initdef_val__; if ~exist(''%s'',''var''), %s=initdef_val__; end',...
     vname,vname)); clear val
  end

  clear global initdef_val__ initdef_flag__

end

