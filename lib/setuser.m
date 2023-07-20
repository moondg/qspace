function u=setuser(ah,varargin)
% function u=setuser(h,'name1',value1, ...)
% function u=setuser(h,var1,var2,...)
%
%    Set or update field <name1> of UserData of given handle(s)
%    to specified value. If the field name is not explicitly specified,
%    the name of the variable in the caller is used (which thus must
%    exist). A mix of the above thow usages is permitted.
%
%    NB! If the UserData of a specified handle is already initialized
%    to a non-structure data, this will be erased without any warning.
%
% Examples
%
%    setuser(gca,a,b,c,'some_name',pi/3)
%
% This serves to assign global variables that survives `clear all'.
% Examples for h:
%
%    0      system root handle
%    groot  graphics root handle
%    gcf    current figure
%    gca    current axis
%
% Wb,Sep02,08 ; Wb,Apr07,22

  if nargin<2 || ~all(ishandle(ah))
     eval(['help ' mfilename]);
     if nargin || nargout, wbdie('invalid usage'), end, return
  end

  i=1; nargs=numel(varargin); vx=cell(floor(nargs/2),2); l=1;
  while i<=nargs,   x=varargin{i};
     if ~ischar(x), v=inputname(i+1);
        if isempty(v), wbdie(...
          'invalid usage (failed to determine field name)'); end
        vx(l,:)={v,x}; i=i+1; l=l+1;
     elseif i<nargs
        vx(l,:)={x, varargin{i+1}}; i=i+2; l=l+1;
     else
        wbdie('invalid usage (name/value mismatch; arg #%d)',i+1);
     end
     if isempty(regexp(vx{end,1},'^[a-zA-Z][\w_]*$'))
        wbdie('invalid variable name ''%s''',vx{end,1});
     end
  end

  n=numel(ah);

  for i=1:size(vx,1), v=vx{i,1}; x=vx{i,2};
     for j=1:n, h=ah(j);
        u=get(h,'UserData');
        if isstruct(u), u=setfield(u,v,x); else u=struct(v,{x}); end
        set(h,'UserData',u);
     end
  end

  if ~nargout, clear u; end

end

