function ah=getahs(ftag)
% Function ah=getahs(ftag)
%
%    retrieve axes handles set by smaxis() from specified figure
%    ftag  figure tag (default: gcf)
%
% Wb,Dec28,07

  if nargin>0
     fh=findall(groot,'type','figure','tag',ftag);

     if isempty(fh), ah=[]; return;
     elseif length(fh)>1, wblog('WRN', ...
      ['found several figures with given handle (%g)\n' ...
       'take first.'],length(fh)); fh=fh(1);
     end

  else fh=gcf;
  end

  s=get(fh,'UserData');

  if isempty(s) || ~isstruct(s)
     wblog('WRN','no access to UserData -> acess handles');
     ah=findall(fh,'type','axes');
  else
     ah=s.ah;
  end

  figure(fh);

end

