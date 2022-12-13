function setahs(ahs)
% Function setahs(ahs)
% Wb,Dec28,07

   s=get(gcf,'UserData');

   if isempty(s), s=struct('ah',ahs);
   elseif isstruct(s), s.ah=ahs;
   else disp(s), wberr('invalid figure user data !??'); end

   set(gcf,'UserData',s);

end

