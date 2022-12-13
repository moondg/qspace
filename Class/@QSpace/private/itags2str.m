function s=itags2str(t)
% function s=itags2str(t)
% Wb,Feb24,16

% outsourced from info.m

 % keep compact by default
 % e.g. info() in compact mode for QSpace array uses avg. width 4
 % per label incl. separator
   if isempty(t), s='';
   elseif iscell(t)
      t=sprintf(',%s',t{:}); s=['{ ' t(2:end) ' }'];
   elseif ischar(t)
      s=['{ ' regexprep(t,'[,;| ]+',',') ' }'];
   else 
      t, error('Wb:ERR','\n   ERR invalid itags');
   end

end

