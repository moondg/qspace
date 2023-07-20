function wberr(varargin)
% function wberr(varargin)
% deprecated; replaced by wbdie.m
% Wb,May18,07 / Wb,Jun01,19

   if nargin, s=regexprep(varargin{1},'\\N',char(11));
      if nargin>1
         s=sprintf(s,varargin{2:end});
         s=regexprep(s,'\\N',char(11));
      end
   else s=''; end

   s=regexprep([10 s],'\n','\n   ERR ');
   s=regexprep(s,'\\n','\n   ERR ');
   s=regexprep(s,char(11),[char(10) '   ']);

  q=beep;
  if isequal(q,'off'), beep on; beep; beep off; end

   if isdesktop>1, e=char(27);
        wesc_={ [e '[31m'], [e '[0m'] };
   else wesc_={'',''}; end

   fprintf(1,[wesc_{1}, s wesc_{2} '\n\n']);

   S=dbstack;
   if numel(S)>1, dispstack(S(2:end)); end

  S=struct('message','','identifier','Wb:ERR', 'stack',S(min(2,end)));
  error(S);

end

