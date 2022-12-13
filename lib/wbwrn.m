function wbwrn(varargin)
% function wbwrn(varargin)
% Wb,Jun01,19

% adapted from wberr.m // Wb,Jun01,19

   if nargin>1,   s=sprintf(varargin{:});
   elseif nargin, s=varargin{1};
   else s=''; end

   s=regexprep([10 s],'\n','\n   WRN ');
   s=regexprep(s,'\\n','\n   WRN ');

   if isdesktop>1, e=char(27);
        wesc_={ [e '[35m'], [e '[0m'] };
   else wesc_={ '','' }; end

   fprintf(1,[wesc_{1} s wesc_{2} '\n\n']);

   S=dbstack;
   if numel(S)>1, dispstack(S(2:end)); end

   warning('Wb:WRN','');

end

