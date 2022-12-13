function s = str2tok(str,varargin)
% s = str2tok(str [,delim])
% returns cell array of string tokens of given input string
% Wb,May10,07

   k=0;
   while 1, k=k+1;
      [s{k},str]=strtok(str,varargin{:});
      if isempty(str), break; end
   end

end

