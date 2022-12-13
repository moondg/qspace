function ff=dir2(varargin)
% Function: ff=dir2([file regexp])
% 
%    Similar to dir() but also deals with regular expressions
%    for files to to be included or excluded.
%
% Wb,Sep08,17

   getopt('init',varargin);
      ipat=getopt('regexp','');
      xpat=getopt('-x','');
   varargin=getopt('get_remaining'); narg=length(varargin);

   ff=dir(varargin{:});
   nf=numel(ff); mark=ones(1,nf);

   if ~isempty(ipat)
      for i=1:nf
         if isempty(regexp(ff(i).name,ipat)), mark(i)=0; end
      end
   end
   if ~isempty(xpat)
      for i=1:nf
         if ~isempty(regexp(ff(i).name,ipat)), mark(i)=0; end
      end
   end
   ff=ff(find(mark));

end

