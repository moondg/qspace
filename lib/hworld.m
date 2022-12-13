function hworld(varargin)
% function hworld(varargin)
% Wb,Jan25,11

   fprintf(1,'\n   Hello world ...\n\n');
   nsleep=0; i=0;

   while i<nargin, i=i+1;
       if isequal(varargin{i},'--sleep') && i<nargin && regexp(varargin{i+1},'^\d+$')
          nsleep=str2num(varargin{i+1});
          fprintf(1,'  %4d %s %s\n',i,varargin{i},varargin{i+1});
          i=i+1;
       else
          fprintf(1,'  %4d %s\n',i,varargin{i});
       end
   end
   if nargin, fprintf(1,'\n'); end

   fprintf(1,'   ismcc: %g\n',ismcc);
   fprintf(1,'   isdeployed: %g\n',isdeployed);
   fprintf(1,'   usejava(awt): %g\n',usejava('awt'));
   fprintf(1,'   usejava(jvm): %g\n',usejava('jvm'));
   fprintf(1,'   isbatch: %g\n',isbatch);

   if nsleep
      fprintf(1,'\n');
      wblog(' * ','Sleeping for %g seconds ... ',nsleep);
      pause(nsleep);
      wblog(' * ','finished.');
      fprintf(1,'\n');
   end

   if exist('helloworld')==3
      helloworld(varargin{:});
   end

   if isbatch || isdeployed, save2tmp; end

end

