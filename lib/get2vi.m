function get2vi(hg,varargin)
% function get2vi(hg [,'-s'])
%
%     List options of given object handle(s) in vi window
%
% Options
%
%     '-s'   expand handle to struct
%
% see also cmd2vi, help2vi
% Wb,Aug04,16

  if nargin<1 || ~all(ishandle(hg(:)))
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     sflag=getopt('-s');
  getopt('check_error');

  ftmp='/tmp/get2vi.tmp';
  line=repmat('=',1,80); tt={'Type','Tag'}; nt=numel(tt);

  fid=fopen(ftmp,'w');
     for i=1:numel(hg), h=hg(i); s=cell(1,nt);
        for j=1:nt, s{j}=get(h,tt{j}); end
        if ~isequal(s{1},s{2}) && ~isempty(s{2})
             s{2}=[' (tag: ' s{2} ')'];
        else s{2}=''; end
        if nargin>1
             fprintf(fid,'\n%s\n%% handle %g/%g: %s\n',line,i,nargin,[s{:}]);
        else fprintf(fid,'\n%s\n%% %s\n',line,[s{:}]);
        end
        if ~sflag
           s=evalc('get(h)'); fprintf(1,'\n');
        else
           if ishghandle(h)
                q='% NB! showing full structure including hidden elements';
           else q=''; end
           s=evalc('x___ = struct(h)');
           s=regexprep(s,'.*x___[ =]+',[q 10]);
        end
        fprintf(fid,'%s',s);
     end
     fprintf(fid,'%s\n\n',line);
  fclose(fid);

  system(['unset LD_LIBRARY_PATH ; ' 10 ...
     sprintf('vim -R -c ''set ignorecase'' %s',ftmp)]);
  system(['/bin/rm ' ftmp]);

end

