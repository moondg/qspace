function QSdocs(varargin)
% function QSdocs()
% Wb,Oct09,18

  getopt('init',varargin);
     if getopt('-v'), vflag=2;
     elseif getopt('-q'), vflag=0; else vflag=1; end
     tflag=getopt('-t');
  getopt('check_error');

  tex='QSdocs.tex';

  fid=fopen(tex,'w');

  fprintf(fid,'\\documentclass[twoside,12pt]{article}\n');
  fprintf(fid,'\\usepackage{QSdocs}\n\n');
  fprintf(fid,'\\begin{document}\n\n');

  fprintf(fid,['\\noindent ' ...
     'Documentation on mex-functions and QSpace member functions\\\\\n' ...
     'Auto-generated via ' mfilename '.m $\\to$ ' tex ' $\\to$ ' ...
           regexprep(tex,'\.tex','.pdf') '\\\\\n' ...
     'A. Weichselbaum \\copyright{} \\today\n\n' ...
     '\\tableofcontents\n\n\\newpage\n\n' ...
  ]);

  N=59; % estimate lines/page
  n=6;  % estimate lines/subsection
  L=0;

% -------------------------------------------------------------------- %
for m=1:2
  if m==1
      fprintf(fid,'\n\\mysec[MEX routines (bin)]{MEX routines (bin)}\n\n');

      ff=dir('../bin/*.mex*');
      ff=sort_files(ff,'--mtime',180,'-x','mpsFull');
      nf=numel(ff);  % compiled within last 180 days
      wblog(' * ','%g mex files',nf);
  else
      fprintf(fid,'\n\\mysec[MEX routines (util)]{MEX routines (util)}\n\n');

      ff=dir('../util/*.mex*');
      ff=sort_files(ff,'--mtime',180,'-x','helloworld');
      nf=numel(ff);  % compiled within last 180 days
      wblog(' * ','%g mex files',nf);
  end

  for i=1:nf
     if vflag, fprintf(1,'  ... %s\n',ff(i).name); end
     f=regexprep(ff(i).name,'\.mex.*','');
     s=regexprep(f,'_','\\_');
     q=skip_trailing_newlines(evalc([f ' -h']));
     q=skip_matlab_formats(q);
     l=numel(find(q==10));

     check_newpage(fid,N,L,l);

     fprintf(fid,'\n\\mymex{%s} %% l=%g/%g\n\n',s,l,L);
     fprintf(fid,'\\begin{verbatim}\n%s\n\\end{verbatim}\n',q);

     L=mod(L+l+n,N);
  end
end

  L=0; fprintf(fid,'\n\\newpage\n\n');

% -------------------------------------------------------------------- %
  ff=dir('../Class/@QSpace/*.m');
  ff=sort_files(ff);
  nf=numel(ff);

  wblog(' * ','%g QSpace routines',nf);

  fprintf(fid,'\n\\mysec[QSpace member functions]{Class/@QSpace member functions}\n\n');

  for i=1:nf
     if vflag, fprintf(1,'  ... %s\n',ff(i).name); end
     f=regexprep(ff(i).name,'\.mex.*','');
     s=regexprep(f,'_','\\_');
     q=skip_trailing_newlines(evalc(['help QSpace/' f]));
     q=skip_matlab_formats(q);
     l=numel(find(q==10));

     check_newpage(fid,N,L,l);

     fprintf(fid,'\n\\myfun{QSpace/%s} %% l=%g/%g\n\n',s,l,L);
     fprintf(fid,'\\begin{verbatim}\n%s\n\\end{verbatim}\n',q);

     L=mod(L+l+n,N);
  end

  fprintf(fid,'\n\\end{document}\n');
  fclose(fid);

  wblog('-->','pdflatex %s',tex);

  pdflatex='unset LD_LIBRARY_PATH; pdflatex -halt-on-error ';
  for i=1:2 % 2 latex passes
     if system([pdflatex tex ' > /dev/null']), break; end
  end
  if i==2 % got at least successful first latex pass
     system('/bin/rm -v QSdocs.aux QSdocs.log QSdocs.out QSdocs.toc');
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function ff=sort_files(ff,varargin)

  getopt('init',varargin);
     xpat =getopt('-x',''); % Wb,Aug28,22
     ndays=getopt('--mtime',1E9);
  getopt('check_error');

  if ~isempty(xpat)
     ff=ff(cellfun(@isempty,regexp({ff.name},xpat)));
  end

  if nargin>1
     tt=[ff.datenum]; i=find(tt>max(tt)-ndays);
     ff=ff(i);
  end

  q={ff.name};
  if 0
     for i=1:numel(q), q{i}=lower(q{i}); end
  else
     for i=1:numel(q)
        if ~isempty(regexp(q{i},'NRG')), q{i}=['NRG' q{i}]; end
     end
  end

  [q,is]=sort(q); ff=ff(is);

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function check_newpage(fid,N,L,l)

   fprintf(fid,'\n\\pagebreak[3]\n\n');
   return

   if l<20
      if L>40 && mod(L+l,N)<10
         fprintf(fid,'\n\\pagebreak[3]\n\n'); % newpage
      end
   elseif l>50
      if L>50
         fprintf(fid,'\n\\pagebreak[3]\n\n'); % newpage
      end
   end
end

% -------------------------------------------------------------------- %
% skip leading or trailing newlines
% -------------------------------------------------------------------- %

function q=skip_trailing_newlines(q)

   l=length(q); q_=q;
   if l>10
      i1=find(q(1:10)~=10,1);
      i2=find(fliplr(q(end-9:end))~=10,1);
   else
      i1=find(q~=10,1);
      i2=find(fliplr(q)~=10,1);
   end

   if ~isempty(i1) && ~isempty(i2), q=q(i1:end-i2+1); end
end

% -------------------------------------------------------------------- %

function q=skip_matlab_formats(q);

   q=regexprep(q,'</*strong>','');
   q=regexprep(q,'<a[^>]*>',''); % Wb,Jun21,19
   q=regexprep(q,'</a>','');

end

% -------------------------------------------------------------------- %

