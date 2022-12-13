
% display host/jobid/mat as footer in lower-right corner of plot
% Wb,Sep29,17

% save/restore current value of s
  if isvar('s'), s_old=s; else s_old=[]; end

  if isvar('mat'), s='';
     if isvar('Imain')
        if isfield(Imain,'host'), s=Imain.host;
        else s=getfield2(Imain,'info','host','--def',''); end
     end
     if ~isempty(s)
        if iscell(s), s=s{1}; end
        if isfield(Imain,'jobid') && ~isempty(Imain.jobid) && iscell(Imain.jobid)
           s={s,Imain.jobid{:}};
           if regexp(s{2},'^\d')
                s(2,1:end-1)={'.'}; s{2,1}=' ('; s{2,end}=')';
           else s(2,1:end-1)={'.'}; s=s(:,2:end);
           end
           s=[s{:}];
        elseif isfield(Imain,'pwd')
           s=[s ' / ' regexprep(Imain.pwd,'.*\/','')];
        end
        s=[s ' :: '];
     else s=''; end

     if iscell(mat)
          s={s sprintf('%s, ',mat{:})};, s=[s{1} untex(s{2}(1:end-2))];
     else s=[s untex(mat)];
     end

     if isset('MAT'), s=[untex(regexprep(MAT,'\.mat','')) 10 s]; end
     header('SE',s,{'FontSize',0.66*get(gca,'FontSize')});
  end

if isvar('HAM')
  s=regexprep(HAM.info.istr,'\n',' ');
  if regexp(s,'with end-spins')
     s=regexprep(s,'(with end-spins)','{\\bf\\color{red}$1}');
     if isset('perBC'), perBC=2; end
  end
  if isset('perBC')
     if perBC>1
          s=[s ' {\bf\color{red} using center geometry}'];
     else s=[s ' {\bf\color{red} using perBC}']; end
  end

  if regexp(s,'Heisenberg ladder.*having SU3')
     if ~isempty(regexp(s,' and qproj=\[2 0\]')) || ...
        size(HAM.ops,1)==1 && isequal(HAM.ops(1).op.Q{1},[2 0])
        s=regexprep(s,'and qproj=\[2 0\]','');
        s=regexprep(s,'\(J=\[1, 1\]\)',...
        '[J_{||}=1, J_\\perp=-\\infty \\Rightarrow q_{rung}=(20)]');

        s=regexprep(s,' with qloc=\[1 0\]','');
     else
        s=regexprep(s,'qloc=\[1 0\]','q_{loc}=(10)');
     end
  else s=regexprep(s,'Heisenberg Hamilton\w*','Heisenberg chain');
  end

  istr2=regexprep(s,'having SU3 spin-operator \(\[1 1\]\)','having S=(11)');

  s=getsym(HAM.info.IS.E);
  if ~isempty(s)
     istr2=[istr2, sprintf(' [%s]',s)]; 
  end

end

  if ~isempty(s_old), s=s_old; clear s_old; else clear s; end

